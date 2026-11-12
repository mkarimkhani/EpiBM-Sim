"""Proof-of-concept, dimensionless stochastic hybrid leukemia–marrow toy model.

WARNING: PROOF OF CONCEPT ONLY. DIMENSIONLESS AND NOT CLINICAL.
Not clinically calibrated, validated, or suitable for clinical/patient decisions.
The agents are illustrative H (healthy), Ls (sensitive leukemia), and Lr
(resistant leukemia), coupled to an illustrative local stromal support field.

Time convention: history[k] is the state at t=k*dt. Step k advances the
half-open interval [t[k-1], t[k]); treatment applies exactly when the interval
start is in [treatment_start, treatment_end). Thus a [18,24) treatment window
uses intervals [18.0,18.1),...,[23.9,24.0) for dt=.1, and records at 18.0
through 23.9 are included in time-mask statistics (the 24.0 snapshot is off).
Mutation-window counters use the interval in which the daughter division occurs.

Dependencies: numpy, matplotlib, scipy. Run: python leukemia_hybrid_resistance_v5.py
"""
from dataclasses import dataclass, replace
from pathlib import Path
import math
import hashlib, json, platform, time, datetime, csv, sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree

WARNING = 'PROOF OF CONCEPT ONLY • DIMENSIONLESS • NOT CLINICAL / NOT FOR PATIENT DECISIONS'

@dataclass
class Params:
    seed: int = 17
    T: float = 32.0
    dt: float = 0.10
    grid_n: int = 32
    initial_H: int = 36
    initial_L: int = 18
    max_cells: int = 5000
    birth_H: float = 0.34
    birth_L: float = 0.48
    death_H: float = 0.10
    death_L: float = 0.12
    move_H: float = 0.0020
    move_L: float = 0.0026
    cell_diameter: float = 0.018
    division_offset: float = 0.021
    support_recovery: float = 0.30
    leukemia_impairment: float = 0.055
    treatment_start: float = 18.0
    treatment_end: float = 24.0
    treatment_rate_H: float = 0.18
    treatment_rate_L: float = 0.75
    initial_resistant_fraction: float = 0.05
    mutation_probability: float = 0.002
    treatment_resistance_factor: float = 0.10

    def validate(self):
        for name in ('T','dt','cell_diameter','division_offset','birth_H','birth_L','death_H','death_L','move_H','move_L','support_recovery','leukemia_impairment','treatment_start','treatment_end','treatment_rate_H','treatment_rate_L','initial_resistant_fraction','mutation_probability','treatment_resistance_factor'):
            value=getattr(self,name)
            if not math.isfinite(float(value)): raise ValueError(name+' must be finite')
        for name in ('grid_n','initial_H','initial_L','max_cells'):
            value=getattr(self,name)
            if isinstance(value,bool) or int(value)!=value: raise ValueError(name+' must be an integer')
        if self.T <= 0 or self.dt <= 0 or self.grid_n < 1 or self.initial_H<0 or self.initial_L<0 or self.max_cells < self.initial_H+self.initial_L:
            raise ValueError('Invalid horizon, population, grid or cap')
        for name in ('birth_H','birth_L','death_H','death_L','move_H','move_L','support_recovery','leukemia_impairment','treatment_rate_H','treatment_rate_L'):
            if getattr(self,name)<0: raise ValueError(name+' must be nonnegative')
        if not math.isclose(self.T/self.dt, round(self.T/self.dt), abs_tol=1e-9):
            raise ValueError('T must be an integer multiple of dt')
        for name in ('initial_resistant_fraction','mutation_probability','treatment_resistance_factor'):
            if not 0 <= getattr(self,name) <= 1: raise ValueError(name+' must lie in [0,1]')
        if not (0 <= self.treatment_start <= self.treatment_end <= self.T): raise ValueError('Invalid treatment window')
        for t in (self.treatment_start,self.treatment_end):
            if not math.isclose(t/self.dt,round(t/self.dt),abs_tol=1e-9): raise ValueError('Treatment boundaries must align to dt')
        if self.cell_diameter <= 0: raise ValueError('cell_diameter must be positive')
        return int(round(self.T/self.dt))

class HybridMarrow:
    """Point-agent toy model; local stromal impairment uses local absolute density."""
    def __init__(self,p):
        self.p=p; self.nsteps=p.validate(); self.rng=np.random.default_rng(p.seed)
        n=p.initial_H+p.initial_L; self.pos=self.rng.random((n,2))
        nr=int(self.rng.binomial(p.initial_L,p.initial_resistant_fraction))
        self.kind=np.r_[np.zeros(p.initial_H,dtype=np.int8),np.ones(p.initial_L-nr,dtype=np.int8),np.full(nr,2,dtype=np.int8)]
        self.S=np.ones((p.grid_n,p.grid_n),float); self.times=np.arange(self.nsteps+1)*p.dt
        self.H_hist=np.zeros(self.nsteps+1,int); self.Ls_hist=np.zeros_like(self.H_hist); self.Lr_hist=np.zeros_like(self.H_hist)
        self.L_hist=np.zeros_like(self.H_hist); self.total_hist=np.zeros_like(self.H_hist)
        self.current_step=0; self.cap_hit=False
        self.attempted_divisions=0; self.accepted_divisions=0; self.spatial_rejected=0; self.cap_rejected=0
        self.mutations_pre=0; self.mutations_during=0; self.mutations_post=0
        self.snapshots={}; self.snapshots_S={}; self._record(0)

    @property
    def mutations(self): return self.mutations_pre+self.mutations_during+self.mutations_post
    def _grid_indices(self,pos): return np.floor((pos%1)*self.p.grid_n).astype(int)%self.p.grid_n
    def _support_at(self,pos):
        ij=self._grid_indices(pos); return self.S[ij[:,0],ij[:,1]]
    def _record(self,step,snapshot=False):
        self.H_hist[step]=np.count_nonzero(self.kind==0); self.Ls_hist[step]=np.count_nonzero(self.kind==1); self.Lr_hist[step]=np.count_nonzero(self.kind==2)
        self.L_hist[step]=self.Ls_hist[step]+self.Lr_hist[step]; self.total_hist[step]=len(self.kind)
        if snapshot: self.snapshots[step]=(self.pos.copy(),self.kind.copy()); self.snapshots_S[step]=self.S.copy()
    def _move_and_repel(self):
        p=self.p; diffusion=np.where(self.kind==0,p.move_H,p.move_L)
        self.pos=(self.pos+self.rng.normal(size=self.pos.shape)*np.sqrt(2*diffusion[:,None]*p.dt))%1
        if len(self.pos)<2: return
        pairs=cKDTree(self.pos,boxsize=1.0).query_pairs(p.cell_diameter,output_type='ndarray')
        if len(pairs)==0: return
        i,j=pairs[:,0],pairs[:,1]; delta=self.pos[i]-self.pos[j]; delta-=np.round(delta)
        dist=np.linalg.norm(delta,axis=1); valid=dist>1e-12
        i,j,delta,dist=i[valid],j[valid],delta[valid],dist[valid]
        corr=(0.25*(p.cell_diameter-dist))[:,None]*delta/dist[:,None]
        shift=np.zeros_like(self.pos); np.add.at(shift,i,corr); np.add.at(shift,j,-corr); self.pos=(self.pos+shift)%1
    def _division_bins(self,positions):
        # Vectorized bin assignment; stable sorting retains original index order
        # within each bin, exactly as the v5 append loop did.
        nb=max(1,int(1/self.p.cell_diameter))
        if len(positions)==0: return {},nb
        ij=np.floor((positions%1)*nb).astype(np.int64)%nb
        code=ij[:,0]*nb+ij[:,1]
        order=np.argsort(code,kind='stable'); sorted_code=code[order]
        starts=np.r_[0,np.flatnonzero(np.diff(sorted_code))+1]
        ends=np.r_[starts[1:],len(order)]
        bins={}
        for a,b in zip(starts,ends):
            c=int(sorted_code[a]); bins[(c//nb,c%nb)]=order[a:b].tolist()
        return bins,nb
    def step(self,step,treatment_action=None):
        p=self.p; interval_start=(step-1)*p.dt
        ij=self._grid_indices(self.pos); Lcount=np.zeros_like(self.S); mask=self.kind!=0
        np.add.at(Lcount,(ij[mask,0],ij[mask,1]),1)
        self.S+=p.dt*(p.support_recovery*(1-self.S)-p.leukemia_impairment*Lcount); np.clip(self.S,0,1,out=self.S)
        self._move_and_repel(); support=self._support_at(self.pos)
        in_treatment=(p.treatment_start<=interval_start<p.treatment_end) if treatment_action is None else bool(treatment_action)
        if in_treatment and not (p.treatment_start<=interval_start<p.treatment_end): raise ValueError("Treatment action outside the original window")
        death=np.where(self.kind==0,p.death_H,p.death_L).astype(float)
        if in_treatment:
            death+=np.where(self.kind==0,p.treatment_rate_H,np.where(self.kind==1,p.treatment_rate_L,p.treatment_rate_L*p.treatment_resistance_factor))
        die=self.rng.random(len(self.kind)) < -np.expm1(-death*p.dt)
        self.pos=self.pos[~die]; self.kind=self.kind[~die]; support=self._support_at(self.pos)
        br=np.where(self.kind==0,p.birth_H*(.45+.55*support),p.birth_L*(.35+.65*support))
        parents=np.flatnonzero(self.rng.random(len(self.kind)) < -np.expm1(-br*p.dt)); self.rng.shuffle(parents)
        self.attempted_divisions+=len(parents); bins,nb=self._division_bins(self.pos); add_pos=[]; add_kind=[]; occupied=[xy.copy() for xy in self.pos]
        for order,i in enumerate(parents):
            if len(self.kind)+len(add_pos)>=p.max_cells:
                rejected=len(parents)-order; self.cap_rejected+=rejected; self.cap_hit=True; break
            theta=self.rng.uniform(0,2*np.pi); candidate=(self.pos[i]+p.division_offset*np.array([np.cos(theta),np.sin(theta)]))%1.
            key=tuple(np.floor(candidate*nb).astype(int)%nb); nearby=[]
            for dx in (-1,0,1):
                for dy in (-1,0,1): nearby.extend(bins.get(((key[0]+dx)%nb,(key[1]+dy)%nb),[]))
            pool=np.asarray([occupied[j] for j in nearby]) if nearby else np.empty((0,2)); d=candidate-pool; d-=np.round(d)
            if len(pool) and np.any(np.sum(d*d,axis=1)<(p.cell_diameter*.88)**2):
                self.spatial_rejected+=1; continue
            parent_kind=self.kind[i]; daughter_kind=parent_kind
            if parent_kind==1 and self.rng.random()<p.mutation_probability:
                daughter_kind=2
                if interval_start<p.treatment_start: self.mutations_pre+=1
                elif interval_start<p.treatment_end: self.mutations_during+=1
                else: self.mutations_post+=1
            add_pos.append(candidate); add_kind.append(daughter_kind); self.accepted_divisions+=1
            occupied.append(candidate); bins.setdefault(key,[]).append(len(occupied)-1)
        if add_pos:
            self.pos=np.vstack((self.pos,np.asarray(add_pos))); self.kind=np.r_[self.kind,np.asarray(add_kind,dtype=np.int8)]
        self.current_step=step; self._record(step,snapshot=(step in self.snapshot_steps))
    def run(self,snapshot_times=(0,18,24,32)):
        self.snapshot_steps=set(int(round(t/self.p.dt)) for t in snapshot_times)
        assert all(np.isclose(self.times[int(round(t/self.p.dt))],t) for t in snapshot_times)
        self._record(0,snapshot=(0 in self.snapshot_steps))
        for s in range(1,self.nsteps+1): self.step(s)
        return self

def correlation_ratio(values, groups):
    """Correlation ratio eta^2 for categorical group labels (0..1)."""
    y=np.asarray(values,dtype=float); g=np.asarray(groups)
    ok=np.isfinite(y)
    y=y[ok]; g=g[ok]
    if len(y)==0 or np.var(y)==0: return 0.0
    _,inv=np.unique(g,return_inverse=True); mean=y.mean()
    between=sum(np.count_nonzero(inv==k)*(y[inv==k].mean()-mean)**2 for k in range(inv.max()+1))
    return float(between/np.sum((y-mean)**2))

def latin_hypercube(n, bounds, seed=0):
    """Deterministic plain-NumPy Latin hypercube on named (lo,hi) bounds."""
    rng=np.random.default_rng(seed); names=list(bounds); out={}
    for name in names:
        lo,hi=bounds[name]; out[name]=lo+(hi-lo)*(rng.permutation(n)+rng.random(n))/n
    return out

def _rankdata(x):
    x=np.asarray(x); order=np.argsort(x,kind='mergesort'); ranks=np.empty(len(x),float)
    i=0
    while i<len(x):
        j=i+1
        while j<len(x) and x[order[j]]==x[order[i]]: j+=1
        ranks[order[i:j]]=(i+j-1)/2+1; i=j
    return ranks

def partial_rank_correlations(X,y):
    """PRCC: rank-transform and correlate residuals after controlling other X."""
    X=np.asarray(X,float); y=np.asarray(y,float); out=[]
    for j in range(X.shape[1]):
        others=[k for k in range(X.shape[1]) if k!=j]
        A=np.column_stack([np.ones(len(y))]+[_rankdata(X[:,k]) for k in others])
        rx=_rankdata(X[:,j]); ry=_rankdata(y)
        ex=rx-A@np.linalg.lstsq(A,rx,rcond=None)[0]; ey=ry-A@np.linalg.lstsq(A,ry,rcond=None)[0]
        denom=np.linalg.norm(ex)*np.linalg.norm(ey)
        out.append(float(ex@ey/denom) if denom else 0.0)
    return np.asarray(out)

def validate_run(m):
    assert np.isfinite(m.pos).all() and np.isfinite(m.S).all()
    assert np.all(m.pos>=0) and np.all(m.pos<1) and np.all(m.S>=0) and np.all(m.S<=1)
    for a in (m.H_hist,m.L_hist,m.Ls_hist,m.Lr_hist,m.total_hist): assert np.all(a>=0)
    assert np.array_equal(m.L_hist,m.Ls_hist+m.Lr_hist)
    assert np.array_equal(m.total_hist,m.H_hist+m.L_hist)
    assert m.total_hist.max()<=m.p.max_cells and len(m.kind)==m.total_hist[-1]
    assert set(np.unique(m.kind)).issubset({0,1,2})
    assert m.attempted_divisions==m.accepted_divisions+m.spatial_rejected+m.cap_rejected
    assert m.mutations==m.mutations_pre+m.mutations_during+m.mutations_post<=m.accepted_divisions
    assert m.snapshots and all(0<=i<=m.nsteps for i in m.snapshots)

def ensemble_stats(runs):
    result={}
    for field in ('H_hist','Ls_hist','Lr_hist','L_hist','total_hist'):
        a=np.stack([getattr(r,field) for r in runs]); result[field]=(a.mean(axis=0),a.std(axis=0,ddof=1))
    return result

def valid_mean_sd(values):
    x=np.asarray(values,float); x=x[np.isfinite(x)]
    return (float(x.mean()),float(x.std(ddof=1)),len(x)) if len(x)>1 else (float(x.mean()) if len(x) else float('nan'),float('nan'),len(x))

def time_mask_summary(run, start, end, field):
    # Statistics use recorded history timestamps and the explicitly half-open [start,end) mask.
    mask=(run.times>=start)&(run.times<end)
    assert mask.any() and run.times[mask][0]>=start and run.times[mask][-1]<end
    return run.times[mask],getattr(run,field)[mask]

def resistant_window_fraction(run,start=None,end=None):
    start=run.p.treatment_start if start is None else start; end=run.p.treatment_end if end is None else end
    _,denom=time_mask_summary(run,start,end,'L_hist'); _,numer=time_mask_summary(run,start,end,'Lr_hist')
    frac=np.full(denom.shape,np.nan); np.divide(numer,denom,out=frac,where=denom>0)
    return float(np.nanmean(frac)) if np.isfinite(frac).any() else float('nan')

def plot_ensemble(base_runs,treated_runs,output):
    p=treated_runs[0].p; bs=ensemble_stats(base_runs); ts=ensemble_stats(treated_runs); times=treated_runs[0].times
    fig=plt.figure(figsize=(15,11),constrained_layout=True); gs=fig.add_gridspec(3,4,height_ratios=(1.0,1.0,1.05))
    metrics=[('H_hist','Healthy H','#2878b5'),('Ls_hist','Sensitive leukemia Ls','#e67e22'),('Lr_hist','Resistant leukemia Lr','#a52a6a'),('total_hist','Total leukemia Ls + Lr','#3a7d44')]
    for i,(field,label,color) in enumerate(metrics):
        ax=fig.add_subplot(gs[i//2,i%2]); ax.axvspan(p.treatment_start,p.treatment_end,color='#8a6bb8',alpha=.12,label='Treatment interval [18, 24)')
        for arm,stats,ls in [('Baseline',bs,'--'),('Treatment',ts,'-')]:
            mean,sd=stats[field]; ci=1.96*sd/np.sqrt(len(base_runs if arm=='Baseline' else treated_runs))
            ax.plot(times,mean,color=color,ls=ls,lw=1.8,label=f'{arm} mean')
            ax.fill_between(times,np.maximum(0,mean-ci),mean+ci,color=color,alpha=.15,linewidth=0,label=f'{arm} 95% CI')
        ax.set(xlabel='Time (dimensionless)',ylabel='Cells',title=label); ax.grid(alpha=.22); ax.set_xlim(0,p.T)
        if i==0: ax.legend(fontsize=8,ncol=2,frameon=False)
    # Representative spatial states are one explicitly identified treatment-arm realization.
    rep=treated_runs[0]; snap_times=(0.,18.,24.,32.); colors={0:'#2878b5',1:'#e67e22',2:'#a52a6a'}
    for k,requested in enumerate(snap_times):
        ax=fig.add_subplot(gs[2,k]); idx=int(round(requested/p.dt)); actual=rep.times[idx]
        assert np.isclose(actual,requested) and idx in rep.snapshots
        pos,kind=rep.snapshots[idx]; field=rep.snapshots_S[idx]
        ax.imshow(field.T,origin='lower',extent=(0,1,0,1),cmap='YlGn',vmin=0,vmax=1,alpha=.52,interpolation='nearest')
        for typ,lab in ((0,'H'),(1,'Ls'),(2,'Lr')):
            q=kind==typ
            if q.any(): ax.scatter(pos[q,0],pos[q,1],s=15,c=colors[typ],edgecolors='white',linewidths=.2,label=lab)
        ax.set(xlim=(0,1),ylim=(0,1),aspect='equal',title=f'Single-run representative state\nt = {actual:g} (dimensionless)',xlabel='x',ylabel='y')
        ax.text(.02,.02,'Treatment active' if p.treatment_start<=actual<p.treatment_end else ('Treatment ended' if actual>=p.treatment_end else 'Before treatment'),transform=ax.transAxes,fontsize=7,color='black',bbox=dict(facecolor='white',alpha=.8,edgecolor='none'))
        if k==0: ax.legend(fontsize=8,loc='upper right')
    fig.suptitle('Hybrid leukemia–marrow toy model | ensemble means and pointwise 95% CIs',fontsize=15,fontweight='bold')
    fig.text(.5,.975,WARNING,ha='center',va='top',fontsize=12,fontweight='bold',color='#9b1c1c')
    fig.text(.5,.005,'Ensemble bands: mean ± 1.96 × SD/√n (pointwise, n runs per arm). Snapshots: first treatment run only; background = dimensionless support S (0–1).',ha='center',fontsize=8)
    fig.savefig(output,dpi=240,bbox_inches='tight'); plt.close(fig)

def run_ensemble(n=10,ensemble_seed=20250308,output=None, sensitivity_reps=3):
    if n<2: raise ValueError('At least two runs per arm are required for SD/CI')
    p=Params(); seedseq=np.random.SeedSequence(ensemble_seed).spawn(2*n+3*sensitivity_reps)
    seeds=[int(s.generate_state(1)[0]) for s in seedseq]
    basep=replace(p,treatment_rate_H=0.,treatment_rate_L=0.)
    base=[]; treated=[]
    for i in range(n):
        b=HybridMarrow(replace(basep,seed=seeds[i])).run(); tr=HybridMarrow(replace(p,seed=seeds[n+i])).run()
        validate_run(b); validate_run(tr); base.append(b); treated.append(tr)
    assert len({r.p.seed for r in base+treated})==2*n
    print(WARNING); print(f'Independent runs per arm: n={n}; ensemble seed={ensemble_seed}; all seeds distinct; dt={p.dt:g}; T={p.T:g}; treatment [18,24)')
    for arm,runs in [('Baseline',base),('Treatment',treated)]:
        print(f'\n{arm} arm:')
        for label,field,at_end in [('Final H','H_hist',True),('Final Ls','Ls_hist',True),('Final Lr','Lr_hist',True),('Final total leukemia','L_hist',True)]:
            vals=[getattr(r,field)[-1] for r in runs]; mean,sd,count=valid_mean_sd(vals); print(f'  {label}: mean={mean:.3f}, SD={sd:.3f}, n={count}')
        vals=[resistant_window_fraction(r) for r in runs]; mean,sd,count=valid_mean_sd(vals); print(f'  Mean resistant fraction over recorded t in [18,24): mean={mean:.4f}, SD={sd:.4f}, n={count}')
        for window_name,lo,hi in [('pre-treatment [0,18)',0,18),('during-treatment [18,24)',18,24)]:
            attr='mutations_pre' if window_name.startswith('pre') else 'mutations_during'
            counts=[getattr(r,attr) for r in runs]; print(f'  De novo mutations {window_name}: total={sum(counts)}, mean/run={np.mean(counts):.3f}, per-run={counts}')
        for attr,label in [('attempted_divisions','attempted'),('accepted_divisions','accepted'),('spatial_rejected','spatial/overlap rejected'),('cap_rejected','cap rejected')]: print(f'  {label} divisions: total={sum(getattr(r,attr) for r in runs)}')
        print(f'  Cap reached: {sum(r.cap_hit for r in runs)}/{n} runs')
    print('\nGrid-resolution sensitivity (3 independent runs/grid; final total leukemia mean ± SD):')
    sens=[]
    for j,grid in enumerate((16,32,48)):
        vals=[]
        for k in range(sensitivity_reps):
            sp=replace(p,grid_n=grid,seed=seeds[2*n+j*sensitivity_reps+k])
            r=HybridMarrow(sp).run(); validate_run(r); vals.append(r.L_hist[-1])
        mean,sd,count=valid_mean_sd(vals); sens.append((grid,mean,sd)); print(f'  grid_n={grid}: {mean:.3f} ± {sd:.3f} (n={count})')
    base_mean=np.mean([x[1] for x in sens if x[0]==32]);
    for grid,mean,sd in sens: print(f'    relative difference from grid_n=32 mean: {100*(mean-base_mean)/max(abs(base_mean),1e-12):+.2f}% (grid {grid})')
    # Verify history/interval alignment explicitly at both boundaries.
    tr=treated[0]; mask=(tr.times>=18)&(tr.times<24)
    assert np.array_equal(np.flatnonzero(mask),np.arange(180,240))
    assert np.isclose(tr.times[np.flatnonzero(mask)[0]],18.) and np.isclose(tr.times[np.flatnonzero(mask)[-1]],23.9)
    assert np.isclose(tr.times[240],24.) and not mask[240]
    for arm,runs in [('Baseline',base),('Treatment',treated)]:
        summary=ensemble_stats(runs)
        print(f'\n{arm} time-averaged history means over [18,24) (history points, all runs):')
        for field in ('H_hist','Ls_hist','Lr_hist','L_hist'):
            vals=np.concatenate([time_mask_summary(r,18,24,field)[1] for r in runs]); m,sd,count=valid_mean_sd(vals)
            print(f'  {field}: mean={m:.3f}, SD={sd:.3f}, pooled points={count}')
    if output: plot_ensemble(base,treated,output)
    return base,treated,sens



# OS-inspired FEEDBACK SCHEDULER: operational policy only, not a biological mechanism.
@dataclass
class SchedulerConfig:
    exposure_budget: float = 6.0  # dimensionless active-time budget, equals historic window duration
    healthy_floor_fraction: float = 0.75  # configurable guard; not a clinical threshold
    min_off_time: float = 0.30  # dimensionless hysteresis interval
    def validate(self):
        if not math.isfinite(self.exposure_budget) or self.exposure_budget < 0: raise ValueError('budget must be finite and nonnegative')
        if not 0 <= self.healthy_floor_fraction <= 1: raise ValueError('healthy floor fraction must be in [0,1]')
        if not math.isfinite(self.min_off_time) or self.min_off_time < 0: raise ValueError('min_off_time must be finite/nonnegative')

class FeedbackScheduler:
    """Healthy-preservation priority, healthy-population watchdog, bounded budget and off-time hysteresis."""
    def __init__(self,p,config=None):
        self.p=p; self.config=config or SchedulerConfig(); self.config.validate()
        self.remaining=float(self.config.exposure_budget); self.active=False; self.off_until=-math.inf
        self.guard_triggers=0; self.hysteresis_blocks=0; self.budget_blocks=0
    def decide(self,t,H):
        # No lookahead: uses only the current recorded state and scheduler state.
        in_window=self.p.treatment_start <= t < self.p.treatment_end
        floor=self.p.initial_H*self.config.healthy_floor_fraction
        guard=in_window and H < floor
        if guard: self.guard_triggers += 1
        if not in_window:
            self.active=False; return 0,guard,'outside_window'
        if guard:
            self.active=False; self.off_until=max(self.off_until,t+self.config.min_off_time)
            return 0,guard,'healthy_watchdog'
        if self.remaining+1e-12 < self.p.dt:
            self.active=False; self.budget_blocks+=1; return 0,guard,'budget_exhausted'
        if self.active:
            return 1,guard,'continue'
        if t < self.off_until-1e-12:
            self.hysteresis_blocks+=1; return 0,guard,'minimum_off_time'
        self.active=True
        return 1,guard,'start_or_resume'
    def consume(self,action):
        if action:
            self.remaining=max(0.,self.remaining-self.p.dt)


def run_feedback(p,config=None):
    model=HybridMarrow(p); sched=FeedbackScheduler(p,config)
    model.snapshot_steps={0}; model._record(0,snapshot=True)
    records=[]
    for step in range(1,model.nsteps+1):
        t=(step-1)*p.dt
        H=int(np.count_nonzero(model.kind==0))
        action,guard,reason=sched.decide(t,H)
        model.step(step,treatment_action=action)
        sched.consume(action)
        records.append(dict(step=step,time=t,healthy_observed=H,leukemia_observed=int(model.L_hist[step-1]),action=action,guard_trigger=int(guard),reason=reason,budget_remaining=sched.remaining))
    model.scheduler=sched; model.scheduler_records=records
    return model

def run_arm(p,arm,config=None):
    if arm=='fixed_historical':
        model=HybridMarrow(p).run(snapshot_times=(0,))
        records=[]
        for s in range(1,model.nsteps+1):
            t=(s-1)*p.dt; a=int(p.treatment_start<=t<p.treatment_end)
            records.append(dict(step=s,time=t,healthy_observed=int(model.H_hist[s-1]),leukemia_observed=int(model.L_hist[s-1]),action=a,guard_trigger=0,reason='fixed_window' if a else 'off',budget_remaining=max(0.,6.-sum(x['action'] for x in records+[{'action':a}])*p.dt)))
        model.scheduler_records=records; model.scheduler=None
        return model
    if arm=='no_treatment':
        q=replace(p,treatment_rate_H=0.,treatment_rate_L=0.)
        model=HybridMarrow(q).run(snapshot_times=(0,)); model.scheduler_records=[]
        for s in range(1,model.nsteps+1):
            t=(s-1)*p.dt; model.scheduler_records.append(dict(step=s,time=t,healthy_observed=int(model.H_hist[s-1]),leukemia_observed=int(model.L_hist[s-1]),action=0,guard_trigger=0,reason='no_treatment',budget_remaining=6.))
        model.scheduler=None; return model
    if arm=='feedback_scheduled': return run_feedback(p,config)
    raise ValueError(arm)

def metrics(model,records):
    h=model.H_hist; total=model.L_hist; sus=model.Ls_hist; res=model.Lr_hist
    actions=np.array([r['action'] for r in records],int)
    return dict(healthy_nadir=int(h.min()),final_healthy=int(h[-1]),final_susceptible=int(sus[-1]),final_resistant=int(res[-1]),final_total_leukemia=int(total[-1]),resistant_fraction=(float(res[-1]/total[-1]) if total[-1] else float('nan')),cumulative_exposure=float(actions.sum()*model.p.dt),guard_trigger_count=int(sum(r['guard_trigger'] for r in records)),treated_intervals=int(actions.sum()))


# Reproducible paired-seed export. Normalize seed values at the boundary: CSV/JSON
# readers often supply strings, while run metadata uses integers.
def normalize_seeds(seed_values):
    seeds=[]
    for value in seed_values:
        if isinstance(value, (bool, np.bool_)):
            raise ValueError('seed must be an integer, not bool')
        try:
            seed=int(value)
        except (TypeError, ValueError, OverflowError):
            raise ValueError(f'invalid integer seed: {value!r}')
        if str(value).strip() not in (str(seed), f'+{seed}') and not isinstance(value, (int, np.integer)):
            raise ValueError(f'seed must be an exact integer: {value!r}')
        if seed < 0 or seed > 2**32-1:
            raise ValueError('seed must be in [0, 2**32-1]')
        seeds.append(seed)
    if not seeds or len(set(seeds)) != len(seeds):
        raise ValueError('seed list must be non-empty and unique')
    return seeds

def _active_windows(records, dt):
    active=[r['time'] for r in records if r['action']]
    if not active: return []
    groups=[]; begin=prev=active[0]
    for t in active[1:]:
        if not math.isclose(t,prev+dt,abs_tol=1e-9):
            groups.append([float(begin),float(prev+dt)]); begin=t
        prev=t
    groups.append([float(begin),float(prev+dt)])
    return groups

def generate_outputs(output_dir='/mnt/data', seed_values=(3101,3102,3103,3104,3105), test_returncode=None, test_output=''):
    """Run paired three-arm simulations and write comparisons, timeseries, plot, manifest."""
    from pathlib import Path
    import hashlib, subprocess
    out=Path(output_dir); out.mkdir(parents=True,exist_ok=True)
    seeds=normalize_seeds(seed_values); p=Params(); arms=('fixed_historical','no_treatment','feedback_scheduled')
    records=[]; timeseries=[]; by_seed={}
    for seed in seeds:
        by_seed[seed]={}
        for arm in arms:
            model=run_arm(replace(p,seed=seed),arm)
            validate_run(model)
            recs=model.scheduler_records; met=metrics(model,recs)
            row=dict(met,seed=seed,arm=arm,active_windows=json.dumps(_active_windows(recs,p.dt),separators=(',',':')))
            by_seed[seed][arm]=row; records.append(row)
            for i,t in enumerate(model.times):
                timeseries.append(dict(seed=seed,arm=arm,time=float(t),healthy=int(model.H_hist[i]),
                    sensitive_leukemia=int(model.Ls_hist[i]),resistant_leukemia=int(model.Lr_hist[i]),
                    total_leukemia=int(model.L_hist[i])))
    # Explicit integer-keyed matching avoids str-vs-int seed lookup failures.
    for seed in seeds:
        fixed=by_seed[seed]['fixed_historical']
        for arm in arms:
            row=by_seed[seed][arm]
            for key in ('healthy_nadir','final_healthy','final_susceptible','final_resistant','final_total_leukemia','resistant_fraction','cumulative_exposure'):
                row['paired_diff_vs_fixed_'+key]=row[key]-fixed[key]
    comp_path=out/'leukemia_hybrid_resistance_v7_os_comparisons.csv'
    ts_path=out/'leukemia_hybrid_resistance_v7_os_timeseries.csv'
    with comp_path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(records[0])); w.writeheader(); w.writerows(records)
    with ts_path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(timeseries[0])); w.writeheader(); w.writerows(timeseries)
    # Compact endpoint panel, summarizes the actual paired simulation outputs.
    fig,axs=plt.subplots(1,3,figsize=(12,4.1),constrained_layout=True)
    labels=['Fixed historical','No treatment','Feedback scheduled']; colors=['#7651a8','#4c9f70','#e18a32']
    fields=[('final_total_leukemia','Final leukemia count'),('final_healthy','Final healthy H count'),('cumulative_exposure','Cumulative exposure')]
    for ax,(field,title) in zip(axs,fields):
        data=[[by_seed[s][a][field] for s in seeds] for a in arms]
        means=[float(np.mean(v)) for v in data]; sds=[float(np.std(v,ddof=1)) if len(v)>1 else 0. for v in data]
        ax.bar(range(3),means,yerr=sds,capsize=4,color=colors,alpha=.86)
        for j,vals in enumerate(data): ax.scatter(np.full(len(vals),j)+np.linspace(-.08,.08,len(vals)),vals,s=18,color='#222',zorder=3)
        ax.set_xticks(range(3),labels,rotation=20,ha='right'); ax.set_title(title); ax.grid(axis='y',alpha=.22)
        ax.set_ylabel('Dimensionless toy-model units')
    fig.suptitle(f'Paired-seed summary (n={len(seeds)} per arm) | proof-of-concept only, not clinical',fontweight='bold')
    fig.savefig(out/'leukemia_hybrid_resistance_v7_os.png',dpi=220,bbox_inches='tight'); plt.close(fig)
    # Test suite is run by the caller; hashes cover code, tests, CSVs and PNG.
    paths={'leukemia_hybrid_resistance_v7_os.py':Path(__file__),
           'leukemia_hybrid_resistance_v7_os_tests.py':out/'leukemia_hybrid_resistance_v7_os_tests.py',
           comp_path.name:comp_path,ts_path.name:ts_path,out.joinpath('leukemia_hybrid_resistance_v7_os.png').name:out/'leukemia_hybrid_resistance_v7_os.png'}
    hashes={name:hashlib.sha256(path.read_bytes()).hexdigest() for name,path in paths.items() if path.exists()}
    summary={}
    for arm in arms:
        summary[arm]={}
        for field in ('healthy_nadir','final_healthy','final_susceptible','final_resistant','final_total_leukemia','resistant_fraction','cumulative_exposure','guard_trigger_count'):
            vals=np.asarray([by_seed[seed][arm][field] for seed in seeds],float)
            summary[arm][field]={'mean':float(vals.mean()),'sd_sample':float(vals.std(ddof=1)) if len(vals)>1 else 0.}
    manifest={'version':'v7_os','disclaimer':WARNING+'; dimensionless proof-of-concept, not clinically calibrated or suitable for decisions.',
      'seeds':seeds,'seed_type':'integer','paired_seed_design':True,'seed_matching':'Rows are matched by normalized integer seed; every seed has exactly one row per arm.',
      'parameters':p.__dict__,'scheduler_config':SchedulerConfig().__dict__,
      'exposure_fairness':'Nominal feedback budget is 6; fixed historical uses 6 and no-treatment 0. Realized exposure differs and comparisons are descriptive only.',
      'runtime_seconds':None,'test_returncode':test_returncode,'test_output':test_output,
      'summary':summary,'sha256':hashes,'output_rows':{'comparisons':len(records),'timeseries':len(timeseries)}}
    (out/'leukemia_hybrid_resistance_v7_os_manifest.json').write_text(json.dumps(manifest,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    return manifest

if __name__=='__main__':
    import argparse
    _parser=argparse.ArgumentParser(description='Generate v7 OS toy-model paired-seed outputs')
    _parser.add_argument('--seeds',nargs='+',default=['3101','3102','3103','3104','3105'])
    _parser.add_argument('--output-dir',default='/mnt/data')
    _args=_parser.parse_args()
    generate_outputs(_args.output_dir,_args.seeds)
