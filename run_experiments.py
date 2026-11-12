#!/usr/bin/env python3
"""Reproducible v7 experiments. Dimensionless toy model; not clinical."""
from dataclasses import replace
import csv, json, hashlib, platform, sys, time
from pathlib import Path
import numpy as np
import leukemia_hybrid_resistance_v7_os as v7
OUT=Path(__file__).resolve().parent
SEEDS=[3101,3102,3103,3104,3105]
FLOORS=[0.75,0.95,1.00]  # preregistered compact stress grid
GRID_SEEDS=[3101,3102,3103]
TOXICITY=[0.5,1.0,2.0]  # H-toxicity stress; rates >=1 provoke watchdog interventions

def write_csv(path, rows):
    if rows:
        with open(path,'w',newline='',encoding='utf-8') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def run_stress():
    rows=[]; base=v7.Params()
    for floor in FLOORS:
        cfg=v7.SchedulerConfig(exposure_budget=6.0,healthy_floor_fraction=floor,min_off_time=.30)
        for seed in SEEDS:
            p=replace(base,seed=seed)
            for arm in ('fixed_historical','feedback_scheduled'):
                m=v7.run_arm(p,arm,cfg if arm=='feedback_scheduled' else None); v7.validate_run(m)
                metrics=v7.metrics(m,m.scheduler_records)
                # guard checks: every trigger iff observed H is below floor while in treatment window
                for r in m.scheduler_records:
                    expected=int(p.treatment_start<=r['time']<p.treatment_end and r['healthy_observed'] < p.initial_H*floor) if arm=='feedback_scheduled' else 0
                    assert r['guard_trigger']==expected
                    # The observation is from the state before this interval (no future leakage).
                    assert r['healthy_observed']==int(m.H_hist[r['step']-1])
                    assert r['time']==(r['step']-1)*p.dt
                assert abs(metrics['cumulative_exposure']-sum(r['action'] for r in m.scheduler_records)*p.dt)<1e-10
                assert metrics['cumulative_exposure']<=6.0+1e-10
                assert metrics['guard_trigger_count']==sum(r['guard_trigger'] for r in m.scheduler_records)
                rows.append(dict(seed=seed,healthy_floor_fraction=floor,arm=arm,**metrics,
                                 guard_activations=sum(r['guard_trigger'] for r in m.scheduler_records),
                                 exposure_accounting_ok=True,no_future_information_check=True,
                                 feedback_changed_vs_fixed=''))
    # explicit paired differences, preserving no-difference rows
    lookup={(r['healthy_floor_fraction'],r['seed'],r['arm']):r for r in rows}
    for r in rows:
        if r['arm']=='feedback_scheduled':
            f=lookup[(r['healthy_floor_fraction'],r['seed'],'fixed_historical')]
            r['feedback_changed_vs_fixed']=any(r[k]!=f[k] for k in ('healthy_nadir','final_healthy','final_susceptible','final_resistant','final_total_leukemia','cumulative_exposure'))
            r['paired_diff_final_leukemia']=r['final_total_leukemia']-f['final_total_leukemia']
            r['paired_diff_exposure']=r['cumulative_exposure']-f['cumulative_exposure']
        else:
            r['feedback_changed_vs_fixed']='NA';r['paired_diff_final_leukemia']=0;r['paired_diff_exposure']=0
    write_csv(OUT/'stress_runs.csv',rows)
    summary=[]
    for floor in FLOORS:
        for arm in ('fixed_historical','feedback_scheduled'):
            a=[r for r in rows if r['healthy_floor_fraction']==floor and r['arm']==arm]
            for metric in ('final_total_leukemia','final_healthy','healthy_nadir','cumulative_exposure','guard_activations'):
                vals=np.array([r[metric] for r in a],float)
                summary.append(dict(healthy_floor_fraction=floor,arm=arm,metric=metric,n=len(vals),mean=float(vals.mean()),sd_sample=float(vals.std(ddof=1))))
    write_csv(OUT/'stress_summary.csv',summary)
    return rows,summary

def run_toxicity_stress():
    rows=[];base=v7.Params()
    for tox in TOXICITY:
      cfg=v7.SchedulerConfig(exposure_budget=6.0,healthy_floor_fraction=.75,min_off_time=.30)
      for seed in SEEDS:
        p=replace(base,seed=seed,treatment_rate_H=tox)
        for arm in ('fixed_historical','feedback_scheduled'):
          m=v7.run_arm(p,arm,cfg if arm=='feedback_scheduled' else None);v7.validate_run(m)
          met=v7.metrics(m,m.scheduler_records)
          for r in m.scheduler_records:
            exp=int(p.treatment_start<=r['time']<p.treatment_end and r['healthy_observed']<p.initial_H*.75) if arm=='feedback_scheduled' else 0
            assert r['guard_trigger']==exp and r['healthy_observed']==int(m.H_hist[r['step']-1])
            assert r['time']==(r['step']-1)*p.dt
          assert abs(met['cumulative_exposure']-sum(r['action'] for r in m.scheduler_records)*p.dt)<1e-10
          assert met['cumulative_exposure']<=6+1e-10
          rows.append(dict(seed=seed,treatment_rate_H=tox,arm=arm,**met,guard_activations=sum(r['guard_trigger'] for r in m.scheduler_records),exposure_accounting_ok=True,no_future_information_check=True))
    lookup={(r['treatment_rate_H'],r['seed'],r['arm']):r for r in rows}
    for r in rows:
      f=lookup[(r['treatment_rate_H'],r['seed'],'fixed_historical')]
      r['paired_diff_final_leukemia']=r['final_total_leukemia']-f['final_total_leukemia']
      r['paired_diff_exposure']=r['cumulative_exposure']-f['cumulative_exposure']
      r['feedback_changed_vs_fixed']= any(r[k]!=f[k] for k in ('healthy_nadir','final_healthy','final_susceptible','final_resistant','final_total_leukemia','cumulative_exposure')) if r['arm']=='feedback_scheduled' else 'NA'
    write_csv(OUT/'toxicity_stress_runs.csv',rows)
    summ=[]
    for tox in TOXICITY:
      for arm in ('fixed_historical','feedback_scheduled'):
        a=[r for r in rows if r['treatment_rate_H']==tox and r['arm']==arm]
        for key in ('final_total_leukemia','final_healthy','healthy_nadir','cumulative_exposure','guard_trigger_count','guard_activations'):
          x=np.array([r[key] for r in a],float); summ.append(dict(treatment_rate_H=tox,arm=arm,metric=key,n=len(x),mean=float(x.mean()),sd_sample=float(x.std(ddof=1))))
    write_csv(OUT/'toxicity_stress_summary.csv',summ)
    return rows,summ

def run_grid():
    rows=[]; base=v7.Params()
    for n in (16,32,48):
      for seed in GRID_SEEDS:
        p=replace(base,grid_n=n,seed=seed)
        m=v7.run_arm(p,'fixed_historical');v7.validate_run(m)
        z=v7.metrics(m,m.scheduler_records)
        rows.append(dict(grid_n=n,seed=seed,arm='fixed_historical',**z,cap_hit=m.cap_hit,attempted_divisions=m.attempted_divisions,accepted_divisions=m.accepted_divisions))
    write_csv(OUT/'grid_resolution_runs.csv',rows)
    summ=[]
    for n in (16,32,48):
      a=[r for r in rows if r['grid_n']==n]
      for metric in ('final_total_leukemia','final_healthy','healthy_nadir','final_susceptible','final_resistant'):
        x=np.array([r[metric] for r in a],float)
        summ.append(dict(grid_n=n,metric=metric,n=len(x),mean=float(x.mean()),sd_sample=float(x.std(ddof=1))))
    write_csv(OUT/'grid_resolution_summary.csv',summ)
    return rows,summ

def main():
    t=time.time();stress,ss=run_stress();tox,ts=run_toxicity_stress();grid,gs=run_grid()
    manifest={'model_source':'leukemia_hybrid_resistance_v7_os.py copied from supplied v7_os implementation',
      'seeds_stress':SEEDS,'stress_floors':FLOORS,'toxicity_grid':TOXICITY,'toxicity_seeds':SEEDS,'grid_seeds':GRID_SEEDS,'grid_n':[16,32,48],
      'physical_unit_limit':'Dimensionless; grid_n alters discretization of support field on unit torus; parameters are not rescaled or calibrated in physical units.',
      'environment':{'python':sys.version,'platform':platform.platform(),'numpy':np.__version__},
      'runtime_seconds':time.time()-t,'stress_rows':len(stress),'toxicity_rows':len(tox),'grid_rows':len(grid),
      'checks':'guard truth table, prior-state-only observation, action window via v7, exposure=sum(action)*dt <= budget, model validate_run; all asserted',
      'input_provenance':'existing supplied v7 source/tests and seeds from v7 manifest; no external or fabricated data'}
    (OUT/'experiment_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'runtime_seconds':manifest['runtime_seconds'],'stress_rows':len(stress),'toxicity_rows':len(tox),'grid_rows':len(grid),'guard_activations_by_floor':{str(f):sum(r['guard_activations'] for r in stress if r['healthy_floor_fraction']==f and r['arm']=='feedback_scheduled') for f in FLOORS},'changed_feedback_rows':sum(r['feedback_changed_vs_fixed'] is True for r in tox if r['arm']=='feedback_scheduled')},indent=2))
if __name__=='__main__':main()
