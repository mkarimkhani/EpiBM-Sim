"""Plain-Python tests for v7 OS-inspired scheduler; python leukemia_hybrid_resistance_v7_os_tests.py"""
import unittest, sys
from dataclasses import replace
import numpy as np
sys.path.insert(0,'/mnt/data')
import leukemia_hybrid_resistance_v7_os as v7
import leukemia_hybrid_resistance_v6 as v6

class V7Tests(unittest.TestCase):
 def small(self,**kw):
  d=dict(T=2.,dt=.1,grid_n=8,initial_H=12,initial_L=8,max_cells=120,seed=77,treatment_start=.6,treatment_end=1.2)
  d.update(kw); return v7.Params(**d)
 def test_half_open_boundary_and_action(self):
  p=self.small(T=25.,dt=.1,treatment_start=18.,treatment_end=24.)
  r=v7.run_arm(p,'fixed_historical'); acts=[x['action'] for x in r.scheduler_records]
  self.assertEqual([i for i,a in enumerate(acts,1) if a],list(range(181,241)))
  self.assertEqual(acts[180],1); self.assertEqual(acts[240],0)
  with self.assertRaises(ValueError): v7.HybridMarrow(p).step(1,treatment_action=1)
 def test_budget_hard_limit(self):
  p=self.small(); s=v7.FeedbackScheduler(p,v7.SchedulerConfig(exposure_budget=.2,healthy_floor_fraction=0,min_off_time=0))
  for t in (.6,.7,.8,.9):
   a,_,_=s.decide(t,12); s.consume(a)
  self.assertAlmostEqual(s.remaining,0); self.assertEqual(s.decide(1.,12)[0],0)
 def test_watchdog_and_hysteresis(self):
  p=self.small(); s=v7.FeedbackScheduler(p,v7.SchedulerConfig(exposure_budget=1,healthy_floor_fraction=.75,min_off_time=.3))
  self.assertEqual(s.decide(.6,12)[0],1)
  a,g,why=s.decide(.7,8); self.assertEqual((a,g,why),(0,True,'healthy_watchdog'))
  self.assertEqual(s.decide(.8,12)[0],0); self.assertEqual(s.decide(.9,12)[0],0)
  self.assertEqual(s.decide(1.,12)[0],1); self.assertEqual(s.guard_triggers,1)
 def test_deterministic_replay(self):
  p=self.small(seed=123); a=v7.run_feedback(p); b=v7.run_feedback(p)
  self.assertTrue(np.array_equal(a.total_hist,b.total_hist)); self.assertTrue(np.array_equal(a.pos,b.pos))
  self.assertEqual(a.scheduler_records,b.scheduler_records)
 def test_fixed_regression_v6(self):
  p=self.small(seed=231)
  a=v6.HybridMarrow(v6.Params(**p.__dict__)).run(snapshot_times=(0,))
  b=v7.HybridMarrow(p).run(snapshot_times=(0,))
  for f in ('H_hist','Ls_hist','Lr_hist','L_hist','total_hist'): self.assertTrue(np.array_equal(getattr(a,f),getattr(b,f)),f)
  self.assertTrue(np.array_equal(a.pos,b.pos)); self.assertTrue(np.array_equal(a.kind,b.kind)); self.assertTrue(np.array_equal(a.S,b.S))
 def test_seed_normalization_and_matching(self):
  seeds=v7.normalize_seeds(['3101','3102',3103])
  self.assertEqual(seeds,[3101,3102,3103])
  paired={seed:{'fixed_historical':seed,'no_treatment':seed,'feedback_scheduled':seed} for seed in seeds}
  self.assertEqual(paired[3102]['feedback_scheduled'],3102)
  for bad in (['1','01'], ['1.5'], [True], [-1], [2**32]):
   with self.subTest(bad=bad):
    with self.assertRaises(ValueError): v7.normalize_seeds(bad)

 def test_invariants_and_no_treatment(self):
  p=self.small(); r=v7.run_arm(p,'feedback_scheduled'); v7.validate_run(r)
  self.assertTrue(np.all(r.total_hist==r.H_hist+r.L_hist))
  self.assertLessEqual(r.scheduler.remaining,6.0+1e-10)
  nt=v7.run_arm(p,'no_treatment'); self.assertEqual(v7.metrics(nt,nt.scheduler_records)['cumulative_exposure'],0)

if __name__=='__main__': unittest.main(verbosity=2)
