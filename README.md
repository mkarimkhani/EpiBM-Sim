# EpiBM-Sim submission package

A reproducible, dimensionless stochastic spatial hybrid **toy simulation** of healthy marrow and leukemia populations with fixed-window and feedback-scheduled treatment arms. It is for software/methodology research only—not clinically calibrated or validated, and not for patient or treatment decisions.

## Quick start

Requirements: Python 3.11+, `numpy`, `scipy`, `matplotlib`.

```bash
python -m py_compile run_experiments.py
python run_experiments.py
python leukemia_hybrid_resistance_v7_os_tests.py
```

The runner performs the complete threshold stress, toxicity stress, and grid-resolution suite. Current verified run: 30 + 30 + 9 = 69 experiments; all embedded assertions passed. The seven supplied plain-Python tests passed.

## Contents

- `leukemia_hybrid_resistance_v7_os.py` — model, scheduler, and validation utilities.
- `run_experiments.py` — corrected, reproducible experiment runner.
- `leukemia_hybrid_resistance_v7_os_tests.py` — unit tests (imports the supplied v6 comparator from `/mnt/data` when run in this workspace; provide that comparator on `PYTHONPATH` elsewhere).
- `EpiBM_Sim_Revised.md` — concise research methodology and results.
- `paper.md` — JOSS-style software paper draft.
- `experiment_manifest.json` — seeds, settings, environment, counts, and checks.
- `stress_runs.csv`, `stress_summary.csv` — floor/threshold stress results.
- `toxicity_stress_runs.csv`, `toxicity_stress_summary.csv` — toxicity stress results.
- `grid_resolution_runs.csv`, `grid_resolution_summary.csv` — grid sensitivity results.

## Experiment design and important results

Five seeds (3101–3105) are used for threshold floors 0.75, 0.95, 1.00 and toxicity values 0.5, 1.0, 2.0; three seeds (3101–3103) are used for grids 16, 32, 48. Outcomes are dimensionless model quantities. Threshold stress activated no watchdogs at baseline toxicity; toxicity stress did activate the guard at rates 1.0 and 2.0. The grid outcomes are non-monotonic and do not establish convergence. See the methodology for interpretation and limitations.

## Clinical-use warning

This is an uncalibrated illustrative model. Its thresholds are not clinical cutoffs, and results must not inform patient care, dosing, safety, or efficacy claims.
