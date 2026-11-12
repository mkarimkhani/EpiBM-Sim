# EpiBM-Sim: a reproducible spatial hybrid toy model for feedback-scheduling experiments

## Summary

EpiBM-Sim is a small Python research-software package for reproducible simulation experiments involving healthy marrow cells, sensitive and resistant leukemia populations, a local support field, and a budget-constrained feedback scheduler. Its purpose is methodological: to explore software behavior and scheduler logic in a stochastic, spatial, dimensionless toy system. It is not a clinical model, has not been clinically calibrated or validated, and is not intended for patient-specific inference or treatment decisions.

The model represents cells as agents on a periodic unit square and the stromal support as a discretized field. Stochastic movement, repulsion, births, deaths, limited mutation, and treatment-dependent death hazards couple cell populations to their local environment. The scheduler uses only the healthy-cell count from the state available before the current time interval; it can pause treatment below a configurable healthy-cell floor and imposes both a minimum off-time and a hard cumulative-exposure budget. The treatment interval is explicitly half-open, `[18,24)`, at `dt=0.1`.

## Statement of need

Small simulation projects are useful for testing whether feedback-control logic respects timing, observability, and resource constraints before considering more elaborate models. However, stochastic agent-based models can obscure update ordering, leak future state into decisions, or accidentally exceed a specified action budget. EpiBM-Sim provides an inspectable implementation and assertion-based experiments that make these properties explicit. Its contribution is not a clinically predictive leukemia model; it is a reproducible software example for evaluating scheduler safeguards and sensitivity to selected toy-model settings.

## Model and software

The core implementation is `leukemia_hybrid_resistance_v7_os.py`. The agent populations are H (healthy), Ls (treatment-sensitive leukemia), and Lr (treatment-resistant leukemia); local support is represented by a periodic grid field. Parameters are dimensionless and specified in the implementation. `run_experiments.py` runs three fixed-seed batches: healthy-floor threshold stress, healthy-cell toxicity stress, and fixed-arm support-grid resolution. Run-level and summary CSVs, plus an experiment manifest, are written alongside the scripts. The unit-test script exercises scheduler behavior and model invariants.

## Verification and example results

The revised runner compiles and completes all 69 experiment runs (30 threshold stress, 30 toxicity stress, and 9 grid-resolution runs). Assertions check that guard decisions use the observed pre-interval healthy state, that actions stay within the allowed window, and that cumulative active time does not exceed 6.0. Seven plain-Python unit tests pass, including deterministic replay and regression against the supplied v6 fixed-treatment model.

At baseline toxicity, the watchdog does not activate at floors 0.75, 0.95, or 1.00 in these five-seed threshold runs, so fixed and feedback outcomes match for the tested seeds. In the toxicity stress, the feedback guard activates at rates 1.0 and 2.0, reducing mean exposure to 3.44 and 1.74 time units, respectively, compared with 6.0 for fixed scheduling. Grid-resolution mean final leukemia counts (n=3) are 146.0, 193.7, and 131.7 for grids 16, 32, and 48. These are descriptive outcomes from a limited simulation ensemble and do not imply clinical benefit or convergence.

## Reproducibility

Requirements: Python 3.11 or later, NumPy, SciPy, and Matplotlib. From the package directory:

```bash
python -m py_compile run_experiments.py
python run_experiments.py
python leukemia_hybrid_resistance_v7_os_tests.py
```

Seeds and experiment settings are fixed in `run_experiments.py` and recorded in `experiment_manifest.json`. CSV summaries report means and sample standard deviations across the specified seeds; run-level CSVs preserve per-seed outcomes. No external or patient data are included.

## Limitations and responsible use

The system is dimensionless, simplified, stochastic, and uncalibrated. Small ensembles are not sufficient for broad statistical claims; the grid comparison is not a formal convergence analysis. Threshold and rate values are model parameters, not medical recommendations. The software is **not for clinical use**, patient decisions, or estimating treatment efficacy or safety.

## Acknowledgements

The package revises and reruns the supplied EpiBM v7 OS-inspired model and preserves its fixed-arm regression test against the supplied v6 implementation.
