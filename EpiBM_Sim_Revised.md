# EpiBM-Sim: revised research methodology

## Scope and intended use

EpiBM-Sim is a reproducible, stochastic, spatial, hybrid **toy model** of healthy marrow cells (H), treatment-sensitive leukemia cells (Ls), treatment-resistant leukemia cells (Lr), and a local stromal-support field. It is dimensionless, illustrative, and **not clinically calibrated or validated**. It must not be used to guide patient care, infer clinical efficacy, or set medical thresholds. The purpose of this computational study is to exercise and transparently evaluate a bounded, state-feedback treatment scheduler in a controlled simulation.

## Model and time convention

Agents move on a two-dimensional unit torus with periodic boundaries and interact with a discretized support field. H and leukemia agents have separate birth, death, and motility parameters. The support field recovers toward one and is locally impaired by leukemia density; local support modifies birth probabilities. Cells undergo stochastic movement, local repulsion, death, division subject to spacing and a population cap, and sensitive-to-resistant mutation at daughter division. Treatment adds type-specific death hazards; the resistant treatment hazard is reduced by a resistance factor.

The supplied parameterization uses horizon `T=32`, `dt=0.1`, a `32 x 32` support grid, initial H=36 and leukemia=18, and a half-open treatment window `[18,24)`. Time steps advance intervals `[t[k-1], t[k])`, so treatment is applied when the interval start is in `[18,24)`. In the default fixed arm this corresponds to 60 intervals and dimensionless cumulative exposure 6.0. Baseline treatment death rates are 0.18 for H and 0.75 for leukemia; the resistant factor is 0.10. The implementation itself is the authoritative specification for all parameters and update ordering.

## Scheduler and comparison arms

The fixed historical arm applies treatment for the entire prespecified window. The feedback-scheduled arm observes only the pre-interval healthy count and scheduler state. During the same window it pauses treatment if observed H is below the configured fraction of initial H, observes a configurable minimum off-time after a pause, and enforces a cumulative active-time budget of 6.0. An action is binary per interval. `healthy_floor_fraction` and `treatment_rate_H` are experimental dimensionless knobs, not biological/clinical thresholds or doses. A no-treatment arm is available in the model but is not part of the three experiment batches in this package.

## Experiments and analysis

All experiments use the supplied implementation, NumPy random generator, fixed seeds, paired arms where applicable, and assertions in `run_experiments.py`:

1. **Threshold stress:** five seeds (3101–3105), floors 0.75, 0.95, and 1.00, and fixed versus feedback arms (30 runs).
2. **Toxicity stress:** five seeds, healthy-cell treatment death-rate multipliers/values 0.5, 1.0, and 2.0 with floor 0.75, and fixed versus feedback arms (30 runs).
3. **Grid resolution:** three seeds (3101–3103), support grids 16, 32, and 48, fixed arm (9 runs). Other parameters are held fixed; there is no physical-unit rescaling or calibration.

Reported summaries are means and sample standard deviations (`ddof=1`) across seeds. CSV run-level files retain individual outcomes and paired differences. Tests verify the scheduler guard truth table, prior-state-only observation, half-open treatment window, exposure accounting/budget, deterministic replay, model invariants, and fixed-arm regression against the supplied v6 implementation.

## Results and interpretation

The full rerun completed 30 threshold-stress, 30 toxicity-stress, and 9 grid-resolution simulations, with all model and scheduler assertions passing. The existing seven plain-Python unit tests also pass.

In the threshold stress at the baseline healthy-cell toxicity, the watchdog activated zero times at every tested floor. Fixed and feedback outcomes consequently matched exactly at the tested seeds; this is a valid null result, not evidence that feedback scheduling is universally equivalent. In toxicity stress, feedback began to alter behavior at H treatment rate 1.0 (mean 19.6 watchdog activations, exposure 3.44, versus fixed exposure 6.0); at rate 2.0, feedback had mean 35.6 activations and exposure 1.74 versus fixed exposure 6.0. These small-seed simulation outcomes do not establish benefit or clinical safety. Grid summaries were not monotonic: mean final total leukemia was 146.0, 193.7, and 131.7 for grid sizes 16, 32, and 48, respectively (n=3 each), illustrating substantial stochastic/discretization variability and no demonstrated convergence.

## Reproducibility and files

Run from this directory with Python 3.11+ and dependencies `numpy`, `scipy`, and `matplotlib`:

```bash
python -m py_compile run_experiments.py
python run_experiments.py
python leukemia_hybrid_resistance_v7_os_tests.py
```

The experiment runner overwrites its run CSVs and manifest deterministically for the fixed seeds. Outputs include `stress_runs.csv` / `stress_summary.csv`, `toxicity_stress_runs.csv` / `toxicity_stress_summary.csv`, `grid_resolution_runs.csv` / `grid_resolution_summary.csv`, and `experiment_manifest.json`. Environment, seeds, run counts, and checks are recorded in the manifest. No external or patient data are used.

## Limitations

This is a stylized, dimensionless stochastic model with a small number of seeds, simplified cell and microenvironment dynamics, and intentionally uncalibrated parameters. The grid experiment is a sensitivity check, not a formal convergence study. The tested thresholds are design values, not clinical safety limits. Statistical summaries are descriptive and no inferential or clinical efficacy claim is made.
