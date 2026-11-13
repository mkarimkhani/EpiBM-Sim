---
title: 'EpiBM-Sim: A Spatial Hybrid Toy Model for Reproducible Feedback-Scheduling Experiments'
tags:
  - Python
  - Agent-based modeling
  - Bone marrow niche
  - Feedback control
  - Software verification
authors:
  - name: Mohsen Karimkhani
    orcid: 0009-0002-3454-4375
    corresponding: true
    email: mohsen.karimkhani@gmail.com
    affiliation: 1
affiliations:
  - name: University of Qom, Iran
    index: 1
date: 13 November 2026
bibliography: paper.bib
repository: https://github.com/mkarimkhani/EpiBM-Sim
---

## Summary

EpiBM-Sim is a Python research prototype for exploring treatment-scheduling rules in a simplified model of the leukemia–bone marrow niche. It represents healthy cells, treatment-sensitive leukemia cells, and treatment-resistant leukemia cells as spatial agents coupled to a discretized local support field. A feedback scheduler can interrupt treatment when a simulated healthy-cell safeguard is triggered, while enforcing a treatment window, minimum off-time, and cumulative-exposure limit. The software also provides fixed-schedule comparisons, seeded experiment batches, run-level outputs, and automated checks.

EpiBM-Sim is intentionally dimensionless and is **not calibrated to patients, treatment regimens, or measured marrow geometry**. Its purpose is to make the behavior and implementation of a bounded scheduler inspectable, not to predict leukemia outcomes or recommend clinical treatment. The source code and replication scripts are hosted publicly in the GitHub repository (`https://github.com/mkarimkhani/EpiBM-Sim`) [@EpiBMSim2026].

## Statement of need

Interactions between leukemic cells and the bone marrow microenvironment are biologically consequential. Experimental work has shown that stromal signaling can affect drug resistance in acute myeloid leukemia, providing a reason to study cell–environment interactions rather than treating cellular populations as isolated compartments [@Palani2024; @Urs2024]. That finding motivates the *question* explored here; it does not validate the particular interactions or parameter values used by EpiBM-Sim.

A separate software problem arises when treatment policies are evaluated in stochastic simulations. An apparently favorable result may depend on when observations are taken, whether a controller has inadvertently used future information, how treatment exposure is counted, or whether different policies were compared under matching random seeds [@Sandve2013; @Wilson2017]. These issues are difficult to assess from final population counts alone. EpiBM-Sim addresses them by exposing the scheduler’s timing rules and constraints to automated tests, and by recording per-run outcomes alongside experiment settings.

The intended users are researchers and students investigating computational modeling, agent-based systems, and feedback-control implementation [@terHoeven2025]. The package offers a deliberately small setting in which to examine **what a scheduling rule does under specified assumptions**. It does not claim that those assumptions reproduce a particular leukemia subtype, patient population, or clinical intervention.

## State of the field

Established platforms address broader modeling needs. PhysiCell provides a physics-based framework for three-dimensional multicellular simulation [@Ghaffarizadeh2018]. PhysiBoSS extends a multicellular modeling approach with intracellular signaling models [@Letort2019]. Mesa 3 supplies general-purpose Python infrastructure for building and analyzing agent-based models across research domains [@terHoeven2025]. In parallel, cellular automata and continuous-mathematics models have explored cell-population dynamics within shared niches [@Saikawa2025; @Zhu2024]. These tools are not interchangeable with EpiBM-Sim, nor does EpiBM-Sim seek to replace their capabilities.

EpiBM-Sim instead uses a compact, application-specific model to examine a narrower software question: whether a feedback treatment scheduler behaves as specified when embedded in a stochastic, spatially interacting system. Its simpler architecture makes the comparison arms, safeguards, and recorded experiments relatively easy to inspect. The trade-off is substantial: it does not provide the detailed mechanics, signaling networks, or biological calibration needed for mechanistic or translational claims.

Transparent descriptions of agent-based models are important because results depend not only on stated parameters but also on the order of events and interactions. The ODD reporting framework emphasizes describing those design choices clearly enough for readers to understand and reproduce a model [@Grimm2020]. EpiBM-Sim’s documentation and tests should be read in that spirit: as an account of a specific computational construction, not as evidence that the construction is biologically established.

## Software design

The model combines individually represented cell agents with a grid-based local support field on a periodic unit square. Its three agent populations are healthy cells ($H$), treatment-sensitive leukemia cells ($L_s$), and treatment-resistant leukemia cells ($L_r$). Spatial representation allows local environmental variation and stochastic differences between runs to affect simulated population trajectories. The support field is a modeling abstraction; it should not be interpreted as a measured stromal quantity.

The core implementation is provided in `leukemia_hybrid_resistance_v7_os.py`, and complete replication experiments can be executed via `run_experiments.py` from the project repository (`https://github.com/mkarimkhani/EpiBM-Sim`) [@EpiBMSim2026]. The scheduler is separated conceptually from the population-update process. For each simulation interval, its treatment decision uses the healthy-cell state available **before** that interval is advanced. This ordering matters: deciding from the state produced after treatment would give the controller information unavailable at decision time. The feedback policy includes a configurable healthy-cell floor, a minimum off-time after interruption, a half-open treatment window $[18,24)$, and a hard cumulative-exposure budget of 6.0 simulation-time units under the supplied settings. The fixed-schedule arm provides a comparison within the same treatment window; it is not an independently optimized clinical standard.

This separation makes individual safeguards testable. The supplied test module (`leukemia_hybrid_resistance_v7_os_tests.py`) contains seven test methods covering, among other properties, treatment-window boundaries, the exposure limit, watchdog behavior, deterministic replay, and a fixed-treatment regression comparison with the supplied v6 implementation. The experiment runner additionally records checks for exposure accounting and use of pre-interval observations. These tests can reveal implementation errors within the tested cases, but passing them cannot establish that the biological model is valid.

The runner organizes three experiment batches. A healthy-floor stress experiment compares fixed and feedback arms at floor fractions of 0.75, 0.95, and 1.00 using five seeds per setting. A toxicity stress experiment compares the same arms at **absolute healthy-cell treatment-rate settings** of 0.5, 1.0, and 2.0, again using five seeds per setting; these numbers replace the model parameter and are **not multipliers** of its default value. A grid-resolution experiment evaluates the fixed arm on grids of 16, 32, and 48 using three seeds per grid. Together, these batches produce 69 recorded runs. Their settings and outputs appear in the experiment manifest, run-level CSV files, and summary CSV files.

Keeping seeds, settings, and per-run results identifiable makes the reported comparisons easier to audit and repeat, following standard recommendations for computational research robustness [@Taschuk2017; @Stodden2016; @Peng2011]. It does not remove uncertainty arising from limited seeds or model structure.

## Research impact statement

EpiBM-Sim’s presently demonstrable contribution is a reproducible example of **scheduler verification**, rather than documented adoption by other research groups. Its experiments show both when a safeguard remains inactive and what changes when it activates. In the healthy-floor stress runs, the feedback safeguard did not activate at any of the three tested floor settings; fixed and feedback runs had the same mean exposure and outcomes under their matched settings. This is an informative negative result: a policy should not be credited with an effect that these runs do not show.

The toxicity runs demonstrate a different behavior. At a healthy-cell treatment rate of 0.5, the safeguard did not activate. At rates of 1.0 and 2.0, mean feedback exposure fell to **3.44** and **1.74** simulation-time units, respectively, versus **6.0** for the corresponding fixed schedules. The feedback arms also had higher mean final healthy-cell counts. Crucially, this came with **higher**, not lower, mean final leukemia counts: **827.6 versus 415.6** at rate 1.0, and **1079 versus 458** at rate 2.0, comparing feedback with fixed scheduling. Reporting both populations is necessary to avoid presenting reduced treatment exposure as an unqualified improvement. These are outcomes of the supplied toy-model runs, not estimates of therapeutic benefit.

The grid experiment yielded mean final leukemia counts of **146.0**, **193.7**, and **131.7** for grids of 16, 32, and 48, respectively. This non-monotonic pattern, obtained with three seeds per grid and without physical rescaling, is **not** a convergence result. It instead identifies grid dependence as a limitation requiring further study.

The current materials do not document external users, a publication enabled by the software, or integration into an established workflow. A plausible research use is to adapt the explicit scheduling checks to other computational experiments; that potential has not yet been demonstrated through independent use. Reusable documentation, installation, testing, and release practices would strengthen that case.

## Limitations and responsible use

All model quantities are dimensionless, and the supplied parameters have not been fitted to clinical observations. The small experiment batches do not support population-level inference, robust optimization, or claims of numerical convergence. The grid experiment does not establish resolution independence. Software checks establish selected implementation properties only; they cannot establish biological realism. **EpiBM-Sim must not be used for diagnosis, prognosis, or treatment decisions.**

## AI usage disclosure

OpenAI GPT-6, accessed through GapGPT, was used to assist with drafting and editing the manuscript text and checking bibliographic information. The author reviewed the AI-assisted text and remains responsible for the accuracy and integrity of the submitted work.
## Acknowledgements

*To be completed by the authors: acknowledge actual contributors and funding, if any. Do not retain this placeholder in the submitted paper.*

