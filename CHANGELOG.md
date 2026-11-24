# Changelog

## 0.2.1 (24 Nov 2026)

Readability and reproducibility, against the coding guide of the project; results are
unchanged except for the time step (see the first item).

- The start vector of the power iteration that sets the time step comes from
  `StableRNG` (StableRNGs.jl) instead of `MersenneTwister`, whose stream may change
  between Julia versions. The step changes at the level of the power-iteration
  tolerance, so new runs are not bitwise equal to those of v0.2.0. The NRMSE values
  are expected not to change at the digits reported; `Pkg.test()` checks the coarse
  intact case against v0.1 with rtol = 1e-2.
- `IntactPipeCase` takes `n_axial` as a keyword, and `INTACT_SERIES` and the notebooks
  use keywords instead of five positional numbers.
- Named constants for the sensor positions of the notched case and for the observation
  time per metre of the intact case.
- `MODEL_COLORS` (Okabe-Ito) for the figures of the notebooks; the notebooks no longer
  define functions that read global variables.
- `notebooks/examples/README.md` with the index of the examples.
- Run logs moved from the project root to `results/logs/`; the one-off regression
  script against v0.1 (`scripts/compare_with_original.jl`, which needed the v0.1 project
  next to this one) removed from the package, its result is in VALIDATION.md.
- `results/notched_pipe.csv` keeps only the last run of the converged notched case
  (607.5 s); the first run of 22 Nov (1013 s) is in `results/logs/notch_L1_log.txt`.

## 0.2.0 (23 Nov 2026)

Refactoring for publication; the numerical method is unchanged.
A one-off regression script checked that the results match v0.1 (VALIDATION.md).

- One mesh builder, `half_pipe_model(...; notch = nothing)`, for the intact and the
  notched pipe (v0.1 had two). Both map index coordinates, so the intact mesh differs
  from v0.1 only by round-off in the node coordinates.
- `run_explicit` reduced to the time march: the whole-field L2 error and peak sampling
  are gone (the validation uses the sensors), and `material`/`geometry` are no longer
  arguments. Restart logic split into `initial_state`, `resumed_state` and `march!`.
- `ExplicitResult.sensor_velocity`, which held displacement in the notched runs, is now
  `sensor_history`, with `recorded_quantity` saying which. Sensors are named
  (`SurfaceSensor`) and record the three components; `cylindrical_components` projects
  them.
- The cases of the paper are functions with cache-free inputs (`IntactPipeCase`,
  `NotchedPipeCase`) instead of scripts with globals; the notched case covers L = 1, 2,
  3 and 10 m.
- The signal comparison (NRMSE, scattered-wave window, cross-correlation delay) is in
  Julia (`signals.jl`, `compare_signal_sets`) instead of a Python script.
- Tests with `Test` in `test/runtests.jl` (`Pkg.test()`), replacing the test scripts.
- Reference signals and costs copied into `data/`; no paths to sibling projects.
- Results written as text with provenance (date, git commit, Julia, threads).
- History and measured results moved from the README to VALIDATION.md.
- Diagnostic scripts of the debugging phase (face tags, end-face integration,
  NRMSE floor) not carried over; their findings are in VALIDATION.md.
- `Random`, `Serialization`, `Dates` declared in Project.toml (v0.1 used the first two
  without declaring them).
- Checkpoints of v0.1 are not compatible (CHECKPOINT_VERSION = 2).

## 0.1.0 (21-22 Nov 2026)

Project `Explicit_FEM_Gridap_UGW`: first version, validation and all the runs reported
in VALIDATION.md.
