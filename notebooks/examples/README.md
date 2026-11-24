# Examples

Each notebook solves or reads one problem, compares it with a reference and ends with
the error. They use the environment of `notebooks/` (Plots, CSV, DataFrames):

```
julia --project=notebooks -e "using Pkg; Pkg.instantiate()"
```

Start Julia with 8 threads (`julia -t 8`, or the Jupyter kernel set to 8 threads).

| Example | Shows | Validated against | Run time |
|---|---|---|---|
| `01_intact_pipe.ipynb` | a small intact pipe solved in the notebook; convergence by direction; cost against COMSOL's explicit solver and against the pipe length | closed-form T(0,1) solution with end reflections (`analytic_surface_velocity`) | ~1 min (section 1 solves; sections 2-3 read `results/`) |
| `02_notched_pipe.ipynb` | the notched pipe of the paper; optional finer meshes and the 3 m pipe | hybrid BEM, full 3D BEM and implicit COMSOL signals in `data/` | < 1 min (reads `results/`) |

## Results the notebooks read

The long runs go through `scripts/`; if a result is missing, the notebook stops with
the command that produces it.

| Notebook section | Command (from the project folder) | Run time, 8 threads |
|---|---|---|
| 01, section 2 | `julia -t 8 --project=. scripts/run_intact_pipe.jl directions` | ~30-40 min |
| 01, section 3 | `julia -t 8 --project=. scripts/run_intact_pipe.jl length` | ~1-2 h |
| 02, sections 1 and 4 | `julia -t 8 --project=. scripts/run_notched_pipe.jl 1 48 1.25` | ~10-20 min |
| 02, section 2 (optional) | `scripts/run_notched_pipe.jl 1 96 1.25` and `scripts/run_notched_pipe.jl 1 48 0.625` | ~30 min each |
| 02, section 3 (optional) | `julia -t 8 --project=. scripts/run_notched_pipe.jl 3 48 1.25` | ~1 h |

Run times depend on the state of the computer (VALIDATION.md, run-time variability):
compare times measured in one session.

## Where each reference is valid

| Reference | Valid for | Expected difference with this code |
|---|---|---|
| closed-form T(0,1) (`analytic_surface_velocity`) | intact pipe, torsional load proportional to r, any time (end reflections included) | NRMSE below 1e-3 % at the default mesh; limited by the polygonal cross-section |
| hybrid BEM (`data/bem/`) | notched pipe, L = 1 m, up to 500 us | 0.06 % downstream, 0.25 % in the echo (u_theta) |
| full 3D BEM (`data/bem/`) | same | as the hybrid BEM, within 0.01 % of it |
| implicit COMSOL (`data/comsol_implicit/`) | same | 0.5-3 %, mostly a delay: the numerical dispersion of the implicit model |

## Adding an example

1. Put the reference (closed form, or signals in `data/` with their source in
   `data/README.md`) in the module or in `data/`, never in the notebook.
2. Copy `01_intact_pipe.ipynb`: a first markdown cell (what it shows, the problem,
   the reference, the run time), a setup cell, one step per cell, and a final cell with
   the error and how to read it.
3. Long runs: add the cases to `src/studies/` and a runner to `scripts/`; the notebook
   only reads their results.
4. Add a row to the tables above.
