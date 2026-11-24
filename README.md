# ExplicitGuidedWaves

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24503606.svg)](https://doi.org/10.5281/zenodo.24503606)

Explicit finite element solver for transient torsional guided waves in pipes with
localized defects, in Julia with [Gridap](https://github.com/gridap/Gridap.jl).

It is the open FEM reference of the hybrid axisymmetric/3D BEM paper: the same pipe,
notch, tone-burst load and sensors as the hybrid BEM and the COMSOL models, so a
reviewer can re-run the FEM side of every comparison.

## Method

| | |
|---|---|
| discretization | continuous Galerkin, displacement formulation, Lagrange hexahedra of order 1 or 2 |
| mesh | structured (r, theta, z) grid of half the pipe, 0 <= theta <= pi; straight-sided cells |
| mass | row-sum lumped (assembled directly as a vector) |
| time integration | central differences, dt = 0.8 x 2 / omega_max, omega_max by power iteration |
| symmetry | antisymmetry on y = 0 (u_x = u_z = 0), exact for T(0,1) |
| defect | notch as a box in polar coordinates, removed cell by cell (exact) |
| load | torsional traction on the end z = L, proportional to r, Hann-windowed tone burst |
| parallelism | sparse product and update on all Julia threads; results independent of the thread count |

Validated against the closed-form T(0,1) solution with end reflections (NRMSE below
1e-3 % at the default mesh) and against the hybrid BEM on the notched pipe (0.06 %
downstream, 0.25 % in the echo). Details and all measured costs: [VALIDATION.md](VALIDATION.md).

## Installation

Julia 1.11. From the project folder:

```
julia --project=. -e "using Pkg; Pkg.instantiate()"
julia -t 8 --project=. -e "using Pkg; Pkg.test()"
```

The tests take about two minutes.

## Use

Long runs go through two scripts. They write a checkpoint every 10 minutes: if the
computer goes off, run the same command again and the run resumes where it stopped.

```
julia -t 8 --project=. scripts/run_intact_pipe.jl length        # intact pipe, 0.5, 1 and 2 m
julia -t 8 --project=. scripts/run_intact_pipe.jl directions    # convergence by direction
julia -t 8 --project=. scripts/run_notched_pipe.jl 1 48 1.25    # notched pipe, L = 1 m
julia -t 8 --project=. scripts/run_notched_pipe.jl 3 48 1.25    # notched pipe, L = 3 m
```

| script argument | meaning | values |
|---|---|---|
| `run_intact_pipe.jl <series>` | series of cases in `INTACT_SERIES` | `convergence`, `directions`, `length`, or a length in m |
| `run_intact_pipe.jl ... --list` | show the cases and which are done | |
| `run_notched_pipe.jl <L> <nc> <dz>` | pipe length (m), sectors on the half circle, axial cell size across the notch (mm) | L = 1, 2, 3, 10; nc multiple of 4 |
| `run_notched_pipe.jl ... --redo` | discard the checkpoint and run again (to re-measure the time) | |

Results go to `results/` (tables `intact_pipe.csv` and `notched_pipe.csv`, signals in
`results/intact/` and `results/notch/<case>/`). The notebooks read them:

| notebook | shows | run time |
|---|---|---|
| `notebooks/examples/01_intact_pipe.ipynb` | a small case solved in the notebook; convergence and cost against COMSOL | ~1 min |
| `notebooks/examples/02_notched_pipe.ipynb` | this code against the hybrid BEM, the full 3D BEM and the implicit COMSOL; 3 m against 1 m | < 1 min (reads results) |

The notebooks use their own environment (Plots, CSV, DataFrames):

```
julia --project=notebooks -e "using Pkg; Pkg.instantiate()"
```

From Julia, one case directly:

```julia
using ExplicitGuidedWaves
result = run_intact_pipe(IntactPipeCase(0.5; n_circumferential = 24, axial_size = 0.005))
result.nrmse_percent        # NRMSE of v_theta at L/5, L/2, 4L/5, in %
```

## Project layout

| folder | content |
|---|---|
| `src/` | the module: `problem.jl` (material, pipe, notch, burst), `mesh.jl`, `analytic.jl`, `model.jl` (assembly), `time_integration.jl`, `sensors.jl`, `signals.jl`, `checkpoint.jl`, `results_io.jl` |
| `src/studies/` | the intact-pipe and notched-pipe cases of the paper |
| `scripts/` | terminal runners of the long cases, and the memory measurement |
| `test/` | `runtests.jl` |
| `data/` | reference signals and costs from COMSOL and the BEM (see `data/README.md`) |
| `notebooks/examples/` | examples that read `results/` |

## Limitations

- Straight pipes of circular cross-section and notches that are boxes in (r, theta, z);
  the grid lines must pass through the notch edges (checked).
- Torsional excitation with antisymmetry on y = 0: the model is exact for defects
  symmetric about that plane. Non-symmetric defects need the full annulus, which the
  mesh builder does not provide.
- The stiffness matrix is assembled, about 4.5 GB per million degrees of freedom; a
  laptop with 12 GB free holds a notched pipe of about 4.5 m.
- No absorbing boundaries: the pipe ends reflect.

See [CHANGELOG.md](CHANGELOG.md) for the changes from v0.1.

## License and citation

MIT (see [LICENSE](LICENSE)). Archived in Zenodo: https://doi.org/10.5281/zenodo.24503606. To cite the code, use [CITATION.cff](CITATION.cff); GitHub shows it as "Cite this repository".
