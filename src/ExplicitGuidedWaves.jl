"""
    ExplicitGuidedWaves

Explicit finite element solver for transient torsional guided waves in pipes with
localized defects, built on Gridap. It is the open, reproducible FEM reference of the
hybrid axisymmetric/3D BEM paper: the same pipe, notch, load and sensors as the hybrid
BEM and the COMSOL models.

Method

  - Continuous Galerkin, displacement formulation, Lagrange hexahedra of order 1 or 2
    on a structured (r, theta, z) grid of HALF the pipe, 0 <= theta <= pi.
  - Consistent stiffness, row-sum lumped mass.
  - Central differences (Newmark beta = 0, gamma = 1/2) with time step
    0.8 * 2 / omega_max, omega_max estimated by power iteration.
  - Antisymmetry on the cut plane y = 0 (u_x = u_z = 0), exact for the torsional
    mode T(0,1).

Reference for the time integration: T.J.R. Hughes, The Finite Element Method: Linear
Static and Dynamic Finite Element Analysis, Dover, 2000, Chapter 9.

Use

  - `run_intact_pipe(IntactPipeCase(...))`: intact pipe, validated against the
    closed-form T(0,1) solution with end reflections (`analytic_surface_velocity`).
  - `run_notched_pipe(NotchedPipeCase(...))`: notched pipe, compared with the hybrid
    BEM and the implicit COMSOL signals in `data/`.

Long runs go through `scripts/` (checkpointed, they survive a shutdown); the notebooks
in `notebooks/examples/` read their results. See README.md.
"""
module ExplicitGuidedWaves

using Gridap
using LinearAlgebra
using SparseArrays
using Printf
using Random
using StableRNGs
using Serialization
using DelimitedFiles
using Dates

include("problem.jl")
include("mesh.jl")
include("analytic.jl")
include("model.jl")
include("sensors.jl")
include("checkpoint.jl")
include("time_integration.jl")
include("signals.jl")
include("results_io.jl")
include("plot_colors.jl")
include("studies/intact_pipe.jl")
include("studies/notched_pipe.jl")
include("studies/warm_up.jl")

# Problem data
export ElasticMaterial, reference_steel
export shear_modulus, lame_first_parameter, shear_wave_speed, longitudinal_wave_speed
export PipeGeometry, reference_pipe, wall_thickness
export NotchGeometry, reference_notch
export ToneBurst, notch_case_burst, intact_case_burst, burst_duration, burst_shape

# Mesh
export axial_grid, uniform_axial_grid, half_pipe_model, polygonal_volume, mesh_volume

# Closed-form solution
export analytic_surface_velocity, arrival_time

# Model and time integration
export ExplicitModel, build_explicit_model
export critical_time_step, symmetric_mul!, central_difference_update!
export ExplicitResult, run_explicit

# Sensors and signals
export SurfaceSensor, sensor_point, node_component_dofs, recorded_dofs
export cylindrical_components, resample_on_grid
export signal_range, nrmse_percent, cross_correlation_delay

# Results on disk
export write_signal_file, read_signal_file, append_csv_row, read_csv_table

# Figures
export MODEL_COLORS

# Studies
export warm_up_solver
export IntactPipeCase, INTACT_SERIES, case_name, intact_sensors, run_intact_pipe, save_intact_result
export NotchedPipeCase, notch_shift, notch_extra_delay, notched_sensors
export run_notched_pipe, save_notched_result, load_notched_signals
export reference_signals, compare_signal_sets, delayed_to_1m_time

end # module
