# Tests of ExplicitGuidedWaves. From the project folder:
#
#   julia -t 8 --project=. -e "using Pkg; Pkg.test()"
#
# About two minutes. Every test compares against an independent result: a closed
# form, a quantity assembled another way, or an uninterrupted run.

using ExplicitGuidedWaves
using Gridap
using LinearAlgebra
using Random
using Test

const MATERIAL = reference_steel()
const INTACT_BURST = intact_case_burst()
const SHORT_PIPE = reference_pipe(0.5)

"Small intact model used by several tests (order 2: 22 725 dofs)."
function small_intact_model(; order::Int = 2)
    mesh, _ = half_pipe_model(SHORT_PIPE, 1, 12, uniform_axial_grid(SHORT_PIPE.length, 50))
    model = build_explicit_model(mesh, MATERIAL, SHORT_PIPE, INTACT_BURST; order = order,
                                 verbose = false)
    return mesh, model
end

@testset "ExplicitGuidedWaves" begin

    @testset "tone burst equals sin * Hann window" begin
        burst = notch_case_burst()
        for t in range(0.0, burst_duration(burst); length = 101)
            frequency = burst.center_frequency
            hann_form = sin(2 * pi * frequency * t) * 0.5 *
                        (1 - cos(2 * pi * frequency * t / burst.n_cycles))
            @test isapprox(burst_shape(burst, t), hann_form; atol = 1e-14)
        end
        @test burst_shape(burst, -1.0e-9) == 0.0
        @test burst_shape(burst, burst_duration(burst) + 1.0e-9) == 0.0
    end

    @testset "closed form before the first reflection is the travelling wave" begin
        wave_speed = shear_wave_speed(MATERIAL)
        amplitude = INTACT_BURST.traction_amplitude / (MATERIAL.density * wave_speed)
        z = 0.25
        for t in range(0.0, 2 * SHORT_PIPE.length / wave_speed * 0.9; length = 50)
            travelling = amplitude * burst_shape(INTACT_BURST, t - (SHORT_PIPE.length - z) / wave_speed)
            reflected = amplitude * burst_shape(INTACT_BURST, t - (SHORT_PIPE.length + z) / wave_speed)
            @test analytic_surface_velocity(SHORT_PIPE, MATERIAL, INTACT_BURST, z, t) ≈
                  travelling + reflected atol = 1e-15
        end
    end

    @testset "axial grid keeps every breakpoint" begin
        notch = reference_notch()
        breakpoints = [0.0, 0.25, notch.z_start, notch.z_end, 0.75, 1.0]
        edges = axial_grid(breakpoints, [0.005, 0.005, 0.00125, 0.005, 0.005])
        @test all(point in edges for point in breakpoints)
        @test issorted(edges)
        @test length(uniform_axial_grid(0.5, 100)) == 101
    end

    @testset "meshed volume equals the polygonal prism" begin
        mesh, removed = half_pipe_model(SHORT_PIPE, 1, 12, uniform_axial_grid(0.5, 50))
        @test removed == 0
        @test mesh_volume(mesh) ≈ polygonal_volume(SHORT_PIPE, 12) rtol = 1e-12

        geometry = reference_pipe(1.0)
        notch = reference_notch()
        edges = axial_grid([0.0, notch.z_start, notch.z_end, 1.0], [0.05, 0.0025, 0.05])
        notched_mesh, notched_removed = half_pipe_model(geometry, 2, 8, edges; notch = notch)
        @test notched_removed == 1 * 2 * 4   # radial x angular x axial notch cells
        @test mesh_volume(notched_mesh) ≈ polygonal_volume(geometry, 8; notch = notch) rtol = 1e-12
        @test_throws ErrorException half_pipe_model(geometry, 2, 6, edges; notch = notch)
    end

    @testset "lumped mass equals the row sums of the consistent mass" begin
        for order in (1, 2)
            mesh, model = small_intact_model(; order = order)
            measure = Measure(Triangulation(mesh), 2 * order)
            mass_form(u, v) = ∫(MATERIAL.density * (u ⋅ v)) * measure
            consistent = assemble_matrix(mass_form, model.space_trial, model.space_test)
            row_sums = vec(sum(consistent; dims = 2))
            @test maximum(abs.(row_sums .- model.lumped_mass)) / maximum(row_sums) < 1e-12
        end
    end

    @testset "threaded symmetric product equals K x" begin
        _, model = small_intact_model()
        x = randn(MersenneTwister(3), model.n_dofs)
        y = similar(x)
        symmetric_mul!(y, model.stiffness, x)
        @test maximum(abs.(y .- model.stiffness * x)) / maximum(abs.(y)) < 1e-12
    end

    @testset "resumed run is bitwise identical" begin
        _, model = small_intact_model()
        sensor = SurfaceSensor("middle", 0.25, 90.0)
        dofs = recorded_dofs(model.space_trial, SHORT_PIPE, [sensor])
        run_case(; options...) = run_explicit(model, INTACT_BURST; final_time = 100.0e-6,
                                              recorded_dofs = dofs, verbose = false, options...)
        uninterrupted = run_case()

        file = joinpath(mktempdir(), "checkpoint.jls")
        @test run_case(; checkpoint_file = file, max_steps_this_call = 150) === nothing
        @test run_case(; checkpoint_file = file, max_steps_this_call = 150) === nothing
        resumed = run_case(; checkpoint_file = file, checkpoint_interval = 0.0)
        rerun = run_case(; checkpoint_file = file)
        @test resumed.time_step == uninterrupted.time_step
        @test resumed.sensor_history == uninterrupted.sensor_history
        @test rerun.sensor_history == uninterrupted.sensor_history
    end

    @testset "intact pipe against the closed form (order 2, coarse mesh)" begin
        result = run_intact_pipe(IntactPipeCase(0.5; n_circumferential = 12, n_axial = 50); verbose = false)
        # Original project, 22 Nov 2026 (results/convergence_sensors.csv): 1.439e-2,
        # 1.745e-2, 1.912e-2 %. That run estimated the time step with an unseeded power
        # iteration (2.1505e-7 s), so rtol = 1e-2 rather than round-off.
        @test isapprox(result.nrmse_percent, [1.439369e-2, 1.744525e-2, 1.912065e-2]; rtol = 1e-2)
    end

    @testset "signal comparison tools" begin
        sample_interval = 1.0e-6
        times = collect(0.0:sample_interval:400.0e-6)
        burst = notch_case_burst()
        reference = [burst_shape(burst, t - 100.0e-6) for t in times]
        delayed = [burst_shape(burst, t - 100.3e-6) for t in times]
        @test cross_correlation_delay(delayed, reference, sample_interval) ≈ 0.3e-6 atol = 0.02e-6
        @test nrmse_percent(reference, reference) == 0.0
        fine_times = collect(0.0:0.25e-6:400.0e-6)
        fine = [burst_shape(burst, t - 100.0e-6) for t in fine_times]
        @test maximum(abs.(resample_on_grid(fine_times, fine, times) .- reference)) < 1e-2
    end
end
