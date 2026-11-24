# Peak memory of the intact-pipe model, without the time loop. From the project folder, one length per process
# (the peak memory of a process only grows):
#
#   julia -t 8 --project=. scripts/measure_intact_memory.jl 2.0
#
# Builds what a run of scripts/run_intact_pipe.jl builds before its time loop (mesh, stiffness, lumped mass, load,
# constraints, sensors and the critical time step), plus the vectors the loop allocates, and reports the peak memory
# of the process (Sys.maxrss, which includes Julia and the compiled code, as in the notched-pipe runs). The time
# loop allocates nothing more, so this is the peak memory of the full run. Appends a row to
# results/intact_memory.csv.

using ExplicitGuidedWaves
using Printf
using Dates

const RESULTS_DIR = joinpath(dirname(@__DIR__), "results")

function main(arguments)
    pipe_length = isempty(arguments) ? 1.0 : parse(Float64, arguments[1])
    case = IntactPipeCase(pipe_length)            # the mesh of the length series: order 2, 96 sectors, 5 mm
    material = reference_steel()
    burst = intact_case_burst()
    geometry = reference_pipe(case.length)
    @printf("%s: building the model (no time loop)\n", case_name(case))
    flush(stdout)

    start = time()
    mesh, _ = half_pipe_model(geometry, case.n_radial, case.n_circumferential,
                              uniform_axial_grid(case.length, case.n_axial))
    model = build_explicit_model(mesh, material, geometry, burst; order = case.order)
    dofs = recorded_dofs(model.space_trial, geometry, intact_sensors(case))
    time_step = 0.8 * critical_time_step(model)
    loop_vectors = [zeros(model.n_dofs) for _ in 1:4]     # u_previous, u_current, u_next, K u of the time loop
    build_time = time() - start

    peak_gb = Sys.maxrss() / 1e9          # peak so far: the loop vectors are already allocated
    @printf("dofs %d (%d recorded), nonzeros of K %d, built in %.1f s, time step %.3e s\n",
            length(loop_vectors[1]), length(dofs), length(model.stiffness.nzval), build_time, time_step)
    @printf("peak memory of the process: %.2f GB (%.2f GB per million dofs)\n", peak_gb, peak_gb / (model.n_dofs / 1e6))
    header = ["case", "length_m", "n_dofs", "nonzeros", "peak_memory_gb", "build_s", "threads", "date"]
    row = [case_name(case), case.length, model.n_dofs, length(model.stiffness.nzval), round(peak_gb; digits = 2),
           round(build_time; digits = 1), Threads.nthreads(), Dates.format(now(), "yyyy-mm-dd")]
    append_csv_row(joinpath(RESULTS_DIR, "intact_memory.csv"), header, row)
    return nothing
end

main(ARGS)
