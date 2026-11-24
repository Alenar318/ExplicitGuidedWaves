# Notched-pipe run, checkpointed. From the project folder:
#
#   julia -t 8 --project=. scripts/run_notched_pipe.jl [length_m] [n_circumferential] [notch_dz_mm] [--redo]
#
#   length_m           1 (default), 2, 3 or 10
#   n_circumferential  sectors on the half circumference, multiple of 4 (default 48)
#   notch_dz_mm        axial cell size across the notch (default 1.25)
#   --redo             discard the checkpoint and run again (e.g. to re-measure the time;
#                      without it, a finished checkpoint returns the old result at once)
#
# Writes the six sensor signals to results/notch/<case>/ and appends the cost to
# results/notched_pipe.csv. A checkpoint is written every 10 minutes to
# results/checkpoints/: if the computer goes off, run the same command again.
# The comparison with the BEM and COMSOL is in notebooks/examples/02_notched_pipe.ipynb.

using ExplicitGuidedWaves
using Printf

const RESULTS_DIR = joinpath(dirname(@__DIR__), "results")

function main(all_arguments)
    arguments = filter(argument -> !startswith(argument, "--"), all_arguments)
    pipe_length = length(arguments) >= 1 ? parse(Float64, arguments[1]) : 1.0
    n_circumferential = length(arguments) >= 2 ? parse(Int, arguments[2]) : 48
    notch_axial_size = length(arguments) >= 3 ? parse(Float64, arguments[3]) / 1000 : 0.00125
    case = NotchedPipeCase(; length = pipe_length, n_circumferential = n_circumferential,
                           notch_axial_size = notch_axial_size)

    println("warm-up (compilation, not timed)...")
    warm_up_solver()
    println("\n=== ", case_name(case))
    flush(stdout)
    checkpoint = joinpath(RESULTS_DIR, "checkpoints", case_name(case) * ".jls")
    if "--redo" in all_arguments
        rm(checkpoint; force = true)
    end
    result = run_notched_pipe(case; checkpoint_file = checkpoint)
    folder = save_notched_result(result, RESULTS_DIR)

    @printf("\ndofs %d, steps %d, assembly %.1f s, loop %.1f s (%.2f h), %d threads\n",
            result.n_dofs, result.n_steps, result.assembly_time, result.loop_time,
            result.loop_time / 3600, result.threads)
    @printf("peak memory of the process: %.2f GB\n", result.peak_memory_gb)
    println("signals written to ", folder)
    return nothing
end

main(ARGS)
