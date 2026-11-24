# Intact-pipe runs, checkpointed. From the project folder:
#
#   julia -t 8 --project=. scripts/run_intact_pipe.jl <series | length> [--list] [--redo]
#
#   series   "convergence", "directions" or "length" (INTACT_SERIES in
#            src/studies/intact_pipe.jl)
#   length   a pipe length in m, with the matched-accuracy mesh (e.g. 1.0)
#   --list   show the cases and which ones are done
#   --redo   run again cases already in results/intact_pipe.csv
#
# Each case appends a row to results/intact_pipe.csv and writes its traces to
# results/intact/. A checkpoint is written every 10 minutes to results/checkpoints/:
# if the computer goes off, run the same command again and the case resumes.

using ExplicitGuidedWaves
using Printf

const RESULTS_DIR = joinpath(dirname(@__DIR__), "results")

function completed_cases()
    table_file = joinpath(RESULTS_DIR, "intact_pipe.csv")
    if !isfile(table_file)
        return Set{String}()
    end
    return Set(string.(read_csv_table(table_file)["case"]))
end

function main(arguments)
    selections = filter(argument -> !startswith(argument, "--"), arguments)
    if isempty(selections)
        error("usage: run_intact_pipe.jl <convergence | directions | length | L in m> [--list] [--redo]")
    end
    selection = selections[1]
    if haskey(INTACT_SERIES, selection)
        cases = INTACT_SERIES[selection]
    else
        cases = [IntactPipeCase(parse(Float64, selection))]
    end
    done = completed_cases()

    if "--list" in arguments
        for case in cases
            status = case_name(case) in done ? "done" : "to run"
            @printf("%-40s %s\n", case_name(case), status)
        end
        return nothing
    end

    println("warm-up (compilation, not timed)...")
    warm_up_solver()
    for case in cases
        name = case_name(case)
        checkpoint = joinpath(RESULTS_DIR, "checkpoints", name * ".jls")
        if name in done && !("--redo" in arguments)
            println("\n$name: already in results/intact_pipe.csv, skipped (--redo to run again)")
            continue
        end
        if "--redo" in arguments
            rm(checkpoint; force = true)
        end
        println("\n=== $name")
        flush(stdout)
        result = run_intact_pipe(case; checkpoint_file = checkpoint)
        save_intact_result(result, RESULTS_DIR)
        @printf("dofs %d, steps %d, assembly %.1f s, loop %.1f s (%d threads)\n",
                result.n_dofs, result.n_steps, result.assembly_time, result.loop_time,
                result.threads)
        @printf("NRMSE of v_theta at L/5, L/2, 4L/5: %.3e / %.3e / %.3e %%\n",
                result.nrmse_percent...)
        flush(stdout)
    end
    return nothing
end

main(ARGS)
