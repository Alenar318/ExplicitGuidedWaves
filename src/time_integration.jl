# Central differences with a lumped mass (Hughes 2000, Chapter 9):
#
#     u_{n+1} = 2 u_n - u_{n-1} + dt^2 M^{-1} ( f(t_n) l - K u_n )
#
# with M diagonal, so every step costs one sparse matrix-vector product. Velocity, when
# recorded, is the centred difference (u_{n+1} - u_{n-1}) / (2 dt).
#
# Stable for dt <= 2 / omega_max, omega_max the highest natural frequency of M^{-1} K;
# the step used is safety_factor * 2 / omega_max, rounded down so that the last step
# lands exactly on the final time.
#
# The two kernels of a step run on all Julia threads (start Julia with `-t 8`). Their
# results do not depend on the number of threads.

"""
    critical_time_step(model; n_iterations = 200)

Stability limit 2 / omega_max of central differences. omega_max^2 is the largest
eigenvalue of M^{-1} K, estimated by power iteration with the M-weighted Rayleigh
quotient. The start vector comes from a StableRNG with a fixed seed, whose stream does
not change between Julia versions, so the step is reproducible.
"""
function critical_time_step(model::ExplicitModel; n_iterations::Int = 200)
    x = randn(StableRNG(20261122), model.n_dofs)
    x[model.constrained_dofs] .= 0.0   # constrained motions are not part of the system
    x ./= norm(x)
    largest_eigenvalue = 0.0
    for _ in 1:n_iterations
        y = (model.stiffness * x) ./ model.lumped_mass
        y[model.constrained_dofs] .= 0.0
        largest_eigenvalue = dot(x, model.stiffness * x) / dot(x, model.lumped_mass .* x)
        y_norm = norm(y)
        if y_norm == 0
            break
        end
        x = y ./ y_norm
    end
    omega_max = sqrt(max(largest_eigenvalue, 0.0))
    return 2.0 / omega_max
end

"""
    symmetric_mul!(y, A, x)

y = A x for a SYMMETRIC sparse matrix in CSC format, computed as y = A' x: entry j is
the dot product of column j with x. Columns are independent, so they are split across
threads without write conflicts (SparseArrays' `mul!` for CSC scatters into y and runs
on one thread). Each y[j] is always summed in the same order, so the result is
bitwise the same with any number of threads.

Valid only for symmetric A; `build_explicit_model` checks that the stiffness is.
"""
function symmetric_mul!(y::Vector{Float64}, A::SparseMatrixCSC{Float64}, x::Vector{Float64})
    column_pointers = A.colptr
    row_indices = A.rowval
    values = A.nzval
    Threads.@threads :static for column in 1:size(A, 2)
        column_sum = 0.0
        # @inbounds: the indices come from the CSC structure itself
        @inbounds for entry in column_pointers[column]:(column_pointers[column + 1] - 1)
            column_sum += values[entry] * x[row_indices[entry]]
        end
        @inbounds y[column] = column_sum
    end
    return y
end

"""
    central_difference_update!(u_next, u_current, u_previous, load_spatial,
                               internal_force, lumped_mass, load_factor, dt)

u_next = 2 u_current - u_previous + dt^2 (load_factor load_spatial - internal_force) / m,
entry by entry, on all threads. `internal_force` is K u_current.
"""
function central_difference_update!(u_next, u_current, u_previous, load_spatial,
                                    internal_force, lumped_mass, load_factor, dt)
    Threads.@threads :static for i in eachindex(u_next)
        # @inbounds: all vectors have length n_dofs
        @inbounds u_next[i] = 2 * u_current[i] - u_previous[i] +
                              dt^2 * (load_factor * load_spatial[i] - internal_force[i]) /
                              lumped_mass[i]
    end
    return u_next
end

"Result of one explicit run."
struct ExplicitResult
    time_step::Float64
    n_steps::Int
    loop_time::Float64                 # s, time loop only, summed over resumed segments
    sensor_times::Vector{Float64}      # s, t_n of every step
    sensor_history::Matrix{Float64}    # n_steps x number of recorded dofs
    recorded_quantity::Symbol          # :displacement (m) or :velocity (m/s)
end

"State of the time march; it is what a checkpoint stores."
mutable struct MarchState
    time_step::Float64
    n_steps::Int
    next_step::Int
    loop_time::Float64
    u_previous::Vector{Float64}
    u_current::Vector{Float64}
    sensor_times::Vector{Float64}
    sensor_history::Matrix{Float64}
end

"""
    run_explicit(model, burst; final_time, recorded_dofs,
                 recorded_quantity = :displacement, safety_factor = 0.8,
                 checkpoint_file = nothing, checkpoint_interval = 600.0,
                 max_steps_this_call = typemax(Int), verbose = true)

Integrate from rest to `final_time` (s), recording `recorded_quantity` (:displacement
or :velocity) of the degrees of freedom `recorded_dofs` at every step.

Restart: with `checkpoint_file`, the state is saved every `checkpoint_interval` seconds
of wall time and at the end. Calling again with the same file resumes where it stopped,
or returns at once if the run had finished. The saved time step is reused, so a resumed
run is bitwise identical to an uninterrupted one. A checkpoint of a different case
(dofs, sensors or final time) is refused.

`max_steps_this_call` stops after that many steps and returns `nothing` (used by the
restart test).

`loop_time` excludes assembly and checkpoint writes, so it is comparable with the
solution time that COMSOL reports.
"""
function run_explicit(model::ExplicitModel, burst::ToneBurst;
                      final_time::Float64,
                      recorded_dofs::Vector{Int},
                      recorded_quantity::Symbol = :displacement,
                      safety_factor::Float64 = 0.8,
                      checkpoint_file = nothing,
                      checkpoint_interval::Float64 = 600.0,
                      max_steps_this_call::Int = typemax(Int),
                      verbose::Bool = true)
    if !(recorded_quantity in (:displacement, :velocity))
        error("recorded_quantity must be :displacement or :velocity")
    end
    run_identity = (n_dofs = model.n_dofs, final_time = final_time,
                    recorded_dofs = recorded_dofs, recorded_quantity = recorded_quantity)

    state = nothing
    if checkpoint_file !== nothing
        state = resumed_state(checkpoint_file, run_identity, verbose)
    end
    if state === nothing
        state = initial_state(model, final_time, safety_factor, length(recorded_dofs), verbose)
    end
    if verbose
        @printf("step used       : %.4e s  (safety %.2f)\n", state.time_step, safety_factor)
        @printf("number of steps : %d\n", state.n_steps)
        @printf("degrees of freed: %d\n", model.n_dofs)
        @printf("threads         : %d\n", Threads.nthreads())
        flush(stdout)
    end

    # Count the steps left first: next_step + max_steps_this_call - 1 overflows for
    # the default typemax(Int).
    steps_left = state.n_steps - state.next_step + 1
    last_step = state.next_step - 1 + min(max_steps_this_call, steps_left)
    checkpoint_write_time = march!(state, model, burst, last_step, recorded_dofs,
                                   recorded_quantity, checkpoint_file, checkpoint_interval,
                                   run_identity, verbose)

    if last_step < state.n_steps
        if verbose
            @printf("stopped after step %d of %d; loop time so far %.1f s\n",
                    last_step, state.n_steps, state.loop_time)
        end
        return nothing
    end
    if verbose
        @printf("loop time       : %.1f s (checkpoint writes %.1f s excluded)\n",
                state.loop_time, checkpoint_write_time)
    end
    return ExplicitResult(state.time_step, state.n_steps, state.loop_time,
                          state.sensor_times, state.sensor_history, recorded_quantity)
end

"State at rest, with the time step from the stability limit."
function initial_state(model::ExplicitModel, final_time::Float64, safety_factor::Float64,
                       n_recorded::Int, verbose::Bool)
    critical_step = critical_time_step(model)
    n_steps = ceil(Int, final_time / (safety_factor * critical_step))
    time_step = final_time / n_steps   # the last step lands exactly on final_time
    if verbose
        @printf("critical step   : %.4e s\n", critical_step)
    end
    return MarchState(time_step, n_steps, 1, 0.0,
                      zeros(model.n_dofs), zeros(model.n_dofs),
                      zeros(n_steps), zeros(n_steps, n_recorded))
end

"State saved in `checkpoint_file`, or `nothing` if there is none."
function resumed_state(checkpoint_file::AbstractString, run_identity, verbose::Bool)
    saved = load_checkpoint(checkpoint_file)
    if saved === nothing
        return nothing
    end
    if saved.version != CHECKPOINT_VERSION
        error("checkpoint $(checkpoint_file) was written by another version of the solver; delete it")
    end
    if saved.run_identity != run_identity
        error("checkpoint $(checkpoint_file) belongs to a different case " *
              "(dofs, sensors or final time differ); delete it to start afresh")
    end
    state = saved.state
    if verbose
        if state.next_step > state.n_steps
            @printf("checkpoint holds a finished run: %d steps, %.1f s of loop time\n",
                    state.n_steps, state.loop_time)
        else
            @printf("RESUMING at step %d of %d (%.1f %%), %.1f s of loop time already done\n",
                    state.next_step, state.n_steps,
                    100 * (state.next_step - 1) / state.n_steps, state.loop_time)
        end
    end
    return state
end

"""
Advance `state` up to `last_step`, saving checkpoints on the way. Returns the time
spent writing checkpoints, which is not counted in `state.loop_time`.
"""
function march!(state::MarchState, model::ExplicitModel, burst::ToneBurst, last_step::Int,
                recorded_dofs::Vector{Int}, recorded_quantity::Symbol,
                checkpoint_file, checkpoint_interval::Float64, run_identity, verbose::Bool)
    if last_step < state.next_step
        return 0.0   # finished run read from a checkpoint: nothing to do
    end
    dt = state.time_step
    u_previous = state.u_previous
    u_current = state.u_current
    u_next = similar(u_current)
    internal_force = zeros(model.n_dofs)   # K u_current

    loop_time_before = state.loop_time
    segment_start = time()
    last_save = segment_start
    checkpoint_write_time = 0.0
    for step in state.next_step:last_step
        t = (step - 1) * dt
        symmetric_mul!(internal_force, model.stiffness, u_current)
        central_difference_update!(u_next, u_current, u_previous, model.load_spatial,
                                   internal_force, model.lumped_mass, burst_shape(burst, t), dt)
        # Antisymmetry on y = 0. Zeroing the updated vector imposes the homogeneous
        # constraint exactly and, unlike a penalty, leaves the stable step unchanged.
        u_next[model.constrained_dofs] .= 0.0

        state.sensor_times[step] = t
        record_sensors!(state.sensor_history, step, recorded_dofs, recorded_quantity,
                        u_previous, u_current, u_next, dt)

        u_previous, u_current, u_next = u_current, u_next, u_previous

        if checkpoint_file !== nothing && time() - last_save > checkpoint_interval
            write_start = time()
            state.next_step = step + 1
            state.loop_time = loop_time_before + (write_start - segment_start) - checkpoint_write_time
            state.u_previous = u_previous
            state.u_current = u_current
            save_checkpoint(checkpoint_file, (version = CHECKPOINT_VERSION,
                                              run_identity = run_identity, state = state))
            checkpoint_write_time += time() - write_start
            last_save = time()
            if verbose
                @printf("  checkpoint at step %d of %d (%.1f %%), loop time %.1f s\n",
                        step, state.n_steps, 100 * step / state.n_steps, state.loop_time)
                flush(stdout)
            end
        end
    end

    state.next_step = last_step + 1
    state.loop_time = loop_time_before + (time() - segment_start) - checkpoint_write_time
    state.u_previous = u_previous
    state.u_current = u_current
    if checkpoint_file !== nothing
        save_checkpoint(checkpoint_file, (version = CHECKPOINT_VERSION,
                                          run_identity = run_identity, state = state))
    end
    return checkpoint_write_time
end

"Store the recorded quantity of every recorded dof at step `step` (time t_n of u_current)."
function record_sensors!(history::Matrix{Float64}, step::Int, recorded_dofs::Vector{Int},
                         recorded_quantity::Symbol, u_previous::Vector{Float64},
                         u_current::Vector{Float64}, u_next::Vector{Float64}, dt::Float64)
    if recorded_quantity === :displacement
        for (column, dof) in enumerate(recorded_dofs)
            history[step, column] = u_current[dof]
        end
    else
        for (column, dof) in enumerate(recorded_dofs)
            history[step, column] = (u_next[dof] - u_previous[dof]) / (2 * dt)
        end
    end
    return nothing
end
