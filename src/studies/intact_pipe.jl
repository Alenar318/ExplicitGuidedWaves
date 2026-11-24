# Intact pipe: the cost series of COMSOL "Elastic Waves, Time Explicit"
# (Comparacion_FEM_explicito), validated against the closed-form T(0,1) solution.
#
#   pipe length L, 30 kHz 5-cycle burst, 1 MPa, t_max = 400 us * L / 0.5 m (reflections
#   from both ends included), sensors on the outer surface at z = L/5, L/2, 4L/5,
#   theta = 90 deg; error: NRMSE of v_theta in % of the range of the exact signal.
#
# The default mesh (order 2, one cell through the wall, 96 sectors on the half circle,
# 5 mm along the pipe) matches the accuracy of COMSOL's quartic DG with hx = 20 mm at
# L = 0.5 m (NRMSE 7e-4 to 9e-4 % against 1.0e-3 to 2.3e-3 %).

"Mesh and length of one intact-pipe run."
struct IntactPipeCase
    length::Float64          # m
    order::Int               # 1 or 2
    n_radial::Int
    n_circumferential::Int   # sectors on the half circumference, even
    n_axial::Int
end

"""
    IntactPipeCase(length; order = 2, n_radial = 1, n_circumferential = 96,
                   axial_size = 0.005, n_axial = round(Int, length / axial_size))

Case with `n_axial` equal cells along the pipe: give either the target cell size
`axial_size` (m) or the number of cells `n_axial`. The defaults are the mesh of the
length series.
"""
function IntactPipeCase(length::Float64; order::Int = 2, n_radial::Int = 1,
                        n_circumferential::Int = 96, axial_size::Float64 = 0.005,
                        n_axial::Int = round(Int, length / axial_size))
    return IntactPipeCase(length, order, n_radial, n_circumferential, n_axial)
end

function case_name(case::IntactPipeCase)
    return @sprintf("intact_L%gm_p%d_nr%d_nc%d_nz%d", case.length, case.order,
                    case.n_radial, case.n_circumferential, case.n_axial)
end

"Observation time per metre of pipe (s/m): 400 us for the 0.5 m pipe of the COMSOL series."
const INTACT_FINAL_TIME_PER_LENGTH = 400.0e-6 / 0.5

"Final time of a case: long enough for the reflections from both ends."
intact_final_time(case::IntactPipeCase) = INTACT_FINAL_TIME_PER_LENGTH * case.length

"Sensors at L/5, L/2 and 4L/5, at theta = 90 deg (a node for any even n_circumferential)."
function intact_sensors(case::IntactPipeCase)
    pipe_length = case.length
    return [SurfaceSensor("z_L_over_5", pipe_length / 5, 90.0),
            SurfaceSensor("z_L_over_2", pipe_length / 2, 90.0),
            SurfaceSensor("z_4L_over_5", 4 * pipe_length / 5, 90.0)]
end

"""
    run_intact_pipe(case; checkpoint_file = nothing, verbose = true)

Solve one intact-pipe case and compare v_theta at the sensors with
`analytic_surface_velocity`. Returns a NamedTuple with the case, the cost (dofs,
steps, time step, assembly and loop times, threads), the signals (`times`,
`numerical`, `exact`, one column per sensor) and `nrmse_percent` per sensor.
"""
function run_intact_pipe(case::IntactPipeCase; checkpoint_file = nothing,
                         verbose::Bool = true)
    material = reference_steel()
    burst = intact_case_burst()
    geometry = reference_pipe(case.length)
    sensors = intact_sensors(case)

    mesh, _ = half_pipe_model(geometry, case.n_radial, case.n_circumferential,
                              uniform_axial_grid(case.length, case.n_axial))
    assembly_start = time()
    model = build_explicit_model(mesh, material, geometry, burst; order = case.order,
                                 verbose = verbose)
    dofs = recorded_dofs(model.space_trial, geometry, sensors)
    assembly_time = time() - assembly_start

    result = run_explicit(model, burst; final_time = intact_final_time(case),
                          recorded_dofs = dofs, recorded_quantity = :velocity,
                          checkpoint_file = checkpoint_file, verbose = verbose)
    components = cylindrical_components(result.sensor_history, sensors)

    n_samples = length(result.sensor_times)
    numerical = zeros(n_samples, length(sensors))
    exact = zeros(n_samples, length(sensors))
    errors = zeros(length(sensors))
    for (column, sensor) in enumerate(sensors)
        numerical[:, column] = components[sensor.name].circumferential
        for (i, t) in enumerate(result.sensor_times)
            exact[i, column] = analytic_surface_velocity(geometry, material, burst,
                                                         sensor.axial_position, t)
        end
        errors[column] = nrmse_percent(numerical[:, column], exact[:, column])
    end
    return (case = case, n_dofs = model.n_dofs, n_steps = result.n_steps,
            time_step = result.time_step, assembly_time = assembly_time,
            loop_time = result.loop_time, threads = Threads.nthreads(),
            sensors = sensors, times = result.sensor_times,
            numerical = numerical, exact = exact, nrmse_percent = errors)
end

"Columns of results/intact_pipe.csv."
const INTACT_TABLE_HEADER = ["case", "length_m", "order", "n_radial", "n_circumferential",
                             "n_axial", "n_dofs", "n_steps", "time_step_s", "assembly_s",
                             "loop_s", "threads", "nrmse_1_percent", "nrmse_2_percent",
                             "nrmse_3_percent", "date"]

"""
    save_intact_result(result, results_dir)

Append the cost and errors to `results_dir/intact_pipe.csv` and write the sensor
traces (time, numerical and exact v_theta at the three sensors) to
`results_dir/intact/<case name>_traces.txt`.
"""
function save_intact_result(result, results_dir::AbstractString)
    case = result.case
    row = [case_name(case), case.length, case.order, case.n_radial, case.n_circumferential,
           case.n_axial, result.n_dofs, result.n_steps, result.time_step,
           round(result.assembly_time; digits = 2), round(result.loop_time; digits = 2),
           result.threads, result.nrmse_percent..., Dates.format(now(), "yyyy-mm-dd")]
    append_csv_row(joinpath(results_dir, "intact_pipe.csv"), INTACT_TABLE_HEADER, row)

    traces_file = joinpath(results_dir, "intact", case_name(case) * "_traces.txt")
    mkpath(dirname(traces_file))
    open(traces_file, "w") do io
        println(io, "% intact pipe $(case_name(case)): v_theta (m/s) at z = L/5, L/2, 4L/5")
        for line in provenance_lines()
            println(io, "% ", line)
        end
        println(io, "% columns: time (s), numerical 1, 2, 3, exact 1, 2, 3")
        writedlm(io, hcat(result.times, result.numerical, result.exact), ' ')
    end
    return traces_file
end

"""
Series of intact-pipe cases run by `scripts/run_intact_pipe.jl <series>`, all at
L = 0.5 m except "length".

  - "convergence": orders 1 and 2, circumference and axis refined together (L = 0.5 m).
  - "directions": order 2, one direction refined at a time. "A": axis fixed at 5 mm,
    circumference refined; "B": circumference fixed at 48 sectors, axis refined.
  - "length": the matched-accuracy mesh at L = 0.5, 1 and 2 m, for the cost scaling
    (about 3, 15 and 60 min with 8 threads).
"""
const INTACT_SERIES = Dict(
    # orders 1 and 2; (sectors, axial cells) refined together
    "convergence" => vcat(
        [IntactPipeCase(0.5; order = 1, n_circumferential = sectors, n_axial = cells)
         for (sectors, cells) in ((12, 50), (24, 100), (48, 200), (96, 400))],
        [IntactPipeCase(0.5; order = 2, n_circumferential = sectors, n_axial = cells)
         for (sectors, cells) in ((12, 50), (24, 100), (48, 200))]),
    # "A": 100 axial cells (5 mm), sectors refined; "B": 48 sectors, axial cells refined
    "directions" => vcat(
        [IntactPipeCase(0.5; n_circumferential = sectors, n_axial = 100) for sectors in (12, 24, 48, 96)],
        [IntactPipeCase(0.5; n_circumferential = 48, n_axial = cells) for cells in (50, 200, 400)]),
    "length" => [IntactPipeCase(pipe_length) for pipe_length in (0.5, 1.0, 2.0)],
)
