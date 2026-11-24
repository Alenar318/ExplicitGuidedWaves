# Notched pipe: the case of the hybrid BEM paper (herramientas/entalla_T01.jl of the BEM
# project, COMSOL model comsol_entalla_T01_3D).
#
#   L = 1 m, T(0,1) excited at z = L by a 40 kHz 5-cycle burst, 1 MPa, t_max = 500 us;
#   notch 2 mm deep, 90 deg wide, 10 mm long at z = 0.51375 ... 0.52375 m;
#   sensors on the outer surface at z = 0.75 m (between load and notch: incident wave
#   and echo) and z = 0.25 m (transmitted wave), theta = 0, 90, 180 deg.
#
# Longer pipes, as in the BEM script for L = 10 m: the notch and the sensors move
# together by `notch_shift` along z, so the load is L - 1 - shift farther away. Until
# waves reflected by the pipe ends arrive, the signals are those of L = 1 m delayed by
# (L - 1 - shift) / c_T, because T(0,1) is not dispersive. That delayed 1 m signal is
# the reference for the long pipes (`delayed_to_1m_time`).
#
# Mesh: order 2, two cells through the wall (a radial node at the notch bottom), 5 mm
# along the pipe and `notch_axial_size` across the notch. Displacement is recorded,
# the quantity the implicit COMSOL and BEM signals give.

"Notch and sensor shift along z (m) for each supported pipe length (m)."
const NOTCH_CASE_SHIFTS = Dict(1.0 => 0.0, 2.0 => 0.5, 3.0 => 1.0, 10.0 => 4.0)   # 10 m: as in the BEM script

"Mesh parameters common to all notched cases."
const NOTCH_CASE_ORDER = 2
const NOTCH_CASE_RADIAL_CELLS = 2          # radial node at r = 28 mm, the notch bottom
const NOTCH_CASE_FAR_AXIAL_SIZE = 0.005    # m, away from the notch

"Sensor positions along z at L = 1 m (m): between load and notch (incident wave and echo), and past the notch."
const ECHO_SENSOR_Z_AT_1M = 0.75
const TRANSMITTED_SENSOR_Z_AT_1M = 0.25

"Final time at L = 1 m (s)."
const NOTCH_CASE_FINAL_TIME_1M = 500.0e-6

"""
Comparison windows, in L = 1 m time (s): the echo of the notch at z = 0.75 m (it
arrives at 224 us) and the transmitted wave at z = 0.25 m (240 us).
"""
const ECHO_WINDOW = (210.0e-6, 360.0e-6)
const TRANSMITTED_WINDOW = (230.0e-6, 380.0e-6)

"""
End of the interval where a longer pipe can be compared with L = 1 m, in L = 1 m time:
after 385 us the 1 m signals contain the notch echo reflected by the loaded end
(z = 0.75 m, front at 390 us) and the transmitted wave reflected by the far end
(z = 0.25 m, front at 399 us), which arrive much later in a longer pipe.
"""
const LONG_PIPE_VALID_END = 385.0e-6

"One notched-pipe run."
struct NotchedPipeCase
    length::Float64             # m, a key of NOTCH_CASE_SHIFTS
    n_circumferential::Int      # sectors on the half circumference, multiple of 4
    notch_axial_size::Float64   # m, axial cell size across the notch
end

"""
    NotchedPipeCase(; length = 1.0, n_circumferential = 48, notch_axial_size = 0.00125)

The defaults are the converged mesh: n_circumferential = 96 and notch_axial_size =
0.625 mm change u_theta by less than 0.02 % (VALIDATION.md).
"""
function NotchedPipeCase(; length::Float64 = 1.0, n_circumferential::Int = 48,
                         notch_axial_size::Float64 = 0.00125)
    if !haskey(NOTCH_CASE_SHIFTS, length)
        error("pipe length must be one of $(sort(collect(keys(NOTCH_CASE_SHIFTS)))) m")
    end
    return NotchedPipeCase(length, n_circumferential, notch_axial_size)
end

function case_name(case::NotchedPipeCase)
    return @sprintf("notch_L%gm_nc%d_dz%.3gmm", case.length, case.n_circumferential,
                    1000 * case.notch_axial_size)
end

notch_shift(case::NotchedPipeCase) = NOTCH_CASE_SHIFTS[case.length]

"Extra travel time from the load to the notch with respect to L = 1 m (s)."
function notch_extra_delay(case::NotchedPipeCase, material::ElasticMaterial)
    extra_distance = case.length - 1.0 - notch_shift(case)
    return extra_distance / shear_wave_speed(material)
end

"Six sensors: echo side (z = 0.75 m + shift) and transmitted side (0.25 m + shift), theta = 0, 90, 180."
function notched_sensors(case::NotchedPipeCase)
    shift = notch_shift(case)
    sensors = SurfaceSensor[]
    for (side, z_at_1m) in (("echo", ECHO_SENSOR_Z_AT_1M), ("transmitted", TRANSMITTED_SENSOR_Z_AT_1M))
        for angle in (0.0, 90.0, 180.0)
            name = @sprintf("%s_theta%d", side, round(Int, angle))
            push!(sensors, SurfaceSensor(name, z_at_1m + shift, angle))
        end
    end
    return sensors
end

"""
    run_notched_pipe(case; checkpoint_file = nothing, verbose = true)

Build the notched mesh, check it (cells removed and volume against the exact
polygonal volume), solve, and return a NamedTuple with the case, cost (dofs, steps,
time step, assembly and loop times, threads), `times` on a 1 us grid and `signals`:
Dict sensor name => (radial, circumferential, axial) displacement (m).
"""
function run_notched_pipe(case::NotchedPipeCase; checkpoint_file = nothing,
                          verbose::Bool = true)
    material = reference_steel()
    burst = notch_case_burst()
    geometry = reference_pipe(case.length)
    shift = notch_shift(case)
    notch = reference_notch(; shift = shift)
    sensors = notched_sensors(case)
    final_time = NOTCH_CASE_FINAL_TIME_1M + notch_extra_delay(case, material)

    echo_z = ECHO_SENSOR_Z_AT_1M + shift
    transmitted_z = TRANSMITTED_SENSOR_Z_AT_1M + shift
    axial_edges = axial_grid([0.0, transmitted_z, notch.z_start, notch.z_end, echo_z, case.length],
                             [NOTCH_CASE_FAR_AXIAL_SIZE, NOTCH_CASE_FAR_AXIAL_SIZE,
                              case.notch_axial_size, NOTCH_CASE_FAR_AXIAL_SIZE,
                              NOTCH_CASE_FAR_AXIAL_SIZE])
    mesh, n_removed = half_pipe_model(geometry, NOTCH_CASE_RADIAL_CELLS, case.n_circumferential,
                                      axial_edges; notch = notch)
    check_notched_mesh(mesh, geometry, notch, case, axial_edges, n_removed, verbose)

    assembly_start = time()
    model = build_explicit_model(mesh, material, geometry, burst; order = NOTCH_CASE_ORDER,
                                 verbose = verbose)
    dofs = recorded_dofs(model.space_trial, geometry, sensors)
    assembly_time = time() - assembly_start
    if verbose
        @printf("assembly            : %.1f s, peak memory so far %.2f GB\n",
                assembly_time, Sys.maxrss() / 1e9)
        flush(stdout)
    end

    result = run_explicit(model, burst; final_time = final_time, recorded_dofs = dofs,
                          recorded_quantity = :displacement,
                          checkpoint_file = checkpoint_file, verbose = verbose)
    components = cylindrical_components(result.sensor_history, sensors)

    output_times = collect(0.0:1.0e-6:result.sensor_times[end])
    signals = Dict{String, NamedTuple}()
    for sensor in sensors
        component = components[sensor.name]
        signals[sensor.name] = (
            radial = resample_on_grid(result.sensor_times, component.radial, output_times),
            circumferential = resample_on_grid(result.sensor_times, component.circumferential, output_times),
            axial = resample_on_grid(result.sensor_times, component.axial, output_times),
        )
    end
    return (case = case, n_dofs = model.n_dofs, n_steps = result.n_steps,
            time_step = result.time_step, assembly_time = assembly_time,
            loop_time = result.loop_time, threads = Threads.nthreads(),
            peak_memory_gb = Sys.maxrss() / 1e9, sensors = sensors,
            times = output_times, signals = signals)
end

"Stop before the long run if the notch is not exactly where it should be."
function check_notched_mesh(mesh, geometry::PipeGeometry, notch::NotchGeometry,
                            case::NotchedPipeCase, axial_edges::Vector{Float64},
                            n_removed::Int, verbose::Bool)
    n_notch_axial_cells = count(i -> notch.z_start <= axial_edges[i] < notch.z_end,
                                1:(length(axial_edges) - 1))
    expected_removed = (NOTCH_CASE_RADIAL_CELLS ÷ 2) * (case.n_circumferential ÷ 4) *
                       n_notch_axial_cells
    volume = mesh_volume(mesh)
    exact_volume = polygonal_volume(geometry, case.n_circumferential; notch = notch)
    relative_difference = abs(volume - exact_volume) / exact_volume
    if verbose
        @printf("notched pipe: L = %.2f m, notch at z = %.5f-%.5f m\n",
                geometry.length, notch.z_start, notch.z_end)
        @printf("mesh: order %d, nr = %d, nc = %d, nz = %d (far %.1f mm, notch %.3g mm)\n",
                NOTCH_CASE_ORDER, NOTCH_CASE_RADIAL_CELLS, case.n_circumferential,
                length(axial_edges) - 1, 1000 * NOTCH_CASE_FAR_AXIAL_SIZE,
                1000 * case.notch_axial_size)
        @printf("notch cells removed : %d (expected %d)\n", n_removed, expected_removed)
        @printf("volume              : %.9e m^3, exact %.9e m^3, difference %.2e\n",
                volume, exact_volume, relative_difference)
        flush(stdout)
    end
    if n_removed != expected_removed || relative_difference > 1e-10
        error("the notch is not where it should be; stopping before the long run")
    end
    return nothing
end

"""
    save_notched_result(result, results_dir)

Write the six sensor signals to `results_dir/notch/<case name>/<sensor>.txt` and
append the cost to `results_dir/notched_pipe.csv`.
"""
function save_notched_result(result, results_dir::AbstractString)
    case = result.case
    folder = joinpath(results_dir, "notch", case_name(case))
    for sensor in result.sensors
        description = [@sprintf("Gridap explicit FEM, %s, sensor %s at z = %.4f m, theta = %g deg",
                                case_name(case), sensor.name, sensor.axial_position, sensor.angle_deg),
                       @sprintf("%d dofs, %d steps, loop %.1f s", result.n_dofs, result.n_steps,
                                result.loop_time),
                       "displacement (m)"]
        write_signal_file(joinpath(folder, sensor.name * ".txt"), result.times,
                          result.signals[sensor.name], description)
    end
    header = ["case", "length_m", "n_circumferential", "notch_axial_size_mm", "n_dofs",
              "n_steps", "time_step_s", "assembly_s", "loop_s", "threads", "peak_memory_gb", "date"]
    row = [case_name(case), case.length, case.n_circumferential, 1000 * case.notch_axial_size,
           result.n_dofs, result.n_steps, result.time_step,
           round(result.assembly_time; digits = 2), round(result.loop_time; digits = 2),
           result.threads, round(result.peak_memory_gb; digits = 2),
           Dates.format(now(), "yyyy-mm-dd")]
    append_csv_row(joinpath(results_dir, "notched_pipe.csv"), header, row)
    return folder
end

"""
    load_notched_signals(case, results_dir)

Signals of a finished run (`save_notched_result`): (times, signals) with signals a
Dict sensor name => (radial, circumferential, axial). Stops with the command that
produces them if they are missing.
"""
function load_notched_signals(case::NotchedPipeCase, results_dir::AbstractString)
    folder = joinpath(results_dir, "notch", case_name(case))
    if !isdir(folder)
        error("no results for $(case_name(case)). Run, from the project folder:\n" *
              @sprintf("  julia -t 8 --project=. scripts/run_notched_pipe.jl %g %d %g",
                       case.length, case.n_circumferential, 1000 * case.notch_axial_size))
    end
    signals = Dict{String, NamedTuple}()
    times = Float64[]
    for sensor in notched_sensors(case)
        signal = read_signal_file(joinpath(folder, sensor.name * ".txt"))
        times = signal.time
        signals[sensor.name] = (radial = signal.radial, circumferential = signal.circumferential,
                                axial = signal.axial)
    end
    return times, signals
end

"""
    reference_signals(source, data_dir)

Reference signals of the 1 m notched case, as (times, signals) like
`load_notched_signals`. `source` is :bem_hybrid, :bem_full_3d or :comsol_implicit;
the files are in `data_dir` (see data/README.md).
"""
function reference_signals(source::Symbol, data_dir::AbstractString)
    file_patterns = Dict(
        :bem_hybrid => ("bem", "bem_hibrido_entalla_x%s_th%d.txt"),
        :bem_full_3d => ("bem", "bem_completo_entalla_x%s_th%d.txt"),
        :comsol_implicit => ("comsol_implicit", "entalla_T01_L1_40kHz_x%s_th%d.txt"),
    )
    if !haskey(file_patterns, source)
        error("unknown reference $source; use one of $(collect(keys(file_patterns)))")
    end
    folder, pattern = file_patterns[source]
    format = Printf.Format(pattern)
    signals = Dict{String, NamedTuple}()
    times = Float64[]
    for sensor in notched_sensors(NotchedPipeCase())
        # the reference files are named by the sensor position at L = 1 m (x0.75, x0.25)
        side_position = @sprintf("%g", sensor.axial_position)
        file = Printf.format(format, side_position, round(Int, sensor.angle_deg))
        signal = read_signal_file(joinpath(data_dir, folder, file))
        times = signal.time
        signals[sensor.name] = (radial = signal.radial, circumferential = signal.circumferential,
                                axial = signal.axial)
    end
    return times, signals
end

"""
    delayed_to_1m_time(case, times, signals)

Signals of a longer pipe moved to L = 1 m time: sampled at t_1m + extra delay for
t_1m = 0, 1, ..., 385 us (`LONG_PIPE_VALID_END`). Returns (times_1m, signals).
"""
function delayed_to_1m_time(case::NotchedPipeCase, times::Vector{Float64}, signals)
    delay = notch_extra_delay(case, reference_steel())
    times_1m = collect(0.0:1.0e-6:LONG_PIPE_VALID_END)
    moved = Dict{String, NamedTuple}()
    for (name, signal) in signals
        moved[name] = (
            radial = resample_on_grid(times, signal.radial, times_1m .+ delay),
            circumferential = resample_on_grid(times, signal.circumferential, times_1m .+ delay),
            axial = resample_on_grid(times, signal.axial, times_1m .+ delay),
        )
    end
    return times_1m, moved
end

"""
    compare_signal_sets(times, signals, reference_times, reference_signals)

Compare two sets of notched-pipe signals, sensor by sensor, on their common time
interval (both on the same 1 us grid starting at 0). Returns one NamedTuple per sensor:

  - `radial`, `circumferential`, `axial`: NRMSE (%) of each component, all normalised
    by the range of the reference u_theta, as in the paper;
  - `window`: NRMSE (%) of u_theta in the scattered-wave window (`ECHO_WINDOW` or
    `TRANSMITTED_WINDOW`), normalised by the reference range in that window;
  - `delay_ns`: delay of u_theta by cross-correlation (positive: later than the reference).
"""
function compare_signal_sets(times::Vector{Float64}, signals,
                             reference_times::Vector{Float64}, reference)
    n_common = min(length(times), length(reference_times))
    if maximum(abs.(times[1:n_common] .- reference_times[1:n_common])) > 1e-9
        error("compare_signal_sets: the two signal sets are not on the same time grid")
    end
    common_times = times[1:n_common]
    sample_interval = common_times[2] - common_times[1]
    rows = NamedTuple[]
    for sensor in notched_sensors(NotchedPipeCase())
        signal = signals[sensor.name]
        reference_signal = reference[sensor.name]
        u_theta = signal.circumferential[1:n_common]
        u_theta_reference = reference_signal.circumferential[1:n_common]
        scale = signal_range(u_theta_reference)

        window = startswith(sensor.name, "echo") ? ECHO_WINDOW : TRANSMITTED_WINDOW
        in_window = (common_times .>= window[1]) .& (common_times .<= window[2])
        push!(rows, (
            sensor = sensor.name,
            circumferential = nrmse_percent(u_theta, u_theta_reference),
            radial = nrmse_percent(signal.radial[1:n_common], reference_signal.radial[1:n_common];
                                   scale = scale),
            axial = nrmse_percent(signal.axial[1:n_common], reference_signal.axial[1:n_common];
                                  scale = scale),
            window = nrmse_percent(u_theta[in_window], u_theta_reference[in_window]),
            delay_ns = 1e9 * cross_correlation_delay(u_theta, u_theta_reference, sample_interval),
        ))
    end
    return rows
end
