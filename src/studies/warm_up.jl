"""
    warm_up_solver()

Compile everything a run uses (mesh with and without notch, assembly, time loop) on a
tiny case of 2 us, so that compilation is not counted in the timings that follow.
"""
function warm_up_solver()
    material = reference_steel()
    burst = notch_case_burst()
    geometry = reference_pipe(1.0)
    notch = reference_notch()
    axial_edges = axial_grid([0.0, notch.z_start, notch.z_end, geometry.length], [0.05, 0.01, 0.05])
    mesh, _ = half_pipe_model(geometry, 2, 8, axial_edges; notch = notch)
    model = build_explicit_model(mesh, material, geometry, burst; order = 2, verbose = false)
    sensors = [SurfaceSensor("warm_up", geometry.length, 90.0)]
    dofs = recorded_dofs(model.space_trial, geometry, sensors)
    for quantity in (:displacement, :velocity)
        run_explicit(model, burst; final_time = 2.0e-6, recorded_dofs = dofs,
                     recorded_quantity = quantity, verbose = false)
    end
    return nothing
end
