# Closed-form T(0,1) response of an intact pipe, used to validate the solver.
#
# T(0,1) is non-dispersive with u_theta proportional to r, so the elastodynamic
# equations reduce to the 1D wave equation along z with speed c_T. A travelling wave
# loaded by the end traction tau(t) has particle velocity
#
#     v_theta(r, z, t) = tau0 (r / ro) / (rho c_T) f(t - (L - z) / c_T)
#
# (shear-impedance relation). Reflections from the free ends are added by the method
# of images.

"Time for the wavefront to travel from the loaded end z = length to the station z."
function arrival_time(geometry::PipeGeometry, material::ElasticMaterial, z::Float64)
    return (geometry.length - z) / shear_wave_speed(material)
end

"""
    analytic_surface_velocity(geometry, material, burst, z, t)

Circumferential velocity v_theta (m/s) on the outer surface at axial position z,
including the reflections from both ends.

Method of images. The load acts on z = L; z = 0 is free. Once the burst is over both
ends are traction free, and a shear wave reflects from a free end with velocity
coefficient +1. With A = tau0 / (rho c_T) at r = ro:

    v_theta(z, t) = A sum_{n >= 0} [ f(t - ((2n+1) L - z) / c_T) + f(t - ((2n+1) L + z) / c_T) ]

Valid while the burst ends before the first reflection returns to the loaded end,
burst_duration < 2 L / c_T (checked).
"""
function analytic_surface_velocity(geometry::PipeGeometry, material::ElasticMaterial,
                                   burst::ToneBurst, z::Float64, t::Float64)
    wave_speed = shear_wave_speed(material)
    pipe_length = geometry.length
    if burst_duration(burst) >= 2 * pipe_length / wave_speed
        error("the burst is still on when the first reflection returns to the loaded end")
    end
    amplitude = burst.traction_amplitude / (material.density * wave_speed)

    total = 0.0
    n_reflection = 0
    delay_direct = (pipe_length - z) / wave_speed   # path from the loaded end
    while delay_direct <= t
        delay_free_end = ((2 * n_reflection + 1) * pipe_length + z) / wave_speed   # one more reflection at z = 0
        total += burst_shape(burst, t - delay_direct)
        total += burst_shape(burst, t - delay_free_end)
        n_reflection += 1
        delay_direct = ((2 * n_reflection + 1) * pipe_length - z) / wave_speed
    end
    return amplitude * total
end
