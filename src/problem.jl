# Problem data: material, pipe, notch, tone-burst excitation and torsional end load.
#
# The reference values are those of the hybrid BEM and COMSOL models of the paper, so
# the three codes solve the same problem.

"Linear isotropic elastic material."
struct ElasticMaterial
    density::Float64          # kg/m^3
    youngs_modulus::Float64   # Pa
    poisson_ratio::Float64    # -
end

"Steel of the reference models: 7850 kg/m^3, 200 GPa, Poisson ratio 0.3."
reference_steel() = ElasticMaterial(7850.0, 200.0e9, 0.3)

"Shear modulus mu = E / (2 (1 + nu))."
function shear_modulus(material::ElasticMaterial)
    return material.youngs_modulus / (2 * (1 + material.poisson_ratio))
end

"First Lame parameter lambda = E nu / ((1 + nu) (1 - 2 nu))."
function lame_first_parameter(material::ElasticMaterial)
    nu = material.poisson_ratio
    return material.youngs_modulus * nu / ((1 + nu) * (1 - 2 * nu))
end

"Shear wave speed c_T. The torsional mode T(0,1) travels at c_T without dispersion."
shear_wave_speed(material::ElasticMaterial) = sqrt(shear_modulus(material) / material.density)

"Longitudinal wave speed c_L, which sets the stability limit of the explicit scheme."
function longitudinal_wave_speed(material::ElasticMaterial)
    lambda = lame_first_parameter(material)
    mu = shear_modulus(material)
    return sqrt((lambda + 2 * mu) / material.density)
end

"Straight hollow cylinder along z, 0 <= z <= length. The load acts on the end z = length."
struct PipeGeometry
    inner_radius::Float64   # m
    outer_radius::Float64   # m
    length::Float64         # m
end

"Pipe of the reference models: outer radius 30 mm, wall thickness 4 mm."
reference_pipe(length::Float64) = PipeGeometry(0.026, 0.030, length)

wall_thickness(geometry::PipeGeometry) = geometry.outer_radius - geometry.inner_radius

"""
    NotchGeometry(depth, half_angle, z_start, z_end)

Notch on the outer surface. It is a box in polar coordinates:

    r     in [outer_radius - depth, outer_radius]
    theta in [0, half_angle]       (half model; the full notch spans -half_angle..half_angle)
    z     in [z_start, z_end]

so a structured (r, theta, z) mesh represents it exactly by removing cells, provided
grid lines pass through its edges. It is symmetric about y = 0, so the antisymmetry
condition on the cut plane stays exact.
"""
struct NotchGeometry
    depth::Float64        # m, from the outer surface
    half_angle::Float64   # rad, from the cut plane y = 0
    z_start::Float64      # m
    z_end::Float64        # m
end

"""
    reference_notch(; shift = 0.0)

Notch of the 1 m reference case (hybrid BEM `herramientas/entalla_T01.jl` and COMSOL
`comsol_entalla_T01_3D`): 2 mm deep (50 % of the wall), 90 deg wide in total, 10 mm
long, at z = 0.51375 ... 0.52375 m. `shift` moves it along z, for the longer pipes.
"""
function reference_notch(; shift::Float64 = 0.0)
    return NotchGeometry(0.002, pi / 4, 0.51375 + shift, 0.52375 + shift)
end

"""
    ToneBurst(center_frequency, n_cycles, traction_amplitude)

Hann-windowed tone burst, written in the three-sine form used by the COMSOL model and
the BEM codes, so the excitation is the same in the three codes:

    f(t) = 0.5 sin(2 pi fc t) - 0.25 sin(2 pi fc (1 + 1/Nc) t) - 0.25 sin(2 pi fc (1 - 1/Nc) t)

for 0 <= t <= Nc / fc, and zero outside. It equals sin(2 pi fc t) 0.5 (1 - cos(2 pi fc t / Nc)).
"""
struct ToneBurst
    center_frequency::Float64     # Hz
    n_cycles::Float64             # -
    traction_amplitude::Float64   # Pa
end

"Burst of the notched-pipe case: 40 kHz, 5 cycles, 1 MPa."
notch_case_burst() = ToneBurst(40.0e3, 5.0, 1.0e6)

"Burst of the intact-pipe cost series of COMSOL (Comparacion_FEM_explicito): 30 kHz, 5 cycles, 1 MPa."
intact_case_burst() = ToneBurst(30.0e3, 5.0, 1.0e6)

burst_duration(burst::ToneBurst) = burst.n_cycles / burst.center_frequency

"Dimensionless burst shape f(t); zero outside [0, n_cycles / center_frequency]."
function burst_shape(burst::ToneBurst, t::Float64)
    if t < 0.0 || t > burst_duration(burst)
        return 0.0
    end
    frequency = burst.center_frequency
    cycles = burst.n_cycles
    return 0.5 * sin(2 * pi * frequency * t) -
           0.25 * sin(2 * pi * frequency * (1 + 1 / cycles) * t) -
           0.25 * sin(2 * pi * frequency * (1 - 1 / cycles) * t)
end

"""
    torsional_traction_field(geometry, burst)

Spatial part of the torsional traction on the loaded end z = length:

    t(x, t) = amplitude (r / outer_radius) f(t) e_theta,    e_theta = (-y, x, 0) / r

The radius cancels, so the spatial part is (amplitude / outer_radius) (-y, x, 0). This
is the T(0,1) mode shape (u_theta proportional to r), so almost no other mode is excited.
"""
function torsional_traction_field(geometry::PipeGeometry, burst::ToneBurst)
    scale = burst.traction_amplitude / geometry.outer_radius
    return x -> VectorValue(-scale * x[2], scale * x[1], 0.0)
end
