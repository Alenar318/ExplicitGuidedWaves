# Semi-discrete system  M a + K u = f(t)  of linear elastodynamics:
#
#     K_ij = int( eps(phi_i) : C : eps(phi_j) ) dV        consistent stiffness
#     M_i  = int( rho phi_i . (1, 1, 1) ) dV              row-sum lumped mass
#     f_i  = f(t) int( t_hat . phi_i ) dA                 end traction, separable in space and time
#
# C is the isotropic elasticity tensor, sigma = lambda tr(eps) I + 2 mu eps.

"Everything the time loop needs, assembled once."
struct ExplicitModel
    space_test                       # Gridap TestFESpace
    space_trial                      # Gridap TrialFESpace
    stiffness                        # sparse K, symmetric (checked)
    lumped_mass::Vector{Float64}     # diagonal of M
    load_spatial::Vector{Float64}    # int( t_hat . phi_i ) dA, scaled by f(t) in the loop
    constrained_dofs::Vector{Int}    # u_x and u_z on the cut plane y = 0
    n_dofs::Int
end

"""
    build_explicit_model(mesh, material, geometry, burst; order, verbose = true)

Assemble stiffness, lumped mass and the spatial part of the load vector on `mesh`
(from `half_pipe_model`) with vector Lagrange elements of degree `order` (1 or 2).

  - Lumped mass, assembled directly as a vector: row i of the consistent mass sums to
    int( rho phi_i . (1, 1, 1) ), because the shape functions of each component add
    up to one. Same numbers as summing the rows of the assembled matrix (test/), without
    building it; at 1.4 million dofs that matrix and its workspace take several GB.
  - The load is separable in space and time, so its spatial part is assembled once.
  - The time loop computes K u as K' u (`symmetric_mul!`), so K must be symmetric. It
    is checked with two random vectors instead of forming K'.
"""
function build_explicit_model(mesh, material::ElasticMaterial, geometry::PipeGeometry,
                              burst::ToneBurst; order::Int, verbose::Bool = true)
    reference_element = ReferenceFE(lagrangian, VectorValue{3, Float64}, order)
    space_test = TestFESpace(mesh, reference_element; conformity = :H1)
    space_trial = TrialFESpace(space_test)
    measure_domain = Measure(Triangulation(mesh), 2 * order)

    labels = get_face_labeling(mesh)
    add_tag_from_tags!(labels, "loaded_end", [LOADED_END_TAG])
    measure_loaded_end = Measure(BoundaryTriangulation(mesh; tags = "loaded_end"), 2 * order)

    lambda = lame_first_parameter(material)
    mu = shear_modulus(material)
    density = material.density
    strain(u) = symmetric_gradient(u)
    stress(strain_tensor) = lambda * tr(strain_tensor) * one(strain_tensor) + 2 * mu * strain_tensor

    stiffness_form(u, v) = ∫(strain(v) ⊙ (stress ∘ strain(u))) * measure_domain
    stiffness = assemble_matrix(stiffness_form, space_trial, space_test)
    check_symmetry(stiffness)

    unit_vector = VectorValue(1.0, 1.0, 1.0)
    lumped_mass = assemble_vector(v -> ∫(density * (v ⋅ unit_vector)) * measure_domain, space_test)
    if any(lumped_mass .<= 0)
        error("the lumped mass has non-positive entries")
    end

    traction = torsional_traction_field(geometry, burst)
    load_spatial = assemble_vector(v -> ∫(traction ⋅ v) * measure_loaded_end, space_test)

    constrained_dofs = antisymmetry_plane_dofs(space_trial, geometry)
    if isempty(constrained_dofs)
        error("no degrees of freedom found on the cut plane y = 0; check the mesh")
    end
    if verbose
        @printf("constrained dofs: %d of %d on the antisymmetry plane\n",
                length(constrained_dofs), length(lumped_mass))
    end
    return ExplicitModel(space_test, space_trial, stiffness, lumped_mass, load_spatial,
                         constrained_dofs, length(lumped_mass))
end

"Stop if a' K b and b' K a differ beyond round-off for two random vectors."
function check_symmetry(stiffness)
    probe_a = randn(StableRNG(1), size(stiffness, 1))
    probe_b = randn(StableRNG(2), size(stiffness, 1))
    forward = dot(probe_a, stiffness * probe_b)
    backward = dot(probe_b, stiffness * probe_a)
    asymmetry = abs(forward - backward) / abs(forward)
    if asymmetry > 1.0e-10
        error("the stiffness matrix is not symmetric (relative asymmetry $asymmetry)")
    end
    return nothing
end

"""
    antisymmetry_plane_dofs(space, geometry)

Degrees of freedom held at zero on the cut plane y = 0: the x and z components at
every node with y = 0. For T(0,1) the displacement there is e_theta = +/- e_y, normal
to the plane, so the condition is exact. COMSOL's Antisymmetry feature imposes the same.

The degrees of freedom are found by interpolating known fields, which works for any
Lagrange order and does not depend on Gridap's entity numbering: interpolating
(y, y, y) gives the y coordinate of the node that owns each degree of freedom, and
(1, 0, 0) and (0, 0, 1) mark the x and z components.
"""
function antisymmetry_plane_dofs(space, geometry::PipeGeometry)
    tolerance = 1.0e-8 * geometry.outer_radius
    y_at_dof = node_coordinate_at_dofs(space, 2)
    is_x_component = component_mask(space, 1)
    is_z_component = component_mask(space, 3)
    on_plane = abs.(y_at_dof) .< tolerance
    return findall(on_plane .& (is_x_component .| is_z_component))
end

"Coordinate `direction` (1, 2, 3) of the node that owns each degree of freedom."
function node_coordinate_at_dofs(space, direction::Int)
    field = x -> VectorValue(x[direction], x[direction], x[direction])
    return get_free_dof_values(interpolate_everywhere(field, space))
end

"True for the degrees of freedom of component `direction` (1, 2, 3)."
function component_mask(space, direction::Int)
    unit = VectorValue(ntuple(j -> j == direction ? 1.0 : 0.0, 3)...)
    return get_free_dof_values(interpolate_everywhere(x -> unit, space)) .> 0.5
end
