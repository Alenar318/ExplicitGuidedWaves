# Structured hexahedral meshes of half a pipe (0 <= theta <= pi), intact or notched.
#
# The mesh is a Cartesian grid in index coordinates (i_r, i_theta, i_z) mapped to
# physical space: uniform in r and theta, piecewise uniform in z (`axial_grid`). The
# map moves the NODES; the cells keep straight faces, so the pipe is an inscribed
# polygonal prism. With n sectors on the half circumference the volume deficit is
# about (pi / n)^2 / 6 (`polygonal_volume` gives the exact meshed volume).
#
# Half the pipe, like the COMSOL and BEM models: a full annulus would double the
# degrees of freedom and the reported cost of the explicit method.
#
# Face labels of the Cartesian grid. Gridap numbers the six faces of the parametric
# box from the LAST coordinate to the first:
#   tag_21 / tag_22 : z = 0 / z = length (the loaded end)
#   tag_23 / tag_24 : theta = 0 / theta = pi (the two halves of the cut plane y = 0)
#   tag_25 / tag_26 : inner / outer surface

"Face label of the loaded end z = length (see the header of mesh.jl)."
const LOADED_END_TAG = "tag_22"

"""
    axial_grid(breakpoints, target_sizes)

Cell edges along the pipe. Every breakpoint becomes a node (the ends, the sensor
positions, the notch edges) and segment i, between breakpoints i and i + 1, is split
into round(length / target_sizes[i]) equal cells, at least one.
"""
function axial_grid(breakpoints::Vector{Float64}, target_sizes::Vector{Float64})
    if length(target_sizes) != length(breakpoints) - 1
        error("axial_grid: need one target size per segment between breakpoints")
    end
    if !issorted(breakpoints)
        error("axial_grid: breakpoints must be sorted")
    end
    edges = [breakpoints[1]]
    for segment in eachindex(target_sizes)
        segment_start = breakpoints[segment]
        segment_end = breakpoints[segment + 1]
        n_cells = max(1, round(Int, (segment_end - segment_start) / target_sizes[segment]))
        for cell in 1:n_cells
            if cell == n_cells
                push!(edges, segment_end)   # exactly the breakpoint, no round-off
            else
                push!(edges, segment_start + (segment_end - segment_start) * cell / n_cells)
            end
        end
    end
    return edges
end

"Cell edges of `n_axial` equal cells over 0 <= z <= length."
uniform_axial_grid(length::Float64, n_axial::Int) = axial_grid([0.0, length], [length / n_axial])

"""
    half_pipe_model(geometry, n_radial, n_circumferential, axial_edges; notch = nothing)

Hexahedral mesh of half the pipe, optionally with a notch removed. Returns
`(model, n_removed_cells)`.

  - `n_radial` cells through the wall and `n_circumferential` sectors on the half
    circumference, both uniform; `axial_edges` from `axial_grid`.
  - With a notch, the cells whose centroid lies inside it are removed. The grid lines
    must pass through the notch edges (checked), so the notch is represented exactly.
    The removed cells leave traction-free faces, and the face labels of the remaining
    model are kept (the loaded end is still `LOADED_END_TAG`).

The model is converted to an unstructured one. A `CartesianDiscreteModel` assumes a
constant Jacobian on the faces normal to the axis; under the polar map that integrates
the loaded end as if every cell had the width of its inner edge, i.e. a load short by
2 ri / (ri + ro) = 0.929 with one radial layer. The unstructured model uses the
general face geometry.
"""
function half_pipe_model(geometry::PipeGeometry, n_radial::Int, n_circumferential::Int,
                         axial_edges::Vector{Float64}; notch = nothing)
    if notch !== nothing
        check_notch_alignment(geometry, notch, n_radial, n_circumferential, axial_edges)
    end
    inner_radius = geometry.inner_radius
    outer_radius = geometry.outer_radius
    n_axial = length(axial_edges) - 1

    # index coordinates (i_r, i_theta, i_z) -> (x, y, z)
    function index_to_physical(index)
        radius = inner_radius + (outer_radius - inner_radius) * index[1] / n_radial
        angle = pi * index[2] / n_circumferential
        cell = clamp(floor(Int, index[3]), 0, n_axial - 1)
        fraction = index[3] - cell
        z = axial_edges[cell + 1] + fraction * (axial_edges[cell + 2] - axial_edges[cell + 1])
        return VectorValue(radius * cos(angle), radius * sin(angle), z)
    end

    index_box = (0.0, Float64(n_radial), 0.0, Float64(n_circumferential), 0.0, Float64(n_axial))
    cartesian = CartesianDiscreteModel(index_box, (n_radial, n_circumferential, n_axial);
                                       map = index_to_physical)
    full_model = Gridap.Geometry.UnstructuredDiscreteModel(cartesian)
    if notch === nothing
        return full_model, 0
    end

    notch_bottom_radius = outer_radius - notch.depth
    kept_cells = Int[]
    for (cell, corners) in enumerate(Gridap.Geometry.get_cell_coordinates(full_model))
        centroid = sum(corners) / length(corners)
        radius = sqrt(centroid[1]^2 + centroid[2]^2)
        angle = atan(centroid[2], centroid[1])
        inside_notch = radius > notch_bottom_radius && angle < notch.half_angle &&
                       notch.z_start < centroid[3] < notch.z_end
        if !inside_notch
            push!(kept_cells, cell)
        end
    end
    n_removed = num_cells(full_model) - length(kept_cells)
    return Gridap.Geometry.DiscreteModelPortion(full_model, kept_cells), n_removed
end

"Stop with a clear message if the grid lines do not pass through the notch edges."
function check_notch_alignment(geometry::PipeGeometry, notch::NotchGeometry,
                               n_radial::Int, n_circumferential::Int,
                               axial_edges::Vector{Float64})
    wall = wall_thickness(geometry)
    radial_index = (wall - notch.depth) / wall * n_radial
    angular_index = notch.half_angle / (pi / n_circumferential)
    if abs(radial_index - round(radial_index)) > 1e-9
        error("no radial node at the notch bottom: n_radial must be a multiple of " *
              "$(round(Int, wall / notch.depth))")
    end
    if abs(angular_index - round(angular_index)) > 1e-9
        error("no node at the notch edge angle: n_circumferential must be a multiple of " *
              "$(round(Int, pi / notch.half_angle))")
    end
    if !(notch.z_start in axial_edges && notch.z_end in axial_edges)
        error("the notch ends are not axial nodes: include them as breakpoints of axial_grid")
    end
    return nothing
end

"""
    polygonal_volume(geometry, n_circumferential; notch = nothing)

Exact volume of the meshed half pipe: the inscribed prism of `n_circumferential` flat
sectors, 0.5 n (ro^2 - ri^2) sin(pi / n) length, minus the notch prism if there is one.
`mesh_volume` must reproduce it to round-off.
"""
function polygonal_volume(geometry::PipeGeometry, n_circumferential::Int; notch = nothing)
    inner_radius = geometry.inner_radius
    outer_radius = geometry.outer_radius
    sector_factor = sin(pi / n_circumferential)
    volume = 0.5 * n_circumferential * (outer_radius^2 - inner_radius^2) * sector_factor *
             geometry.length
    if notch !== nothing
        notch_bottom_radius = outer_radius - notch.depth
        n_notch_sectors = round(Int, notch.half_angle / (pi / n_circumferential))
        volume -= 0.5 * n_notch_sectors * (outer_radius^2 - notch_bottom_radius^2) *
                  sector_factor * (notch.z_end - notch.z_start)
    end
    return volume
end

"Volume of a mesh, integrated with Gridap."
mesh_volume(model) = sum(∫(1.0) * Measure(Triangulation(model), 2))
