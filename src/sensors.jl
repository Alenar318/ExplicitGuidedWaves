# Sensors on the outer surface. Each sensor sits on a mesh node, so the recorded
# signal has no interpolation error; the three Cartesian components are recorded and
# projected afterwards on (radial, circumferential, axial).

"Point on the outer surface, at axial position z and angle theta from the x axis."
struct SurfaceSensor
    name::String
    axial_position::Float64   # m
    angle_deg::Float64        # deg; the cut plane y = 0 is at 0 and 180
end

function sensor_point(sensor::SurfaceSensor, geometry::PipeGeometry)
    radius = geometry.outer_radius
    angle = sensor.angle_deg
    return VectorValue(radius * cosd(angle), radius * sind(angle), sensor.axial_position)
end

"""
    node_component_dofs(space, points; tolerance = 1.0e-7)

For each point, the degrees of freedom of the x, y and z components of the mesh node
located there, as a (number of points) x 3 matrix. Stops if a point is not a node.
"""
function node_component_dofs(space, points; tolerance::Float64 = 1.0e-7)
    x_at_dof = node_coordinate_at_dofs(space, 1)
    y_at_dof = node_coordinate_at_dofs(space, 2)
    z_at_dof = node_coordinate_at_dofs(space, 3)
    masks = [component_mask(space, direction) for direction in 1:3]

    dofs = zeros(Int, length(points), 3)
    for (point_index, point) in enumerate(points)
        at_node = (abs.(x_at_dof .- point[1]) .< tolerance) .&
                  (abs.(y_at_dof .- point[2]) .< tolerance) .&
                  (abs.(z_at_dof .- point[3]) .< tolerance)
        for direction in 1:3
            hits = findall(masks[direction] .& at_node)
            if length(hits) != 1
                error("point $(point) is not a mesh node ($(length(hits)) matches for " *
                      "component $direction): use an even n_circumferential and an axial " *
                      "grid that contains the sensor positions")
            end
            dofs[point_index, direction] = hits[1]
        end
    end
    return dofs
end

"""
    recorded_dofs(space, geometry, sensors)

Degrees of freedom to record for `sensors`, sensor by sensor:
[x1, y1, z1, x2, y2, z2, ...]. Column 3 (k - 1) + c of the history is component c of
sensor k.
"""
function recorded_dofs(space, geometry::PipeGeometry, sensors::Vector{SurfaceSensor})
    points = [sensor_point(sensor, geometry) for sensor in sensors]
    return vec(permutedims(node_component_dofs(space, points)))
end

"""
    cylindrical_components(history, sensors)

Project the recorded Cartesian components on the cylindrical basis at each sensor:

    radial          =  x cos(theta) + y sin(theta)
    circumferential = -x sin(theta) + y cos(theta)
    axial           =  z

Returns a Dict sensor name => (radial, circumferential, axial) vectors.
"""
function cylindrical_components(history::Matrix{Float64}, sensors::Vector{SurfaceSensor})
    components = Dict{String, NamedTuple}()
    for (sensor_index, sensor) in enumerate(sensors)
        first_column = 3 * (sensor_index - 1)
        x_component = history[:, first_column + 1]
        y_component = history[:, first_column + 2]
        z_component = history[:, first_column + 3]
        angle = sensor.angle_deg
        components[sensor.name] = (
            radial = x_component .* cosd(angle) .+ y_component .* sind(angle),
            circumferential = .-x_component .* sind(angle) .+ y_component .* cosd(angle),
            axial = z_component,
        )
    end
    return components
end
