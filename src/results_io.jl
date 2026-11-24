# Results on disk, as plain text so they outlive this code and this Julia version.
#
# Signal files use the column layout of the COMSOL and BEM exports of the paper,
#     time (s)   u_r   u_theta   u_z
# with header lines starting with '%', so the three sources are read the same way.

"Short git commit of the code, or \"unknown\" outside a git repository."
function git_commit()
    try
        # stderr to devnull: outside a repository git prints "fatal: not a git repository"
        return readchomp(pipeline(`git -C $(@__DIR__) rev-parse --short HEAD`; stderr = devnull))
    catch exception
        if exception isa ProcessFailedException || exception isa Base.IOError
            return "unknown"
        end
        rethrow()
    end
end

"Lines describing how a result was produced: date, code version, Julia, threads."
function provenance_lines()
    return ["date: $(Dates.format(now(), "yyyy-mm-dd HH:MM"))",
            "code: ExplicitGuidedWaves, git commit $(git_commit())",
            "julia: $(VERSION), threads: $(Threads.nthreads())"]
end

"""
    write_signal_file(path, times, components, description)

Write one sensor: columns time (s), radial, circumferential, axial. `components` is a
NamedTuple (radial, circumferential, axial); `description` is a vector of lines put in
the header after the provenance.
"""
function write_signal_file(path::AbstractString, times::AbstractVector{Float64},
                           components::NamedTuple, description::Vector{String})
    mkpath(dirname(path))
    open(path, "w") do io
        for line in vcat(description, provenance_lines())
            println(io, "% ", line)
        end
        println(io, "% columns: time (s), u_r, u_theta, u_z")
        writedlm(io, hcat(times, components.radial, components.circumferential, components.axial), ' ')
    end
    return path
end

"""
    read_signal_file(path)

Read a signal file of this project, of COMSOL or of the BEM codes (header lines start
with '%'). Returns (time, radial, circumferential, axial).
"""
function read_signal_file(path::AbstractString)
    if !isfile(path)
        error("signal file not found: $path")
    end
    table = readdlm(path; comments = true, comment_char = '%')
    if size(table, 2) != 4
        error("$path: expected 4 columns (time, u_r, u_theta, u_z), found $(size(table, 2))")
    end
    return (time = Float64.(table[:, 1]), radial = Float64.(table[:, 2]),
            circumferential = Float64.(table[:, 3]), axial = Float64.(table[:, 4]))
end

# The tables are written with `join` rather than CSV.jl, so that the solver has no
# table dependency; the notebooks read them with CSV.jl and DataFrames.

"""
    append_csv_row(path, header, values)

Append one row to a CSV table, writing `header` first if the file is new. Stops if the
file exists with a different header.
"""
function append_csv_row(path::AbstractString, header::Vector{String}, values::Vector)
    if length(header) != length(values)
        error("append_csv_row: $(length(header)) column names for $(length(values)) values")
    end
    mkpath(dirname(path))
    if isfile(path)
        existing_header = split(readline(path), ',')
        if existing_header != header
            error("$path has columns $(existing_header), not $(header)")
        end
    else
        open(path, "w") do io
            println(io, join(header, ','))
        end
    end
    open(path, "a") do io
        println(io, join(string.(values), ','))
    end
    return path
end

"""
    read_csv_table(path)

Read a table written by `append_csv_row`. Returns a Dict column name => vector.
"""
function read_csv_table(path::AbstractString)
    if !isfile(path)
        error("table not found: $path")
    end
    values, header = readdlm(path, ','; header = true)
    return Dict(String(name) => values[:, column] for (column, name) in enumerate(vec(header)))
end
