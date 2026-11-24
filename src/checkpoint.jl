# Restart files for long runs. A run that is interrupted (shutdown, power cut) and
# started again with the same checkpoint file resumes where it stopped and gives
# bitwise the same result as an uninterrupted run (test/runtests.jl).
#
# Serialization (.jls) depends on the Julia version and on this code: checkpoints are
# disposable caches, not results. Results are written as text (results_io.jl).

"Bumped whenever the content of a checkpoint changes."
const CHECKPOINT_VERSION = 2

"""
    save_checkpoint(path, content)

Write `content` atomically: first to `path * ".tmp"`, then rename over `path`. A power
cut during the write leaves the previous checkpoint intact.
"""
function save_checkpoint(path::AbstractString, content)
    folder = dirname(path)
    if !isempty(folder)
        mkpath(folder)
    end
    temporary = path * ".tmp"
    open(temporary, "w") do io
        serialize(io, content)
    end
    mv(temporary, path; force = true)
    return nothing
end

"""
    load_checkpoint(path)

The saved content, or `nothing` if there is no checkpoint. Falls back to the `.tmp`
file in case the power went off between removing the old file and renaming the new
one; a truncated `.tmp` (the write itself was cut) is ignored.
"""
function load_checkpoint(path::AbstractString)
    if isfile(path)
        return open(deserialize, path)
    end
    temporary = path * ".tmp"
    if isfile(temporary)
        try
            return open(deserialize, temporary)
        catch exception
            if exception isa EOFError
                @warn "ignoring a truncated temporary checkpoint" temporary
            else
                rethrow()
            end
        end
    end
    return nothing
end
