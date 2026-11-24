# Durbin inversion of the BEM results of the notched pipe (hybrid and full 3D) to
# plain text, so they can be compared with Gridap and COMSOL outside the BEM code.
#
# Uses the BEM project's own code (src/Hibrido.jl) and its default environment, so
# run it WITHOUT --project:
#   julia scripts/export_bem_notch.jl
#
# Same inversion parameters as herramientas/figuras_entalla.jl. Writes, per sensor,
# results/notch/bem_<model>_x<x>_th<deg>.txt with columns t, u_r, u_theta, u_x,
# sampled every 1 us like the implicit COMSOL files.

using Printf, Serialization

const BEM = normpath(joinpath(@__DIR__, "..", "..", "BEM_HIBRIDO_AXI_3D_JULIA"))
include(joinpath(BEM, "src", "Hibrido.jl"))
const BA = Hibrido.BEMAxi

const FREQUENCY, N_CYCLES, T_MAX = 40.0e3, 5, 500.0e-6
const XS = (0.75, 0.25)
const THS = (0.0, π / 2, π)
const DIRS = (:θ, :x, :r)
const SENS = [(x, θ, d) for x in XS for θ in THS for d in DIRS]   # row order of R.U

T0, ND = BA.sugerir_laplace(T_MAX, BA.f_max_tone_burst(FREQUENCY, N_CYCLES))
SK, AD = BA.frecuencias_laplace(T0, ND)
burst_s = [BA.tone_burst_s(FREQUENCY, N_CYCLES, s) for s in SK]
times = collect(0.0:1.0e-6:T_MAX)
out_dir = joinpath(@__DIR__, "..", "results", "notch")
mkpath(out_dir)
@printf("N = %d values of s, T0 = %.3e s\n", ND, T0)

for model in ("hibrido_entalla", "completo_entalla")
    file = joinpath(BEM, "resultados", "entalla_L1_40kHz", "$(model)_N$(ND).jls")
    if !isfile(file)
        println("not found: ", file)
        continue
    end
    R = deserialize(file)
    @printf("%s: %.1f min in the BEM run\n", model, R.tiempo / 60)
    for x in XS, θ in THS
        signal(d) = real.(BA.durbin(times, R.U[findfirst(==((x, θ, d)), SENS), :] .* burst_s, T0, AD))
        u_theta, u_x, u_r = signal(:θ), signal(:x), signal(:r)
        name = @sprintf("bem_%s_x%g_th%d.txt", model, x, round(Int, rad2deg(θ)))
        open(joinpath(out_dir, name), "w") do io
            println(io, "% BEM ($model), Durbin inversion with N = $ND; columns t, u_r, u_theta, u_x (m)")
            for i in eachindex(times)
                @printf(io, "%.6e  %.10e  %.10e  %.10e\n", times[i], u_r[i], u_theta[i], u_x[i])
            end
        end
    end
    println("   written to results/notch/bem_$(model)_*.txt")
end
