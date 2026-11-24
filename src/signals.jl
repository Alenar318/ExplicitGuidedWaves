# Signal comparison: resampling, normalised RMS error and delay.
#
# The error measure is the one of the paper: the RMS of the difference in % of the
# range (max - min) of the reference signal,
#
#     NRMSE = 100 sqrt( mean( (u - u_ref)^2 ) ) / (max u_ref - min u_ref)

"""
    resample_on_grid(times, values, new_times)

Linear interpolation of `values` (given at increasing `times`) at `new_times`, which
must lie inside [times[1], times[end]].
"""
function resample_on_grid(times::Vector{Float64}, values::Vector{Float64},
                          new_times::AbstractVector{Float64})
    if first(new_times) < times[1] || last(new_times) > times[end] * (1 + 1e-12)
        error("resample_on_grid: the new grid goes outside the signal")
    end
    resampled = similar(values, length(new_times))
    interval = 1
    for (i, t) in enumerate(new_times)
        while interval < length(times) - 1 && times[interval + 1] < t
            interval += 1
        end
        weight = (t - times[interval]) / (times[interval + 1] - times[interval])
        resampled[i] = (1 - weight) * values[interval] + weight * values[interval + 1]
    end
    return resampled
end

signal_range(signal::AbstractVector{Float64}) = maximum(signal) - minimum(signal)

"""
    nrmse_percent(signal, reference; scale = signal_range(reference))

RMS of `signal - reference` in % of `scale`. Pass the range of another signal as
`scale` to normalise several components by the same quantity (the paper normalises
u_r and u_x by the range of u_theta).
"""
function nrmse_percent(signal::AbstractVector{Float64}, reference::AbstractVector{Float64};
                       scale::Float64 = signal_range(reference))
    if length(signal) != length(reference)
        error("nrmse_percent: signals of different length")
    end
    return 100 * sqrt(sum(abs2, signal .- reference) / length(reference)) / scale
end

"""
    cross_correlation_delay(signal, reference, sample_interval)

Delay of `signal` relative to `reference` (s; positive: `signal` arrives later), from
the peak of their cross-correlation, refined to a fraction of a sample with a parabola
through the peak and its two neighbours. Both signals are sampled at the same instants.
"""
function cross_correlation_delay(signal::AbstractVector{Float64},
                                 reference::AbstractVector{Float64},
                                 sample_interval::Float64)
    a = signal .- sum(signal) / length(signal)
    b = reference .- sum(reference) / length(reference)
    lags = -(length(b) - 1):(length(a) - 1)
    correlation = zeros(length(lags))
    for (k, lag) in enumerate(lags)
        total = 0.0
        for i in eachindex(b)
            j = i + lag
            if 1 <= j <= length(a)
                total += a[j] * b[i]
            end
        end
        correlation[k] = total
    end
    peak = argmax(correlation)
    shift = Float64(lags[peak])
    if 1 < peak < length(correlation)
        left = correlation[peak - 1]
        centre = correlation[peak]
        right = correlation[peak + 1]
        shift += 0.5 * (left - right) / (left - 2 * centre + right)
    end
    return shift * sample_interval
end
