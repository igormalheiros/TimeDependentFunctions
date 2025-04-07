function get_speed_interval(S::Vector{Float64}, Z::Int, s::Float64)
    if s > S[end]
        error("s = $s is larger than the time horizon.")
    end
    return something(findfirst(i -> S[i] <= s < S[i+1], 1:Z), Z)
end

function next_time_threshold(S::Vector{Float64}, s::Float64)
    if s > S[end]
        error("s = $s is larger than the time horizon.")
    end
    if s == S[end]
        return S[end]
    end
    return minimum([s_ for s_ in S if s_ > s])
end
