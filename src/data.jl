"""
    TimeDependentData

A struct that contains the data needed for the time dependent functions.
The notation follows tha same as the paper.

# Fields
- `S::Vector{Float64}`: The set of time intervals.
- `Z::Int`: The number of time intervals.
- `H::Float64`: The time horizon.
- `κ::Dict{Tuple{Int, Int}, Int}`: The speed profile by arc.
- `K::Int`: The number of speed profiles.
- `v::Matrix{Float64}`: The speed by each profile and time interval.
- `c::Matrix{Float64}`: The distance by arc.
"""

struct TimeDependentData
    S::Vector{Float64}
    Z::Int
    H::Float64
    κ::Dict{Tuple{Int,Int},Int}
    K::Int
    v::Matrix{Float64}
    c::Matrix{Float64}
end

function TimeDependentData(
    S::Vector{Float64},
    κ::Dict{Tuple{Int,Int},Int},
    v::Matrix{Float64},
    c::Matrix{Float64},
)
    Z = length(S) - 1
    H = S[end]
    K = size(v, 1)
    return TimeDependentData(S, Z, H, κ, K, v, c)
end
