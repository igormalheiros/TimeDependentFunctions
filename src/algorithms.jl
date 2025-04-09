"""
    IGP(data::TimeDependentData, s0::Float64, (i, j)::Tuple{Int, Int}) -> Float64

Computes the **travel time** on a time-dependent arc `(i, j)` when departing at time `s0`, 
based on the piecewise constant speed model proposed by Ichoua et al. (2003).

# Arguments
- `data::TimeDependentData`: A structure containing time-dependent data:
- `s0::Float64`: Departure time.
- `(i, j)::Tuple{Int, Int}`: The origin-destination pair representing the arc.

# Returns
- `Float64`: The **travel time** required to traverse arc `(i, j)` starting at time `s0`.

# Time Complexity
- O(h), where `h` is the number of time intervals in the speed profile.

# Reference
Ichoua, S., Gendreau, M., & Potvin, J. Y. (2003). Vehicle dispatching with time-dependent travel times. *European Journal of Operational Research*.

"""
function IGP(data::TimeDependentData, s0::Float64, (i, j)::Tuple{Int,Int})
    S = data.S
    Z = data.Z
    v = data.v

    s = s0
    h = something(findfirst(i -> S[i] <= s0 < S[i+1], 1:Z), Z)
    c = data.c[i, j]
    k = data.κ[(i, j)]

    s_ = s + (c / v[k, h])

    while s_ > S[h+1]
        c = c - v[k, h] * (S[h+1] - s)
        s = S[h+1]
        h = h + 1
        if h > Z
            break
        end
        s_ = s + (c / v[k, h])
    end
    return s_ - s0
end

"""
    Φ(data::TimeDependentData, s0::Float64, (i, j)::Tuple{Int, Int}) -> Float64

Computes the **arrival time** at the end of arc `(i, j)` when departing at time `s0`,
using a time-dependent travel time model with piecewise constant speeds. The implementation
follows the algorithm showed in Vidal et al. (2020).

# Arguments
- `data::TimeDependentData`: A structure containing time-dependent data:
- `s0::Float64`: Departure time.
- `(i, j)::Tuple{Int, Int}`: The origin-destination pair representing the arc.

# Returns
- `Float64`: The **arrival time** at the destination of arc `(i, j)`.

# Time Complexity
- O(h), where `h` is the number of time intervals in the speed profile.

# References
- Ichoua, S., Gendreau, M., & Potvin, J. Y. (2003). Vehicle dispatching with time-dependent travel times. *European Journal of Operational Research*.
- Vidal, T., Martinelli, R., Pham, T. A., & Hà, M. H. (2020). Arc Routing with Time-Dependent Travel Times and Paths. *Transportation Science*.

"""
function Φ(data::TimeDependentData, s0::Float64, (i, j)::Tuple{Int,Int})
    S = data.S
    Z = data.Z
    v = data.v

    s = s0
    h = something(findfirst(x -> x > s0, S), Z + 1) - 1
    c = data.c[i, j]
    k = data.κ[(i, j)]

    s_ = s + (c / v[k, h])

    while s_ > S[h+1]
        c = c - v[k, h] * (S[h+1] - s)
        s = S[h+1]
        h = h + 1
        if h > Z
            break
        end
        s_ = s + (c / v[k, h])
    end
    return s_
end

"""
    Φ_inv(data::TimeDependentData, s0::Float64, (i, j)::Tuple{Int, Int}) -> Float64

Computes the **departure time** from the origin of arc `(i, j)` such that arrival occurs 
at time `s0`, using a time-dependent model with piecewise constant speeds. The implementation
follows the algorithm showed in Vidal et al. (2020).

# Arguments
- `data::TimeDependentData`: A structure containing time-dependent data:
- `s0::Float64`: Departure time.
- `(i, j)::Tuple{Int, Int}`: The origin-destination pair representing the arc.

# Returns
- `Float64`: The **departure time** such that arrival at the of arc `(i, j)` occurs at time `s0`.

# Time Complexity
- O(h), where `h` is the number of time intervals in the speed profile.

# References
- Ichoua, S., Gendreau, M., & Potvin, J. Y. (2003). Vehicle dispatching with time-dependent travel times. *European Journal of Operational Research*.
- Vidal, T., Martinelli, R., Pham, T. A., & Hà, M. H. (2020). Arc Routing with Time-Dependent Travel Times and Paths. *Transportation Science*.

"""
function Φ_inv(data::TimeDependentData, s0::Float64, (i, j)::Tuple{Int,Int})
    S = data.S
    v = data.v

    s = s0
    h = something(findlast(x -> x < s0, S), 1)
    c = data.c[i, j]
    k = data.κ[(i, j)]

    s_ = s - (c / v[k, h])

    while s_ < S[h]
        c = c - (v[k, h] * (s - S[h]))
        s = S[h]
        h = h - 1
        if h < 1
            break
        end
        s_ = s - (c / v[k, h])
    end
    return s_
end

"""
    travel_time_breakpoints(data::TimeDependentData, (i, j)::Tuple{Int, Int}) -> Vector{Tuple{Float64, Float64}}

Computes the **breakpoints** of the piecewise linear **time-dependent travel time function** for arc `(i, j)`,
based on the time-dependent speed model proposed by Ichoua et al. (2003).

This algorithm is adapted from the method proposed by Vidal et al. (2020) for computing breakpoints of
the **arrival time function**, and modified here to compute breakpoints of the **travel time function** instead.

# Arguments
- `data::TimeDependentData`: A structure containing time-dependent network data, including:
- `(i, j)::Tuple{Int, Int}`: The origin-destination pair representing the arc.

# Returns
- `Vector{Tuple{Float64, Float64}}`: A list of breakpoints `(s, t)` where `s` is the departure time
  and `t` is the corresponding travel time on arc `(i, j)`.

# Time Complexity
- `O(h²)`, where `h` is the number of time intervals in the speed profile of arc `(i, j)`.

# References
- Ichoua, S., Gendreau, M., & Potvin, J. Y. (2003). Vehicle dispatching with time-dependent travel times. *European Journal of Operational Research*, 144(2), 379–396.
- Vidal, T., Martinelli, R., Pham, T. A., & Hà, M. H. (2020). Arc Routing with Time-Dependent Travel Times and Paths. *Transportation Science*.

"""
function travel_time_breakpoints(data::TimeDependentData, (i, j)::Tuple{Int,Int})
    @show "000000000"
    S = data.S
    Z = data.Z
    v = data.v

    h_r = findlast(x -> x <= Φ_inv(data, S[end], (i, j)), S)
    h_l = findfirst(x -> x >= Φ(data, 0.0, (i, j)), S)

    B = Tuple{Float64,Float64}[]

    push!(B, (0.0, Φ(data, 0.0, (i, j))))
    for h = 1:h_r
        @show S[h]
        push!(B, (S[h], Φ(data, S[h], (i, j))))
    end
    for h = h_l:Z
        push!(B, (Φ_inv(data, S[h], (i, j)), S[h]))
    end
    push!(B, (Φ_inv(data, S[end], (i, j)), S[end]))

    unique!(B)
    sort!(B, by = x -> x[1])
    B = [(x[1], x[2] - x[1]) for x in B]
    return B
end

function linear_piecewise_affine_t(s::Float64, segments::Vector{AffineSegment})
    for i = 1:length(segments)
        # If not the last segment, check if s is before the next segment's start.
        if i < length(segments)
            if s < segments[i+1].x_min
                return segments[i].m * s + segments[i].y
            end
        else
            # Last segment: assume s belongs here if s >= x_min.
            if s >= segments[i].x_min
                return segments[i].m * s + segments[i].y
            end
        end
    end
    error("s = $s is outside the domain of the piecewise function.")
end

function bs_piecewise_affine_t(s::Float64, segments::Vector{AffineSegment})
    # Extract x_min values for binary search
    x_mins = getfield.(segments, :x_min)  # Corrected from getindex to getfield

    # Find the rightmost segment where x_min <= s using binary search
    i = searchsortedlast(x_mins, s)

    if i == 0 || i > length(segments)
        error("s = $s is outside the domain of the piecewise function.")
    end

    # Evaluate the function at s
    return segments[i].m * s + segments[i].y
end
