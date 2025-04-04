function IGP(data::TimeDependentData, s0::Float64, (i, j)::Tuple{Int,Int})
    S = data.S
    Z = data.Z
    v = data.v
    @show "XXXXXXXXX"
    s = s0
    h = get_speed_interval(S, Z, s0)
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

function travel_time_breakpoints(data::TimeDependentData, (i, j)::Tuple{Int,Int})
    S = data.S
    Z = data.Z
    H = data.H
    v = data.v
    k = data.κ[(i, j)]

    B = Tuple{Float64,Float64}[]
    s = 0.0
    L = IGP(data, s, (i, j))
    push!(B, (s, L))

    while s < S[end-1]
        Sl = next_time_threshold(S, s)
        Sr = next_time_threshold(S, s + L)

        vl = v[k, get_speed_interval(S, Z, s)]
        vr = v[k, get_speed_interval(S, Z, s + L)]

        rho_l = Sl - s
        rho_r = (vl / vr) * (Sr - s - L)
        rho = min(rho_l, rho_r)
        L = L + (vl / vr - 1) * rho

        s += rho
        t = IGP(data, s, (i, j))
        push!(B, (s, t))
    end
    t = IGP(data, s, (i, j))
    push!(B, (H, t))
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
