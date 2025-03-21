"""
    AffineSegment

A struct that represents an affine segment.

# Fields
- `m::Float64`: The slope of the segment.
- `y::Float64`: The y-intercept of the segment.
- `x_min::Float64`: The minimum x value of the segment.
"""
struct AffineSegment
    m::Float64
    y::Float64
    x_min::Float64
end

function build_segments(B::Vector{Tuple{Float64,Float64}})
    segments = AffineSegment[]
    for i = 1:length(B)-1
        (x1, y1) = B[i]
        (x2, y2) = B[i+1]
        # Calculate slope and intercept for the segment
        m = (y2 - y1) / (x2 - x1)
        y = y1 - m * x1
        push!(segments, AffineSegment(m, y, x1))
    end
    return segments
end
