module TimeDependentFunctions

include("data.jl")
include("utils.jl")
include("piecewise_affine.jl")
include("algorithms.jl")

export TimeDependentData,
    AffineSegment,
    build_segments,
    IGP,
    travel_time_breakpoints,
    linear_piecewise_affine_t,
    bs_piecewise_affine_t

end # module TimeDependentFunctions
