# TimeDependentFunctions.jl

A Julia package of auxiliary functions for computing **time-dependent travel times** on arcs whose speed varies over the course of the day. It provides ready-to-use, tested implementations for researchers and practitioners working on time-dependent routing and scheduling problems, including:

- The time-dependent travel time given the departure time and speed model (Ichoua et al., 2003).
- The time-dependent arrival time given the departure time and speed model (Vidal et al., 2020).
- The departure time given the arrival time and speed model (Ichoua et al., 2003 / Vidal et al., 2020).
- The closed-form travel time breakpoint representation given the speed model (Vidal et al., 2020 / Malheiros et al., 2026).
- The travel time given the departure time and breakpoint representation in linear time (Malheiros et al., 2026).
- The travel time given the departure time and breakpoint representation in logarithmic time (Malheiros et al., 2026).

## Installation

The package is not yet registered. Add it directly from the local path or a Git URL:

```julia
using Pkg
Pkg.develop(path = "path/to/TimeDependentFunctions")
# or
Pkg.add(url = "https://github.com/igormalheiros/TimeDependentFunctions.jl")
```

## Usage

```julia
using TimeDependentFunctions

# Shared time-of-day breakpoints (e.g. minutes from midnight)
S = [0.0, 5.0, 10.0, 20.0]

# Each arc (i, j) is assigned a speed profile index κ[(i, j)]
κ = Dict((1, 2) => 1, (1, 3) => 1, (2, 3) => 2)

# v[k, h] = speed of profile k during time interval h
v = [
    4.0 8.0 2.0
    1.0 5.0 6.0
]

# c[i, j] = distance of arc (i, j)
c = [
    0.0 10.0 2.0
    10.0 0.0 5.0
    2.0 5.0 0.0
]

data = TimeDependentData(S, κ, v, c)

# Arrival time when leaving node 1 towards node 2 at time 0.0
arrival = Φ(data, 0.0, (1, 2))

# Travel time for the same trip
duration = Φ_t(data, 0.0, (1, 2))

# Departure time needed to arrive at `arrival`
departure = Φ_inv(data, arrival, (1, 2))  # ≈ 0.0

# Piecewise-linear approximation of the travel time function for arc (1, 2)
breakpoints = travel_time_breakpoints(data, (1, 2))
segments = build_segments(breakpoints)

# Evaluate the approximation at any departure time
linear_piecewise_affine_t(3.7, segments)   # linear scan
bs_piecewise_affine_t(3.7, segments)       # binary search, faster for many segments
```

## API overview

| Function | Description | Complexity |
|---|---|---|
| `Φ(data, s0, (i, j))` | Arrival time departing `(i, j)` at `s0` | O(h) |
| `Φ_t(data, s0, (i, j))` | Travel time departing `(i, j)` at `s0` | O(h) |
| `Φ_inv(data, s0, (i, j))` | Departure time to arrive at `s0` | O(h) |
| `travel_time_breakpoints(data, (i, j))` | Breakpoints of the piecewise-linear travel time function | O(h²) |
| `build_segments(breakpoints)` | Converts breakpoints into `AffineSegment`s | O(h) |
| `linear_piecewise_affine_t(s, segments)` | Evaluates the piecewise-linear function via linear scan | O(n) |
| `bs_piecewise_affine_t(s, segments)` | Evaluates the piecewise-linear function via binary search | O(log n) |

where `h` is the number of time intervals in the speed profile of the arc being queried, and `n` is the number of breakpoint segments.

## Running the tests

```julia
using Pkg
Pkg.test("TimeDependentFunctions")
```

## References

- Ichoua, S., Gendreau, M., & Potvin, J. Y. (2003). [Vehicle dispatching with time-dependent travel times](https://doi.org/10.1016/S0377-2217(02)00147-9). *European Journal of Operational Research*, 144(2), 379–396.
- Vidal, T., Martinelli, R., Pham, T. A., & Hà, M. H. (2020). [Arc Routing with Time-Dependent Travel Times and Paths](https://doi.org/10.1287/trsc.2020.1035). *Transportation Science*.
