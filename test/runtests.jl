using Test, TimeDependentFunctions

import TimeDependentFunctions: get_speed_interval, next_time_threshold

@testset "Util Functions" begin
    S1 = [0.0, 2.5, 4.0, 8.0]
    Z1 = length(S1) - 1

    @test get_speed_interval(S1, Z1, 0.0) == 1
    @test get_speed_interval(S1, Z1, 2.4) == 1
    @test get_speed_interval(S1, Z1, 2.5) == 2
    @test get_speed_interval(S1, Z1, 3.0) == 2
    @test get_speed_interval(S1, Z1, 4.0) == 3
    @test get_speed_interval(S1, Z1, 7.0) == 3
    @test get_speed_interval(S1, Z1, 8.0) == 3
    @test_throws ErrorException get_speed_interval(S1, Z1, 100.0)

    @test next_time_threshold(S1, 0.0) == 2.5
    @test next_time_threshold(S1, 2.5) == 4.0
    @test next_time_threshold(S1, 4.0) == 8.0
    @test next_time_threshold(S1, 5.0) == 8.0
    @test next_time_threshold(S1, 8.0) == 8.0
    @test_throws ErrorException next_time_threshold(S1, 8.1)
    @test_throws ErrorException next_time_threshold(S1, 100.0)

    S2 = [0.0, 5.0]
    Z2 = length(S2) - 1

    @test get_speed_interval(S2, Z2, 0.0) == 1
    @test get_speed_interval(S2, Z2, 4.9) == 1
    @test get_speed_interval(S2, Z2, 5.0) == 1
    @test_throws ErrorException get_speed_interval(S2, Z2, 9.9)

    @test next_time_threshold(S2, 0.0) == 5.0
    @test next_time_threshold(S2, 4.9) == 5.0
    @test next_time_threshold(S2, 5.0) == 5.0
    @test_throws ErrorException next_time_threshold(S2, 5.1)
    @test_throws ErrorException next_time_threshold(S2, 100.0)
end

@testset "TimeDependent Functions" begin
    S1 = [0.0, 5.0, 10.0, 20.0]
    κ1 = Dict((1, 2) => 1, (1, 3) => 1, (2, 3) => 2)
    v1 = [
        4.0 8.0 2.0
        1.0 5.0 6.0
    ]
    c1 = [
        0.0 10.0 2.0
        10.0 0.0 5.0
        2.0 5.0 0.0
    ]
    data1 = TimeDependentData(S1, κ1, v1, c1)

    @test isapprox(IGP(data1, 0.0, (1, 2)), 2.5; atol = 1e-1)
    @test isapprox(IGP(data1, 0.0, (1, 3)), 0.5; atol = 1e-1)
    @test isapprox(IGP(data1, 0.0, (2, 3)), 5.0; atol = 1e-1)

    @test isapprox(IGP(data1, 3.7, (1, 2)), 1.9; atol = 1e-1)
    @test isapprox(IGP(data1, 3.7, (1, 3)), 0.5; atol = 1e-1)
    @test isapprox(IGP(data1, 3.7, (2, 3)), 2.04; atol = 1e-1)

    B_12 = travel_time_breakpoints(data1, (1, 2))
    B_13 = travel_time_breakpoints(data1, (1, 3))
    B_23 = travel_time_breakpoints(data1, (2, 3))

    @test B_12 ==
          [(0.0, 2.5), (2.5, 2.5), (5.0, 1.25), (8.75, 1.25), (10.0, 5.0), (20.0, 5.0)]
    @test B_13 ==
          [(0.0, 0.5), (4.5, 0.5), (5.0, 0.25), (9.75, 0.25), (10.0, 1.0), (20.0, 1.0)]

    affine_B_12 = build_segments(B_12)
    affine_B_13 = build_segments(B_13)
    affine_B_23 = build_segments(B_23)

    @test isapprox(linear_piecewise_affine_t(0.0, affine_B_12), 2.5; atol = 1e-1)
    @test isapprox(linear_piecewise_affine_t(3.7, affine_B_12), 1.9; atol = 1e-1)
    @test isapprox(bs_piecewise_affine_t(0.0, affine_B_12), 2.5; atol = 1e-1)
    @test isapprox(bs_piecewise_affine_t(3.7, affine_B_12), 1.9; atol = 1e-1)

    for arc in keys(κ1)
        B = travel_time_breakpoints(data1, arc)
        affine_B = build_segments(B)
        s = 0.0
        step = 0.1
        while s < S1[end-1] + 2
            @test isapprox(
                linear_piecewise_affine_t(s, affine_B),
                IGP(data1, s, arc);
                atol = 1e-2,
            )
            @test isapprox(
                bs_piecewise_affine_t(s, affine_B),
                IGP(data1, s, arc);
                atol = 1e-2,
            )
            s += step
        end
    end

    S2 = [0.0, 3.0, 10.0]
    κ2 = Dict((1, 2) => 1, (1, 3) => 2, (2, 3) => 3)
    v2 = [
        12.0 7.5
        3.0 6.3
        9.1 4.2
    ]
    c2 = [
        0.0 4.0 7.0
        4.0 0.0 11.2
        7.0 11.2 0.0
    ]
    data2 = TimeDependentData(S2, κ2, v2, c2)
    for arc in keys(κ2)
        B = travel_time_breakpoints(data2, arc)
        affine_B = build_segments(B)
        s = 0.0
        step = 0.1
        while s < S2[end-1] + 2
            @test isapprox(
                linear_piecewise_affine_t(s, affine_B),
                IGP(data2, s, arc);
                atol = 1e-2,
            )
            @test isapprox(
                bs_piecewise_affine_t(s, affine_B),
                IGP(data2, s, arc);
                atol = 1e-2,
            )
            s += step
        end
    end
end
