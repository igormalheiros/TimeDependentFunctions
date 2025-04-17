using Test, TimeDependentFunctions

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

    @test isapprox(Φ(data1, 0.0, (1, 2)), 2.5; atol = 1e-1)
    @test isapprox(Φ(data1, 0.0, (1, 3)), 0.5; atol = 1e-1)
    @test isapprox(Φ(data1, 0.0, (2, 3)), 5.0; atol = 1e-1)

    @test isapprox(Φ_inv(data1, 2.5, (1, 2)), 0.0; atol = 1e-1)
    @test isapprox(Φ_inv(data1, 0.5, (1, 3)), 0.0; atol = 1e-1)
    @test isapprox(Φ_inv(data1, 5.0, (2, 3)), 0.0; atol = 1e-1)

    @test isapprox(Φ(data1, Φ_inv(data1, 2.5, (1, 2)), (1, 2)), 2.5; atol = 1e-1)
    @test isapprox(Φ(data1, Φ_inv(data1, 0.5, (1, 3)), (1, 3)), 0.5; atol = 1e-1)
    @test isapprox(Φ(data1, Φ_inv(data1, 5.0, (2, 3)), (2, 3)), 5.0; atol = 1e-1)

    @test isapprox(Φ(data1, 3.7, (1, 2)), 5.6; atol = 1e-1)
    @test isapprox(Φ(data1, 3.7, (1, 3)), 4.2; atol = 1e-1)
    @test isapprox(Φ(data1, 3.7, (2, 3)), 5.74; atol = 1e-1)

    @test isapprox(Φ_inv(data1, 5.6, (1, 2)), 3.7; atol = 1e-1)
    @test isapprox(Φ_inv(data1, 4.2, (1, 3)), 3.7; atol = 1e-1)
    @test isapprox(Φ_inv(data1, 5.74, (2, 3)), 3.7; atol = 1e-1)

    @test isapprox(Φ(data1, Φ_inv(data1, 5.6, (1, 2)), (1, 2)), 5.6; atol = 1e-1)
    @test isapprox(Φ(data1, Φ_inv(data1, 4.2, (1, 3)), (1, 3)), 4.2; atol = 1e-1)
    @test isapprox(Φ(data1, Φ_inv(data1, 5.74, (2, 3)), (2, 3)), 5.74; atol = 1e-1)

    B_12 = travel_time_breakpoints(data1, (1, 2))
    B_13 = travel_time_breakpoints(data1, (1, 3))
    B_23 = travel_time_breakpoints(data1, (2, 3))

    @test B_12 ==
          [(0.0, 2.5), (2.5, 2.5), (5.0, 1.25), (8.75, 1.25), (10.0, 5.0), (15.0, 5.0)]
    @test B_13 ==
          [(0.0, 0.5), (4.5, 0.5), (5.0, 0.25), (9.75, 0.25), (10.0, 1.0), (19.0, 1.0)]

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
                Φ_t(data1, s, arc);
                atol = 1e-2,
            )
            @test isapprox(
                bs_piecewise_affine_t(s, affine_B),
                Φ_t(data1, s, arc);
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
                Φ_t(data2, s, arc);
                atol = 1e-2,
            )
            @test isapprox(
                bs_piecewise_affine_t(s, affine_B),
                Φ_t(data2, s, arc);
                atol = 1e-2,
            )
            s += step
        end
    end

    c3 = Matrix{Float64}(undef, 1, 1)
    c3[1, 1] = 380
    k3 = Dict((1, 1) => 1)
    arc = (i, j) = (1, 1)
    S3 = [0.0, 480.0, 720.0, 1680.0, 1920.0, 2400.0]
    v3 = [1.16667 0.666667 1.33333 0.833333 1.0]

    data3 = TimeDependentData(S3, k3, v3, c3)
    B = travel_time_breakpoints(data3, arc)
    affine_B = build_segments(B)

    s = 0.0
    step = 1
    while s < S3[end-1] + 2
        @test isapprox(
            linear_piecewise_affine_t(s, affine_B),
            Φ_t(data3, s, arc);
            atol = 1e-2,
        )
        @test isapprox(bs_piecewise_affine_t(s, affine_B), Φ_t(data3, s, arc); atol = 1e-2)
        s += step
    end
end
