using Test
using LinearAlgebra: norm
using GeometryTypes: Point
using Random: seed!

using GeometryReferenceGenericSolver
using GeometryReferenceGenericSolver: find_minimum_distance

r = find_minimum_distance((Point(-1,0), Point(1,0)), (Point(0,-1), Point(0, 1)))
@test r.d_sq < 100eps()
@test norm(r.p1 - [0, 0]) < 100eps()
@test norm(r.p2 - [0, 0]) < 100eps()

xy_plane = 1e6 .* (Point(-1,0,0), Point(1,0,0), Point(0, 1, 0))

seed!(1234)

@testset "projection xy" begin
    # make all points first, to facilitate test debugging and to improve performance measure
    testpoints = [1e5 .* rand(Point{3, Float64}) for i = 1:1000]
    t_start = nothing
    for i = 1:length(testpoints)
        if i == 2
            t_start = time() # don't measure compilation time in first iteration
        end
        p = testpoints[i]
        # project p onto xy_plane
        r = find_minimum_distance(xy_plane, (p,))
        @test isapprox(p[3]^2, r.d_sq)
        @test norm(r.p1[1:2] - p[1:2]) < 1e-6
    end
    # cheap-o performance benchmark
    @show time() - t_start
end

# with ~16 digits of signifance, and some x * x terms, having a `1e8` value is tough.
xy_plane_bad = 1e8 .* (Point(-1,0,0), Point(1,0,0), Point(0, 1, 0))
ps_bad = [
    Point(15099.411397006435, 63329.17689456945, 489.2398433903766),
    Point(56270.48118496456, 53948.45819773468, 34.889982642960504),
    Point(37517.76357686068, 84362.38878880453, 233.50717418757495)
]

@testset "projection xy numerical issue" begin
    for p in ps_bad
        r = find_minimum_distance(xy_plane_bad, (p,))
        @test isapprox(p[3]^2, r.d_sq, rtol=1e-4)
        @test norm(r.p1[1:2] - p[1:2]) < 1e-6
    end
end
