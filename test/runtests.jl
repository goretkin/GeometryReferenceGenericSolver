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

xy_plane = 1e8 .* (Point(-1,0,0), Point(1,0,0), Point(0, 1, 0))

seed!(1234)
@testset "project xy" begin
    # make all points first, to facilitate test debugging
    testpoints = [1e5 .* rand(Point{3, Float64}) for i = 1:1000]

    for i = 970:1000
        @show i
        p = testpoints[i]
        @show p
        # project p onto xy_plane
        r = find_minimum_distance(xy_plane, (p,))
        @show r
        @test isapprox(p[3]^2, r.d_sq)
        @test norm(r.p1[1:2] - p[1:2]) < 1e-6
    end
end

# Somep oints caused an MathOptInterface.ITERATION_LIMIT with
#=     model = Model(OSQP.Optimizer(
        max_iter=400000,
        warm_start=false,   # for deterministic results
        polish=true,
        eps_prim_inf=1e-14,
        eps_dual_inf=1e-14,
        eps_abs=1e-14,
        eps_rel=1e-6)
    )
    MOI.set(model.optimizer, OSQPSettings.AdaptiveRhoInterval(), 25) #
=#
# looks like they're points with low z coordinate.

p = Point(15099.411397006435, 63329.17689456945, 489.2398433903766)
p = Point(56270.48118496456, 53948.45819773468, 34.889982642960504)
p = Point(37517.76357686068, 84362.38878880453, 233.50717418757495)
