module GeometryReferenceGenericSolver

using LinearAlgebra: norm

using Parametron
using Parametron:setdirty!
using OSQP
using OSQP.MathOptInterfaceOSQP: OSQPSettings

import MathOptInterface
const MOI = MathOptInterface

function _constrain_convex_combination_coefficients(model, coeffs)
    M = length(coeffs)
    for i = 1:M
        @constraint(model, coeffs[i] >= 0)
    end
    s = sum(coeffs)
    @constraint(model, s == 1)
end

function _convex_combination(coeffs, points, m)
    s = @expression coeffs[1] * points[:, 1]
    for i = 2:m
        s = @expression s + coeffs[i] * points[:, i]
    end
    return s
end

function squared_distance(p1, p2, n)
    s = @expression (p1[1] - p2[1])^2
    for i = 2:n
        s = @expression s + (p1[i] - p2[i])^2
    end
    return s
end

struct SolverError{T} <: Exception
    status::T
end

function find_minimum_distance(convex1::NTuple{M1, PointT}, convex2::NTuple{M2, PointT}) where {M1, M2, PointT}
    model = Model(OSQP.Optimizer(
        max_iter=4000,
        warm_start=false,   # for deterministic results
        polish=true,
        eps_prim_inf=1e-8,
        eps_dual_inf=1e-8,
        eps_abs=1e-8,
        eps_rel=1e-6,
        scaling=0)
    )

    # from https://github.com/tkoolen/Parametron.jl/blob/d824c9f9ff4ba7bc38a9407129504bb8ecca31ee/test/model.jl#L22
    #MOI.set(model.optimizer, OSQPSettings.AdaptiveRhoInterval(), 25) # required for deterministic behavior
    #MOI.set(model.optimizer, OSQPSettings.Verbose(), false)
    #MOI.set(model.optimizer, OSQPSettings.Scaling(), false)
    N = length(PointT)      # TODO works for GeometryType.Point, general trait?
    # points
    c1 = Parameter(model, val=zeros(N, M1))
    c2 = Parameter(model, val=zeros(N, M2))

    # TODO use broadcast machinery to set ]up problem
    # coefficients of convex combination
    a1 = [Variable(model) for _ in 1:M1]
    a2 = [Variable(model) for _ in 1:M2]

    _constrain_convex_combination_coefficients(model, a1)
    _constrain_convex_combination_coefficients(model, a2)

    p1 = _convex_combination(a1, c1, M1)
    p2 = _convex_combination(a2, c2, M2)
    d_sq = squared_distance(p1, p2, N)

    @objective(model, Minimize, d_sq)

    # set parameters
    for i = 1:M1
        c1.val[][:,i] = convex1[i]
    end
    for i = 1:M2
        c2.val[][:,i] = convex2[i]
    end

    setdirty!(c1)
    setdirty!(c2)
    #=data_before_solve = deepcopy(
     unsafe_load_data(
        unsafe_load(
            unsafe_load(model.optimizer.inner.workspace).data))
    )
    =#
    data_before_solve = "syke"
    solve!(model)
    status =  MOI.get(model.optimizer, MOI.TerminationStatus())

    if status != MOI.OPTIMAL
        #throw(SolverError(status))
    end

    # TODO remove after done debugging
    g = value(model, a1[1])

    # extract decision varibles and re-compute expressions
    # TODO avoid expressing convex combinations and norm twice
    v_a1 = value.(model, a1)
    v_a2 = value.(model, a2)

    v_p1 = sum(v_a1 .* convex1)
    v_p2 = sum(v_a2 .* convex2)
    v_d_sq = sum((v_p1 - v_p2).^2)
    return (d_sq = v_d_sq, p1=v_p1, p2=v_p2, model=model,
        model_parts = (
            a1=a1, a2=a2, c1=c1, c2=c2,
        ),
        data_before_solve=data_before_solve
    )
end

end # module
