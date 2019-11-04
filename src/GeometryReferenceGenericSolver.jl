module GeometryReferenceGenericSolver

using LinearAlgebra: norm
using LinearAlgebra: kron, I
using SparseArrays: sparse

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

function find_minimum_distance_parametron(convex1::NTuple{M1, PointT}, convex2::NTuple{M2, PointT}) where {M1, M2, PointT}
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

"""
produce U from A such that `x' * U * x == x' * A * x` ∀ x
overwrite A in-place
"""
function upper_triangularize_quadratic!(A)
  s = size(A)
  n = first(s)
  @assert s == (n, n)

  # iterate below diagonal
  for i=2:n, j=1:(i-1)
    A[j,i] += A[i,j]
    A[i,j] = 0
  end
  return nothing
end

"""
Make QP matrices for the problem of finding the minimum distance between two convex hulls.
"""
function setup_qp_data(c1, c2)
  # order all a1 and then all a2
  m1 = length(c1)
  m2 = length(c2)
  n = length(c1[1]) # assume all same dimension

  M = [m1, m2]
  M_offset = [0, m1]

  _C = [c1, c2]

  C = zeros(2*n, m1 + m2)
  for point_set_i = 1:2
    for _j = 1:M[point_set_i]
      j = _j + M_offset[point_set_i]
      n_offset = (point_set_i-1) * n
      C[(1:n) .+ n_offset, j] .= _C[point_set_i][_j]
    end
  end

  p1_minus_p2_op = kron([1, -1]', I(n))

  # maps convex combination coefficients to the difference
  p1_minus_p2__Ta = p1_minus_p2_op * C

  QP_P = 2 * (p1_minus_p2__Ta' * p1_minus_p2__Ta)

  n_constraints = (m1 + m2 + length(M))
  n_variables = m1 + m2

  A = zeros(n_constraints, n_variables)
  l = zeros(n_constraints)
  u = zeros(n_constraints)

  A[1:n_variables, 1:n_variables] = I(n_variables)
  l[1:n_variables] .= 0
  u[1:n_variables] .= Inf

  A[n_variables+1, (1:m1) ] .= ones(m1)
  A[n_variables+2, (1:m2) .+ m1 ] .= ones(m2)
  l[(n_variables+1):end] .= 1
  u[(n_variables+1):end] .= 1
  #upper_triangularize_quadratic!(QP_P)
  return (
    A=A, l=l, u=u, P=QP_P, q=zeros(n_variables),
    parts=(
        m1=m1,
        m2=m2,
        p1_minus_p2__Ta=p1_minus_p2__Ta,
        C = C,
    )
  )
end

quad(x) = x' * x

function setup_model(data)
  model = OSQP.Optimizer().inner

  #model = OSQP.Model()
  #MOI.set(model.optimizer, OSQPSettings.Scaling(), false)
  #OSQP.update_settings!(model; scaling=0)

  OSQP.setup!(model;
          P=sparse(data.P),
          q=data.q,
          A=sparse(data.A),
          l=data.l,
          u=data.u,

          verbose=false,
          max_iter=4000,
          warm_start=false,   # for deterministic results
          polish=true,
          eps_prim_inf=1e-8,
          eps_dual_inf=1e-8,
          eps_abs=1e-8,
          eps_rel=1e-6,
          scaling=10,
          )
  model
end


function find_minimum_distance(c1, c2)
  qp_problem_data = setup_qp_data(c1, c2)

  model = setup_model(qp_problem_data)
  result = OSQP.solve!(model)
  coeffs = result.x
  a1 = coeffs[1:qp_problem_data.parts.m1 .+ 0]
  a2 = coeffs[(1:qp_problem_data.parts.m2) .+ qp_problem_data.parts.m1]
  p1 = sum(c1 .* a1)
  p2 = sum(c2 .* a2)

  return (
    d_sq=result.info.obj_val, p1=p1, p2=p2,
    model=model, result=result)
end


end # module
