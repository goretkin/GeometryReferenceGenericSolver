using SparseArrays: SparseMatrixCSC
using OSQP

# undo https://github.com/oxfordcontrol/OSQP.jl/blob/2e1dba7be28d7e25bac0427ac34a5f65651438bf/src/types.jl#L33

function unsafe_convert(::Type{SparseMatrixCSC}, c::OSQP.Ccsc)
  m = c.m
  n = c.n
  nzmax = c.nzmax
  nzval = [unsafe_load(c.x, i) for i=1:nzmax]
  rowval = [unsafe_load(c.i, i) for i=1:nzmax] .+ 1
  colptr = [unsafe_load(c.p, i) for i=1:nzmax] .+ 1
  SparseMatrixCSC(m, n, colptr, rowval, nzval)
end

function unsafe_load_data(d::OSQP.Data)
  A = unsafe_convert(SparseMatrixCSC, unsafe_load(d.A))
  P = unsafe_convert(SparseMatrixCSC, unsafe_load(d.P))
  l = [unsafe_load(d.l, i) for i=1:d.m]
  u = [unsafe_load(d.u, i) for i=1:d.m]
  q = [unsafe_load(d.u, i) for i=1:d.n]

  (A=A, P=P, l=l, u=u, q=q)
end

unsafe_load_data(m::OSQP.Model) = unsafe_load_data(unsafe_load(unsafe_load(r_m.workspace).data))

using LinearAlgebra: kron, I
using SparseArrays: sparse

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
    A=A, l=l, u=u, P=QP_P, q=zeros(n_variables)
  )
end

quad(x) = x' * x

function setup(c1, c2)
  data = setup_qp_data(c1, c2)
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
          #scaling = false,

          max_iter=4000,
          warm_start=false,   # for deterministic results
          polish=true,
          eps_prim_inf=1e-8,
          eps_dual_inf=1e-8,
          eps_abs=1e-8,
          eps_rel=1e-6,
          scaling=0,
          )
  model
end


function find_minimum_distance2(c1, c2)
  model = setup(c1, c2)
  OSQP.solve!(model)
  return (model=model,)
end


using GeometryReferenceGenericSolver
using GeometryTypes

xy_plane = 1e8 .* (Point(-1,0,0), Point(1,0,0), Point(0, 1, 0))

p = Point(15099.411397006435, 63329.17689456945, 489.2398433903766)
p = Point(56270.48118496456, 53948.45819773468, 34.889982642960504)
p = Point(37517.76357686068, 84362.38878880453, 233.50717418757495)

c1 = xy_plane
c2 = (p, )


using GeometryReferenceGenericSolver
#r_p = GeometryReferenceGenericSolver.find_minimum_distance(c1, c2).model.optimizer.inner
#r_m = find_minimum_distance2(c1, c2).model
