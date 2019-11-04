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

using LinearAlgebra: kron I
using SparseArrays: sparse

function setup(c1, c2)
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

  n_constraints = (m1 + 1 + m2 + 1)
  n_variables = m1 + m2

  A = zeros(n_constraints, n_variables)
  l = zeros(n_constraints)
  u = zeros(n_constraints)

  for point_set_i = 1:2
    m = M[point_set_i]
    _M_offset = M_offset[point_set_i]
    n_offset = (point_set_i-1) * n
    @show point_set_i
    @show (1:m) .+ _M_offset
    @show (m + 1) + _M_offset
    _mrange = (1:m) .+ (_M_offset + (point_set_i-1))
    A[_mrange, (1:m) .+ _M_offset] = I(m)
    A[last(_mrange) + 1, (1:m) .+ n_offset] = ones(1, m)
    l[_mrange] .= 0
    u[_mrange] .= Inf
    l[last(_mrange) + 1] = 1
    u[last(_mrange) + 1] = 1
  end

  model = OSQP.Model()
  #MOI.set(model.optimizer, OSQPSettings.Scaling(), false)
  #OSQP.update_settings!(model; scaling=0)

  OSQP.setup!(model;
          P=sparse(QP_P),
          q=zeros(n_variables),
          A=sparse(A),
          l=l,
          u=u,
          #scaling = false
          )
  model
end


function find_minimum_distance2(c1, c2)
  model = setup(c1, c2)
  OSQP.solve!(model)
  return (model=model)
end
