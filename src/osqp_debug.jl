using SparseArrays: SparseMatrixCSC
using OSQP

# undo https://github.com/oxfordcontrol/OSQP.jl/blob/2e1dba7be28d7e25bac0427ac34a5f65651438bf/src/types.jl#L33

function unsafe_convert(::Type{SparseMatrixCSC}, c::OSQP.Ccsc)
  m = c.m
  n = c.n
  nzmax = c.nzmax
  nzval = [unsafe_load(c.x, i) for i=1:nzmax]
  rowval = [unsafe_load(c.i, i) for i=1:nzmax] .+ 1
  colptr = [unsafe_load(c.p, i) for i=1:(n+1)] .+ 1
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




using GeometryReferenceGenericSolver
using GeometryTypes

xy_plane = 1e6 .* (Point(-1,0,0), Point(1,0,0), Point(0, 1, 0))

p = Point(15099.411397006435, 63329.17689456945, 489.2398433903766)
p = Point(56270.48118496456, 53948.45819773468, 34.889982642960504)
p = Point(37517.76357686068, 84362.38878880453, 233.50717418757495)

c1 = xy_plane
c2 = (p, )


using GeometryReferenceGenericSolver
#r_p = GeometryReferenceGenericSolver.find_minimum_distance(c1, c2).model.optimizer.inner
#r_m = find_minimum_distance2(c1, c2).model
