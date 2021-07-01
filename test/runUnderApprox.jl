# Underapproximation with CubicSOS
using CubicSOS, Test
using Mosek, MosekTools
using JuMP

m = x -> (1 - 2 * x + x^2)
n = x -> (1 + x + x^2)
r = x -> m(x) / n(x)
nderiv = x -> (1 + 2*x)
mderiv = x -> (-2 + 2 * x)

interp_pts = -2:1:2
num_pts = length(interp_pts)

model = Model()
set_optimizer(model, () -> Mosek.Optimizer(
    MSK_IPAR_INTPNT_SOLVE_FORM = MSK_SOLVE_DUAL))
@variable(model, l)
@variable(model, f_vals[1:num_pts])
@variable(model, f_derivs[1:num_pts])
@objective(model, Min, l)
# 4th derivatives of m & n are zero, otherwise use sup over interval as RHS
# 4th derivative of spline is always 0
#TODO: What about 4th derivative of n with spline??
p = CubicSpline(interp_pts, f_vals, f_derivs) # cubic spline approximation
nspline = CubicSpline(interp_pts, n.(interp_pts), nderiv.(interp_pts))
mspline = CubicSpline(interp_pts, m.(interp_pts), mderiv.(interp_pts))
# Test we constructed the spline right
xx = -2:.01:2
mapprox = evaluate_cubic.(Ref(mspline), xx)
mtrue = m.(xx)
norm(mapprox - mtrue)
constrain_spline_nonnegative!(model, nspline * l - mspline + p * nspline)
constrain_spline_nonnegative!(model, nspline * l + mspline - p * nspline)
constrain_spline_nonnegative!(model, mspline - nspline * p)
#constrain_spline_nonnegative!(model, 2-p)
optimize!(model)
pvals = value.(evaluate_cubic.(Ref(p), xx))

plot(xx, pvals)
plot(xx, r.(xx))

@test maximum(abs, value.(evaluate_cubic.(Ref(p), xx)) .- r.(xx)) / maximum(abs, r.(xx)) < 1e-2
all(value.(evaluate_cubic.(Ref(mspline - nspline * p), xx)) .> 0 )
