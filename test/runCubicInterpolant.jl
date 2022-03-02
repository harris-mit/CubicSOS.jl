# Unit tests for the interpolatingCubic.jl methods
using CubicSOS, Test
using Convex, Mosek, MosekTools
using JuMP
using SumOfSquares, DynamicPolynomials, MathOptInterface
using HCubature, SpecialFunctions

"""
 Tests:
1. It should fit known data in an optimization problem.
2. Check some basic arithmetic on the spline objects.
3. We should be able to enforce nonnegativity
"""

# Test that the 2x2 SOCP is equivalent to being SDP.
x_vals = 1:.1:4
my_func = x -> 2 * sin(x) + cos(x/10)
my_func_deriv = x -> 2 * cos(x) - sin(x/10)/10
true_vals = my_func.(x_vals)
true_derivs = my_func_deriv.(x_vals)

# Using Jump
model = Model()
set_optimizer(model, () -> Mosek.Optimizer())
@variable(model, var_y_vals[1:length(x_vals)])
@variable(model, var_deriv_vals[1:length(x_vals)])
cs = CubicInterpolant(x_vals, var_y_vals, var_deriv_vals)

x_vals_coarse = 1:.2:4

# Test that evaluation of the cubic at endpoints gives back the desired nodes
cs2 = CubicInterpolant(x_vals, true_vals, true_derivs)
@test norm((evaluate_cubic.(Ref(cs2), x_vals) .- cs2.y_vals)) == 0

@constraint(model, cs.y_vals .== true_vals)
@constraint(model, cs.deriv_vals .== true_derivs)
optimize!(model)

# Test that the optimization constraints were met:
@test norm(value.(cs.y_vals) - true_vals) == 0
@test norm(value.(cs.deriv_vals) - true_derivs) == 0

xx = minimum(x_vals):(maximum(x_vals)/1000):maximum(x_vals)
true_y = my_func.(xx)
spline_y = value.(evaluate_cubic.(Ref(cs), xx))
# calculate the true upper bound
f4bnd = 2 # Calculated analytically
delta = (x_vals[2] - x_vals[1]) ^ 4 / 384 * f4bnd
# Test that the spline is as correct as we can expect.
# (In practice, it looks like it just about meets this bound.)
@test maximum(abs, true_y - spline_y) < delta
#evaluate deriatives:
spline_y_deriv = value.(evaluate_cubic_derivative.(Ref(cs), xx))
true_derivs = my_func_deriv.(xx)
@test maximum(abs, true_derivs - spline_y_deriv) < 1e-4

# finer mesh constraint, test adding and resampling...
cs_fine = CubicInterpolant(xx, spline_y, spline_y_deriv)
@test maximum(abs, value.((cs_fine + cs).y_vals) - 2 * true_y) < 1e-4
@test maximum(abs, value.((cs_fine + cs).deriv_vals) - 2 * true_derivs) < 1e-4
@test maximum(abs, value.((cs_fine - cs).y_vals)) < 1e-15
@test maximum(abs, value.((cs_fine - cs).deriv_vals)) < 1e-14

# test 2x2 psd constraint
model = Model()
set_optimizer(model, () -> Mosek.Optimizer())
@variable(model, x)
@variable(model, y)
c = 5 # any random number here will do
mat = [x, c, y]
constrain_2x2_psd!(model, mat)
@objective(model, Min, x+y)
optimize!(model)
# We know the analytical solution to the problem xy >= 1, min x+y
@test value(x) == c
@test value(y) == c

# Now test the nonnegativity constraint:
model = Model()
set_optimizer(model, () -> Mosek.Optimizer())
@variable(model, var_y_vals[1:length(x_vals)])
@variable(model, var_deriv_vals[1:length(x_vals)])
cs = CubicInterpolant(x_vals, var_y_vals, var_deriv_vals)
constrain_interpolant_nonnegative!(model, cs)

true_derivs = my_func_deriv.(x_vals)
@constraint(model, cs.deriv_vals .== true_derivs)
@variable(model, err)
@constraint(model, [err; cs.y_vals .- true_vals] in SecondOrderCone())
@objective(model, Min, err)
set_optimizer(model, () -> Mosek.Optimizer(
    MSK_DPAR_INTPNT_CO_TOL_DFEAS = 1e-10,
    MSK_DPAR_INTPNT_CO_TOL_PFEAS = 1e-16,
    ))
optimize!(model)

# This should be >= 0 everywhere, but should be reasonably close in values...
xx = minimum(x_vals):(maximum(x_vals)/1000):maximum(x_vals)
spline_y = value.(cs.(xx))
@test all(spline_y .> 0) # This is the solver toleance.

# We can plot the functions with the below codes so we know they only become
# negative after 3.5. Is it reasonably close before that?

# plot(xx, spline_y)
# plot(xx, my_func.(xx))

xx = minimum(x_vals):(maximum(x_vals)/1000):3.5
spline_y = value.(evaluate_cubic.(Ref(cs), xx))
@test maximum(abs,spline_y - my_func.(xx)) < 1e-6
# There is no theoretical answer for this because we're just minimizing the
# total interpolation error, so it could force larger errors elsewhere in the domain.
