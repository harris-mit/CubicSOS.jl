# Simple demonstration of package
using CubicSOS
using JuMP, MosekTools
N = 5; mu4p = 10 # bound on 4th derivative of p
mu4f = N * (N-1) * (N-2)
knots = range(-1, 1, length = 20)
xN = knots .^ N # true values of x^N
xNp = N * knots .^ (N-1) # true values of x^N's derivative

model = Model(); set_optimizer(model, () -> Mosek.Optimizer())
@variable(model, a[1:N]) # coefficients of approximation
approx = [sum([a[i+1] * knots[j] ^ i for i = 0:N-1]) for j = 1:lastindex(knots)] # value of approximation
approxderiv = [sum([a[i+1] * i * knots[j] ^ (i-1) for i = 1:N-1]) for j = 1:lastindex(knots)] # value of approximation's derivative
p = CubicInterpolant(knots, approx, approxderiv)
f = CubicInterpolant(knots, xN, xNp)

@variable(model, error)
@variable(model, delta)
@constraint(model, delta >= maximum(diff(knots))^4 * (mu4p + mu4f) / 384)
constrain_interpolant_nonnegative!(model, p - f + error - delta) # certifies p - f <= t
constrain_interpolant_nonnegative!(model, f - p + error - delta) # certifies f - p <= t
@objective(model, Min, error); optimize!(model)