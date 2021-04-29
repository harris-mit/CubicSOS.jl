# Unit tests for the interpolatingCubic.jl methods

"""
 Tests:
It should fit known data in an optimization problem.
We should be able to enforce the PSD constraint.
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
cs = CubicSpline(x_vals, var_y_vals, var_deriv_vals)

# Test that evaluation of the cubic at endpoints gives back the desired nodes
@test norm(value.(evaluate_cubic.(Ref(cs), x_vals) .- cs.y_vals)) == 0

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


# BUG in MOSEK TOOLS?? Try to reproduce here

using MathProgBase, SparseArrays
using MathOptInterface
const MOI = MathOptInterface
const MOIU = MOI.Utilities
const MOIB = MOI.Bridges
T = Float64
optimizer = SCS.Optimizer()
model = MOIB.full_bridge_optimizer(
    MOIU.CachingOptimizer(
        MOIU.UniversalFallback(MOIU.Model{T}()),
        optimizer
    ),
    T
)
id_to_variables, conic_constr_to_constr, conic_constraints, var_to_indices, constraint_indices = Convex.load_MOI_model!(model, problem);


# We do not need to copy names because name-related operations are handled by `m.model_cache`
indexmap = MOI.copy_to(model.model.optimizer, model.model.model_cache, copy_names=false)
# MOI does not define the type of index_map, so we have to copy it into a
# concrete container. Also load the reverse map.
for k in keys(indexmap)
    model.model.model_to_optimizer_map[k] = indexmap[k]
    model.model.optimizer_to_model_map[indexmap[k]] = k
end

#MOI.optimize!(model) # if you run this, it clears out the problem info

problem_description = model.model.optimizer
cone = problem_description.cone
m = problem_description.data.m
n = problem_description.data.n
A = sparse(problem_description.data.I, problem_description.data.J, problem_description.data.V)
b = problem_description.data.b
c = problem_description.data.c
c, A, b, cones, var_to_ranges, vartypes, conic_constraints = Convex.conic_problem(problem)


Q1s = Variable(3 * length(x_vals)-3)
Q2s = Variable(length(x_vals) + 1)
cs = CubicSpline(x_vals, Q1s, Q2s)


problem = satisfy()
var_vals_and_derivs = get_vals_and_derivs(cs)
desired_vals_and_derivs = [y_vals y_derivs]
for i = 1:size(var_vals_and_derivs, 1)
    j = 2
    #for j = 1:size(var_vals_and_derivs, 2)
        problem.constraints += [var_vals_and_derivs[i,j] == desired_vals_and_derivs[i,j]]
    #end
end
#opt = () -> Mosek.Optimizer()
opt = () -> SCS.Optimizer()
solve!(problem, opt)

coefs = get_poly_coefficients(cs)
coefs = evaluate.(coefs)
coef_vals = zeros(size(coefs))
for i = 1:size(coefs,1)
    for j = 1:size(coefs,2)
        coef_vals[i,j] = coefs[i,j][1]
    end
end
