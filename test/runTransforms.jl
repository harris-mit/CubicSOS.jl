# Test the Gegenbauer expansion of a known function
using Test, JuMP, LinearAlgebra
d = 4
deltax = .1
xs = -1:deltax:1
num_xs = length(xs)
f = x -> x^2 + sin(x)
fp = x -> 2 * x + cos(x)
p = CubicInterpolant(xs, f.(xs), fp.(xs))
num_gegenbauer = 20
Gcoefs = zeros(num_gegenbauer, num_xs - 1, 4)
populate_gegenbauer_transform!(Gcoefs, d, xs, 1e-2)
C = get_hermite_basis_coefficients(p)
gk = zeros(num_gegenbauer)
for k = 1:num_gegenbauer
    gk[k] = sum(Gcoefs[k,:,:] .* value.(C))
end

xx = -1:.01:1
@test norm(f.(xx) .- compute_from_expansion.(Ref(d), Ref(gk), xx)) < 1e-4


# Test the symbolic computation matches the numerical one:
xs = [-1; .5]
num_xs = 2
num_refined_gegenbauer  = 3
Gcoefs = zeros(num_refined_gegenbauer, num_xs - 1, 4)
Gcoefs_analytic = zeros(num_refined_gegenbauer, num_xs - 1, 4)

populate_gegenbauer_transform!(Gcoefs, d, xs, 1e-5)
populate_gegenbauer_transform_analytic!(Gcoefs_analytic, d, xs)
@test norm(Gcoefs - Gcoefs_analytic) < 1e-4

d = 3
populate_gegenbauer_transform!(Gcoefs, d, xs)
populate_gegenbauer_transform_d3!(Gcoefs_analytic, d, xs)
@test norm(Gcoefs - Gcoefs_analytic) < 1e-8


# Test the expansion of a polynomail

p = x -> 4 * x^3 - x^2 - 3x + 2
xx = -1:.001:1
@test norm(p.(xx) - evaluate_polynomial.(Ref([2, -3, -1, 4]), xx, 0)) < 1e-8

pderiv = x -> 12 * x^2 - 2 * x - 3
@test norm(pderiv.(xx) .- evaluate_polynomial.(Ref([2, -3, -1, 4]), xx, 1)) < 1e-8
