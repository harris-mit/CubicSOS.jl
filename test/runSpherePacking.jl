# First solve sphere packing in a lattice, as solved by Cohn et al
# We precompute the value of the quadrature for each interval and for each b_j(x)

# Define grids
xs = 0:0.05:2 # Larger intervals is larger numbers in matrices, but approx accuracy shrinks
ys = 0:0.1:20
# Needs to be fine enough to allow interpolation
# & long enough to approximate infinite length
# If too fine sometimes it hurts solver's ability to determine status as well
num_xs = length(xs)
num_ys = length(ys)
# We store the integrals and their associated errors in a table that is
# indexed by [yvalue, xinterval, which basis element]
Fcoefs = zeros(num_ys, num_xs - 1, 4)
Fcoefserr = zeros(num_ys, num_xs - 1, 4) # These are estimated errors norm(I - I') in an absolute sense
Fpcoefs = zeros(num_ys, num_xs - 1, 4) # D[F] = Fp
Fpcoefserr = zeros(num_ys, num_xs - 1, 4)
F4bnd = zeros(num_xs - 1, 4) # Calculate an upper bound on the fourth derivative over the interval

n = 2 # number of dimensions in problem
omega = 2 * pi^((n - 1) / 2) / gamma((n - 1) / 2) # radial func notes surface area of sphere
# omega = pi^(n/2) / gamma(n/2 + 1) # Cohn sphere packing notes for VOLUME -- not what we want
# we need the following integral for the fourth derivative bound
trig_func = x -> cos(x)^4 * sin(x)^(n - 2)
trig_integral = hquadrature(trig_func, 0, pi)[1]
populate_fourier_integrals!(n, xs, ys, omega, Fcoefs, Fpcoefs, F4bnd, Fcoefserr, Fpcoefserr)
F4bnd = omega * trig_integral * 16 * pi^4 * F4bnd;

# Improve conditioning of matrix by zeroing out non-contributors
# thresh_tol = 1e-12
# Fcoefs[abs.(Fcoefs) .< thresh_tol] .= 0
# Fpcoefs[abs.(Fpcoefs) .< thresh_tol] .= 0
# F4bnd[abs.(F4bnd) .< thresh_tol] .= 0
rad = 1. # This means we want to constrain f \leq 0 whenever |x| > rad
f, fhat, delta = solve_sphere_packing_socp(rad, Fcoefs, Fpcoefs, F4bnd, xs, ys)
fmax, fhatmin = check_feasibility(rad, f, fhat, xs, ys)
@show fhatmin - value(delta)
# The problem changes from feasible to infeasible at 1.0 to 1.1.
# At higher resolution, changes from 1.05 to 1.1...
# Needed to increase the domain of \hat{y}
# The true optimal density corresponds to a radius of 1.07.
# Shrink domain of x to improve values/conditioning...
@show get_opt_radius(n)


ENV["MPLBACKEND"]="qt5agg"
using PyPlot
pygui(true)
num_plot_pts = 100
fig, ax = get_function_plots(f, fhat, xs, ys, rad, num_plot_pts)


# Solve the Delsarte sphere packing in a sphere
Amax = .5
d = 3
deltax = .1
xs = -1:deltax:1
num_xs = length(xs)
num_refined_gegenbauer = 100
Gcoefs = zeros(num_refined_gegenbauer, num_xs - 1, 4)
populate_gegenbauer_transform_d3!(Gcoefs, d, xs)

num_gegenbauer = 28
f = solve_delsarte_socp(Amax, Gcoefs[1:num_gegenbauer,:,:], xs)
value(f.y_vals[end])
C = get_hermite_basis_coefficients(f)
gk = zeros(num_refined_gegenbauer)
for k = 1:size(Gcoefs, 1)
    gk[k] = value(sum(Gcoefs[k, :, :] .* C))
    @show k, gk[k]
end

# integrate(b * scaled_gegenbauer(3, 30), .9, 1) blows up at k >= 30
xx = -1:.01:1
yy = value.(evaluate_cubic.(Ref(f), xx))
plot(xx,yy)
scatter(1:30, gk[1:30], legend = false)
plot(xx, scaled_gegenbauer(3,50).(xx))
# When k > unrefined it goes from -1e-9 to -1e-3 :(
# They proved 13, we have 18 here, an upper bound. Looks good
# Decreasing delta too small made infeasible? Integrals too small?
# bound increases when we increase k. Need a better sense
# for how many coefficients. But right order of magnitude.

# 21.8 for k < 30
# 22.5 for k < 35

# Check truncated feasibility
xx = -1:.01:Amax
gktrunc = gk[1:num_gegenbauer]
gktrunc[gktrunc .< 0] .= 0
gkyy = compute_from_expansion.(Ref(d), Ref(gktrunc), xx)
fyy = value.(evaluate_cubic.(Ref(f), xx))
sum(gkyy .> 0*1e-8)
sum(fyy .> 0*1e-8)
sum(gktrunc .< -1e-8)
sum(gktrunc .< 0)
compute_from_expansion(d, gktrunc, 1)
value(evaluate_cubic(f, 1))
gk[1]

Amax = cos(pi / 3)
xs = -1:.005:Amax
# Warning: we need an "exact hit" where interval length matches Amax
d = 9
kmax = 50
fk = solve_delsarte_socp2(d, kmax, xs)
gktrunc = value.(fk)[1:7]
#gktrunc[gktrunc .< 0] .= 0 # !!! YOU CAN'T JUST SET SOME OF THEM TO ZERO. BETTER TO TRUNCATE!!
gktrunc[1] = 1
compute_from_expansion(d, gktrunc, 1)
xx = -1:.01:Amax
yy = compute_from_expansion.(Ref(d), Ref(gktrunc), xx)
maximum(yy)

yexpanded = compute_from_expansion.(Ref(d), Ref(value.(fk[1:7])), xx)
ycubic = value.(evaluate_cubic.(Ref(f), xx))
maximum(abs, yexpanded - ycubic) # This is less than delta = .007!
plot(xx, yexpanded)
scatter!(xs, compute_from_expansion.(Ref(d), Ref(value.(fk[1:7])), xs))
plot!(xx, ycubic)
