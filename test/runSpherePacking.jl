# This file is a driver for SpherePacking.jl, the main second order implementations
# of the sphere packing problems. It is not use for any testing but rather is
# a demonstration of how the various functions can be called.

# First solve sphere packing in a lattice, as solved by Cohn et al
# Second we solve sphere packing on a compact set, as done for Spherical codes by Delsarte

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
rad = 1. # This means we want to constrain f \leq 0 whenever |x| > rad
f, fhat, delta = solve_lattice_sphere_packing_socp(rad, Fcoefs, Fpcoefs, F4bnd, xs, ys)
fmax, fhatmin = check_feasibility(rad, f, fhat, xs, ys)
@show fhatmin - value(delta)
# Look at where the problem changes from feasible to infeasible at 1.0 to 1.1.
@show get_opt_radius(n)


ENV["MPLBACKEND"]="qt5agg"
using PyPlot
pygui(true)
num_plot_pts = 100
fig, ax = get_function_plots(f, fhat, xs, ys, rad, num_plot_pts)


# Solve the Delsarte sphere packing in a sphere.
# This version of the method has some numerical issues because we need to
# compute the coefficients which are very large.
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



# Now we compute Delsarte's lower bounds with a more stable method
# This only involves polynomial evaluation, so we don't deal with large coefficients.
# The method sets the Gegenbauer coefficients to be the decision variables.
# We use a low degree cubic spline to certify the nonpositivity condition.
Amax = cos(pi / 3)
xs = range(-1, Amax, length =75)
d = 9
kmax = 50
fk = solve_delsarte_socp2(d, kmax, xs)
# Testing our certificate:
gktrunc = value.(fk)[1:14] # Truncating off the ones that are basically 0 (and possibly slightly negative.)
gktrunc[1] = 1
all(gktrunc .> 0 )
compute_from_expansion(d, gktrunc, 1)
xx = range(-1, Amax, length =1000)
yy = compute_from_expansion.(Ref(d), Ref(gktrunc), xx)
maximum(yy)
