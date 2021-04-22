
# We precompute the value of the quadrature for each interval and for each b_j(x)

# Define grids
xs = 0:0.1:3
ys = 0:.025:15
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
populate_fourier_integrals!(n, Fcoefs, Fpcoefs, F4bnd, Fcoefserr, Fpcoefserr)
F4bnd = omega * trig_integral * 16 * pi^4 * F4bnd;
rad = 1.4# This means we want to constrain f \leq 0 whenever |x| > rad
Q1s, Q2s, X1s, X2s = check_bound_feasibility(rad, Fcoefs, Fpcoefs, F4bnd, num_xs, num_ys);


ENV["MPLBACKEND"]="qt5agg"
using PyPlot
pygui(true)

num_plot_pts = 1000;
fig, ax = get_function_plots(Q1s, Q2s, X1s, X2s, xs, ys, num_plot_pts);
