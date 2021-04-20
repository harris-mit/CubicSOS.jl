# # @test
#
# # We precompute the value of the quadrature for each interval and for each b_j(x)
# # First determine intervals
# xs = 0 : 0.1 : 2
# ys = 0 : 0.1 : 15
# num_xs = length(xs)
# num_ys = length(ys)
# # We store the integrals and their associated errors in a table that is
# # indexed by [yvalue, xinterval, which basis element]
# Fcoefseff = zeros(num_ys, num_xs - 1, 4)
# Fcoefserreff = zeros(num_ys, num_xs - 1, 4) # These are estimated errors norm(I - I') in an absolute sense
# Fpcoefseff = zeros(num_ys, num_xs - 1, 4) # D[F] = Fp
# Fpcoefserreff = zeros(num_ys, num_xs - 1, 4)
# F4bndeff = zeros(num_xs - 1, 4) # Calculate an upper bound on the fourth derivative over the interval
#
# n = 3 # number of dimensions in problem
# omega = (2 * pi)^((n-1)/2) / gamma((n-1)/2) # radial func notes
# #omega = pi^(n/2) / gamma(n/2 + 1) # Cohn sphere packing notes
# trig_func = x -> cos(x)^4 * sin(x)^(n-2)
# trig_integral = hquadrature(trig_func, 0, pi)[1]
#
# populate_fourier_integrals(n, Fcoefseff, Fpcoefseff, F4bndeff, Fcoefserreff, Fpcoefserreff)
# F4bndeff = omega * trig_integral * 16 * pi^4 * F4bndeff;

# Q1s, Q2s, X1s, X2s = check_bound_feasibility2(1.9, Fcoefseff, Fpcoefseff, F4bndeff, num_xs, num_ys);

@test 2 == 2
