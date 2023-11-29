# Run a polynomial approximation semi-infinite example
using ApproxFun
using Remez # For best minimax approximation

# An approximation of sin(x) * exp(x)
ts = range(-1, 1, length = 5)
f = x -> sin(x) * exp(x)
fp = x -> cos(x) * exp(x) + sin(x) * exp(x)
fcnf4thdiv = 12
b = Fun(f, ApproxFun.Chebyshev()).coefficients
K = 7#length(b)
x = get_cheb_approx(ts, f, fp, fcnf4thdiv, K)
# Given fourth derivative of f = -4sin(x)e^x < 12 on [0,1]
@show norm(b[1:K] - x)/norm(x)# pretty close
# with K larger the 4th derivative estimate gets worse


# An approximation of erf(x)
f = x -> erf(x)
fp = x -> 2/sqrt(pi) * exp(-x^2)
fcnf4thdiv = 4.41
fcnsecdiv = .968
# Old way to compare parameters
#b = Fun(f, ApproxFun.Chebyshev()).coefficients
K = 14 #length(b)
# Get minimax optimal approximation via Remez
N,D,E,X = ratfn_minimax(f, [-1,1], K, 0);
mmopt = x -> sum([x^(i-1) * N[i] for i = 1:lastindex(N)]) / D[1]
# Using both LP and SOCP
tlens = 2 .^ (2:6) # larger t goes beyond numerical precision of Mosek defaults...
errors_socp = zeros(length(tlens))
errors_lp = zeros(length(tlens))
# beyond this error doesn't really make sense...
for i = 1:length(tlens)
    ts = range(-1, 1, length = tlens[i])
    x_socp, err_socp = get_cheb_approx(ts, f, fp, fcnf4thdiv, K)
    x_lp, err_lp = get_cheb_approx_lp(ts, f, fp, fcnsecdiv, K)
    socp_test = eval_cheb.(Ref(x_socp), testgrid, Ref(0))
    lp_test = eval_cheb.(Ref(x_lp), testgrid, Ref(0))
    errors_socp[i] = err_socp
    errors_lp[i] = err_lp
end
[tlens errors_lp errors_socp]
# Approximate erf(x) with gaussian "bump" functions (not truly bumps...)
K = 10
bumps = [t -> exp(-(t-i)^2) for i = range(-1,1,length=K)]
bumpsderivs = [t -> -2*(t-i)*exp(-(t-i)^2) for i = range(-1,1,length=K)]
tlens = 2 .^ (2:6) 
errors_socp = zeros(length(tlens))
errors_lp = zeros(length(tlens))
for i = 1:lastindex(tlens)
    ts = range(-1, 1, length = tlens[i])
    #x_socp, err_socp = get_bump_approx(ts, f, fp, fcnf4thdiv, bumps, bumpsderivs)
    x_lp, err_lp = get_bump_approx_lp(ts, f, fcnsecdiv, bumps)
    errors_socp[i] = err_socp
    errors_lp[i] = err_lp
end
[tlens errors_lp errors_socp]

# Certify the Newman bound for rational approximation to |x|
h_lp = []
h_socp = []
for n = 4:12
    println(n)
    push!(h_socp, get_num_required_point(n, certify_newman_bound))
    println("h_socp", h_socp)
    push!(h_lp, get_num_required_point(n, certify_newman_bound_lp))
    println("h_lp", h_lp)
end
certify_newman_bound(4, 1/2)
certify_newman_bound_lp(n, 1/2)

# Find best Newman bound
h = .0001
t_lp = []
t_socp = []
for n = 4:12
    println(n)
    push!(t_socp, CubicSOS.optimize_newman_bound(n, h))
    println("t_socp", t_socp)
    push!(t_lp, CubicSOS.optimize_newman_bound_lp(n, h))
    println("t_lp", t_lp)
end
# LP can't do it with this many grid points...

# Now estimate the best rational approximation of degree n
# Some test cases:
n = 10
h = .1
t = .01
C2M = get_cheb_basis_mat(n, 0)

issolved, pcoefs, qcoefs = find_best_rational_approx(n, h, t)
p = C2M * pcoefs
q = C2M * qcoefs  
# Print as strings
join(string.(p), ",")
join(string.(q), ",")

# Produce curves for paper:
# Find best rational approx guarantees with varying number of grid points
n = 8
# 20k grid points 
hs = [.1, .01, .001, .0001]
t_opt = zeros(length(hs))
for hi = 1:lastindex(hs)
    h = hs[hi]
    println(t_opt)
    t_opt[hi] = find_best_rational_approx(n, h)
end

rt_opt = zeros(length(ns)) # rational optimal t
nt_opt = zeros(length(ns)) # Newman optimal t
ns = 4:7
h = .01
for ni = 1:lastindex(ns)
    println(ni)
    @show rt_opt
    @show nt_opt
    n = ns[ni]
    rt_opt[ni] = find_best_rational_approx(n+1, h)
    # If Newman is degree n, then the actual degree of rational approximation is n+1.
    nt_opt[ni] = CubicSOS.optimize_newman_bound(n, h)
end

n = 7; h = .001;
t = find_best_rational_approx(n, h)
solved, p, q = find_rational_approx(n, h, t + 1e-3)
C2M = get_cheb_basis_mat(n, 0)
p = C2M * p
q = C2M * q  
# Print as strings
join(string.(p), ",")
join(string.(q), ",")