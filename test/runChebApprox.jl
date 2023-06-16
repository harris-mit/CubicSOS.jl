# Run a polynomial approximation semi-infinite example
using ApproxFun

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
fcnf4thdiv = 4.1
fcnsecdiv = .968
b = Fun(f, ApproxFun.Chebyshev()).coefficients
K = 20#length(b)
# Using both LP and SOCP
tlens = 2 .^ (2:9)
errors_socp = zeros(length(tlens))
errors_lp = zeros(length(tlens))
for i = 1:length(tlens)
    ts = range(-1, 1, length = tlens[i])
    x_socp = get_cheb_approx(ts, f, fp, fcnf4thdiv, K)
    x_lp = get_cheb_approx_lp(ts, f, fp, fcnsecdiv, K)
    errors_socp[i] =  norm(b[1:K] - x_socp)/norm(b[1:K])
    errors_lp[i] = norm(b[1:K] - x_lp)/norm(b[1:K])
end

# Approximate erf(x) with bump functions
K = 19
ts = range(-1,1,length = 32)
bumps = [t -> if (abs(t-i) < 1) exp(-1/(1 - (t-i)^2)) else 0 end for i = range(-1,1,length=K)]
bumpsderivs = [t -> if (abs(t-i) < 1) 2 * (i-t) * exp(-1/(1 - (t-i)^2))/(-1 + (t-i)^2) else 0 end for i = range(-1,1,length=K)]
x_socp = get_bump_approx(ts, f, fp, fcnf4thdiv, bumps, bumpsderivs)
x_lp = get_bump_approx_lp(ts, f, fcnsecdiv, bumps)
print(string(x_socp)[2:(end - 1)])
print(string(x_lp)[2:(end - 1)])