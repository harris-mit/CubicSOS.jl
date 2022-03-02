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

