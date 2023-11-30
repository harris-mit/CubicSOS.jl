using MathOptInterface

# Implements a polynomial approximation to a given function
# in the Chebyshev or sup norm. We also use a Chebyshev basis
# for fun!

# We test against the commented out method. Match within 10^-14
# # From TRANSFORMATION OF CHEBYSHEV–BERNSTEIN POLYNOMIAL BASIS
# # by ABEDALLAH RABABAH, Thm 1
# function cheb2bernstein(K)
#     C2B = zeros(K+1, K+1)
#     for j in 0:K
#         for k in 0:K
#             C2B[j+1,k+1] = 1/binomial(K,j) * sum([(-1)^(k-i) * binomial(2*k,2*i) * 
#             binomial(K-k,j-i) for i = max(0,j+k-K):min(j,k)])
#         end
#     end
#     return C2B
# end

"""
    cheb2bernstein(K, deriv_order)
Returns the matrix that takes Chebyshev coefficients to
coefficients in the Bernstein basis
"""
function cheb2bernstein(K, deriv_order)
    C2M = get_cheb_basis_mat(K, deriv_order);
    B2M = get_bernstein_basis_mat(K);
    return inv(B2M) * C2M
end

"""
    get_bernstein_basis_mat(K)
Gets a matrix for Bernstein to monomials
"""
function get_bernstein_basis_mat(K)
    B2M = zeros(K, K)
    xvar = Polynomials.Polynomial([0,1])
    for i in 0:K-1
        thisp = binomial(K-1, i)/2^(K-1) * (1+xvar)^i * (1-xvar)^(K-i-1)
        B2M[:,i+1] = thisp.coeffs
    end
    return B2M
end
"""
    get_cheb_basis_mat(K)
Gets a matrix for Cheb to monomials
"""
function get_cheb_basis_mat(K, deriv_order)
    C2M = zeros(K, K)
    for i in 1:K
        ind = zeros(i)
        ind[i] = 1
        thisp = derivative(convert(Polynomials.Polynomial, Polynomials.ChebyshevT(ind)), deriv_order)
        C2M[1:max((i - deriv_order),1),i] = thisp.coeffs
    end
    return C2M
end

"""
    get_cheb_basis(K, deriv_order)
Get a list of the K first Chebyshev polynomials.
deriv_order is the derivative order of the polynomials
"""
function get_cheb_basis(K, deriv_order)
    polys = []
    for i in 1:K
        ind = zeros(i)
        ind[i] = 1
        push!(polys, derivative(convert(Polynomials.Polynomial, Polynomials.ChebyshevT(ind)), deriv_order))
    end
    return polys
end

"""
    eval_cheb(x, pt, deriv_order)
Take a variable of coefficients x and evaluate the
Chebyshev polynomial at point pt.
deriv_order is the derivative order of the Chebyshev polynomials to use
"""
function eval_cheb(x, pt, deriv_order)
    polys = get_cheb_basis(length(x), deriv_order)
    return sum(x .* [p(pt) for p in polys])
end

"""
    get_cheb_fourth_div_bound(K)
Returns a naive fourth derivative bound on the
Chebyshev coefficients.
"""
function get_cheb_fourth_div_bound(K)
    polys = get_cheb_basis(K, 4)
    return [sum(abs, p.coeffs) for p in polys]
end

"""
    get_cheb_second_div_bound(K)
Returns a naive second derivative bound on the
Chebyshev coefficients.
"""
function get_cheb_second_div_bound(K)
    polys = get_cheb_basis(K, 2)
    return [sum(abs, p.coeffs) for p in polys]
end

"""
    get_cheb_approx(ts, f, fp, fcnf4thdiv)
To approximate a function f with the best estimate in the
sup norm over [ts[1], ts[end]]. Derivative of function is fp. fcnfthdiv is a bound
on the fourth derivative of f. ts are the points used for the spline.
K is the number of polynomials (in this case Chebyshev) to use.
"""
function get_cheb_approx(ts, f, fp, fcnf4thdiv, K)
    model = Model()
    C42B = cheb2bernstein(K, 4)
    set_optimizer(model, () -> Mosek.Optimizer(
        MSK_IPAR_INTPNT_SOLVE_FORM = MSK_SOLVE_DUAL))
    set_silent(model)
    @variable(model, x[1:K])
    @variable(model, err)
    @objective(model, Min, err)
    pspline = CubicInterpolant(ts, eval_cheb.(Ref(x), ts, Ref(0)),
                        eval_cheb.(Ref(x), ts, Ref(1)))
    fspline = CubicInterpolant(ts, f.(ts), fp.(ts))
    @variable(model, delta)
    @variable(model, chebfrthdivbound)
    deltax = maximum(diff(ts))
    p4 = C42B * x
    @constraint(model, chebfrthdivbound .>= p4)
    @constraint(model, chebfrthdivbound .>= -p4)
    # The below bound is more naive. Some slight improvement switching the Bernstein bound.
    #@constraint(model, [chebfrthdivbound; x .* get_cheb_fourth_div_bound(K)]
    #                    in MathOptInterface.NormOneCone(length(x)+1))
    @constraint(model, delta >= deltax^4 / 384 * (fcnf4thdiv + chebfrthdivbound))
    constrain_interpolant_nonnegative!(model, fspline - pspline + err - delta)
    constrain_interpolant_nonnegative!(model, pspline - fspline + err - delta)
    optimize!(model)
    return (value.(x), value(err))
end


"""
    get_cheb_approx_lp(ts, f, fp, fcnf4thdiv)
To approximate a function f with the best estimate in the
sup norm over [ts[1], ts[end]]. Derivative of function is fp. fcnsecdiv is a bound
on the second derivative of f. ts are the points used for linear program.
K is the number of polynomials (in this case Chebyshev) to use.
"""
function get_cheb_approx_lp(ts, f, fp, fcnsecdiv, K)
    model = Model()
    C22B = cheb2bernstein(K, 2)
    set_optimizer(model, () -> Mosek.Optimizer(
        MSK_IPAR_INTPNT_SOLVE_FORM = MSK_SOLVE_DUAL))
    @variable(model, x[1:K])
    @variable(model, err)
    @objective(model, Min, err)
    @variable(model, chebsecdivbound)
    @variable(model, delta)
    deltax = maximum(diff(ts))
    p2 = C22B * x
    @constraint(model,  chebsecdivbound .>= p2)
    @constraint(model,  chebsecdivbound .>= -p2)
    # The below bound is more naive -- just the l1 norm of the coefficients.
    # Looks like a significant improvement here (at least an order of magnitude?)
    #@constraint(model, [chebsecdivbound; x .* get_cheb_second_div_bound(K)]
    #                    in MathOptInterface.NormOneCone(length(x)+1))
    @constraint(model, delta >= deltax^2 / 8 * (fcnsecdiv + chebsecdivbound))
    fvals = f.(ts)
    pvals = eval_cheb.(Ref(x), ts, Ref(0))
    @constraint(model, fvals - pvals .+ err .>= delta)
    @constraint(model, pvals - fvals .+ err .>= delta)
    optimize!(model)
    return (value.(x), value(err))
end

# Now we repeat the above but for bump functions
"""
    get_bump_approx(ts, f, fp, fcnf4thdiv, bumps, bumpsderivs)
To approximate a function f with the best estimate in the
sup norm over [ts[1], ts[end]]. Derivative of function is fp. fcnfthdiv is a bound
on the fourth derivative of f. ts are the points used for the spline.
Bumps are a list of the bump functions we should use (doesn't actually have to be a bump...
i.e. no need for compact support)
"""
function get_bump_approx(ts, f, fp, fcnf4thdiv, bumps, bumpsderivs)
    model = Model()
    set_optimizer(model, () -> Mosek.Optimizer(
        MSK_IPAR_INTPNT_SOLVE_FORM = MSK_SOLVE_DUAL))
    @variable(model, x[1:length(bumps)])
    @variable(model, err)
    @objective(model, Min, err)
    pspline = CubicInterpolant(ts, hcat([b.(ts) for b in bumps]...)*x,
                                   hcat([b.(ts) for b in bumpsderivs]...)*x)
    fspline = CubicInterpolant(ts, f.(ts), fp.(ts))
    @variable(model, delta)
    @variable(model, chebfrthdivbound)
    deltax = maximum(diff(ts))
    # The below bound is more naive than Bernstein used above.
    # Each bumps' 4th derivative is bounded by 12
    @constraint(model, [chebfrthdivbound; x .* 12]
                        in MathOptInterface.NormOneCone(length(x)+1))
    @constraint(model, delta >= deltax^4 / 384 * (fcnf4thdiv + chebfrthdivbound))
    constrain_interpolant_nonnegative!(model, fspline - pspline + err - delta)
    constrain_interpolant_nonnegative!(model, pspline - fspline + err - delta)
    optimize!(model)
    return value.(x), value(err)
end


"""
    get_bump_approx_lp(ts, f, fcnsecdiv, bumps)
To approximate a function f with the best estimate in the
sup norm over [ts[1], ts[end]]. Derivative of function is fp. fcnsecdiv is a bound
on the second derivative of f. ts are the points used for linear program.
"""
function get_bump_approx_lp(ts, f, fcnsecdiv, bumps)
    model = Model()
    #set_optimizer(model, () -> Mosek.Optimizer(
    #    MSK_IPAR_INTPNT_SOLVE_FORM = MSK_SOLVE_DUAL))
    set_optimizer(model, () -> SCS.Optimizer())
    @variable(model, x[1:length(bumps)])
    @variable(model, err)
    @objective(model, Min, err)
    @variable(model, chebsecdivbound)
    @variable(model, delta)
    deltax = maximum(diff(ts))
    #The below bound is more naive -- just the l1 norm of the coefficients.
    # Each bumps second derivative is bounded by 2, max occurs at origin
    @constraint(model, [chebsecdivbound; x .* 2]
                        in MathOptInterface.NormOneCone(length(x)+1))
    @constraint(model, delta >= deltax^2 / 8 * (fcnsecdiv + chebsecdivbound))
    fvals = f.(ts)
    pvals = hcat([b.(ts) for b in bumps]...)*x
    @constraint(model, fvals - pvals .+ err .>= delta)
    @constraint(model, pvals - fvals .+ err .>= delta)
    optimize!(model)
    return value.(x), value(err)
end

"""
    certify_newman_bound(n)
Certify the order N newman bound for rational approximation to |x|
"""
function certify_newman_bound(n, h)
    ts = -1:h:1
    @assert 0 in ts "We need 0 in ts for |x| to be differentiable on each subinterval"
    model = Model()
    set_optimizer(model, () -> Mosek.Optimizer(
        MSK_IPAR_INTPNT_SOLVE_FORM = MSK_SOLVE_DUAL))
    @variable(model, pspline1err)
    h = maximum(diff(ts))
    z = exp(-n^(-1/2))
    xvar = Polynomials.Polynomial([0,1])
    p = prod([xvar + z^k for k = 0:(n-1)])
    pp = derivative(p, 1)
    # Constraint that p(x) + p(-x) >= 0
    pspline1 = CubicInterpolant(ts, p.(ts) + p.(-ts), pp.(ts) - pp.(-ts))
    @constraint(model, pspline1err >= h^4 / 384 * (2 * n^4 * (2^(n-4))))
    constrain_interpolant_nonnegative!(model, pspline1 - pspline1err)

    # Now we constrain that ||x| - (p(x) - p(-x)) * x / (p(x) + p(-x))| <= t
    t = 3 * exp(-sqrt(n))
    # Need t * (p(x) + p(-x)) - |x| *(p(x) + p(-x)) + x * (p(x) - p(-x)) >= 0
    @variable(model, pspline2err)
    @constraint(model, pspline2err >= h^4 / 384 * ((2 * t + 4) * n^4 * (2^(n-4))+4*4*n^(4-1)*2^(n-3)))
    f2p = (x -> if x>= 0 
        -(t - 2* x) * pp(-x) + t * pp(x) - 2*p(-x)
    else 
        (-t * pp(-x) + (t + 2 * x) * pp(x) + 2 * p(x)) end)
    f3p = x -> (-abs(x) * pp(-x) + abs(x)* pp(x) + (sign(x) *p(-x))+
                (sign(x) *p(x)) - t* pp(-x) + t* pp(x) - x *pp(-x) 
                - x *pp(x) + p(-x) - p(x) )
    pspline2 = CubicInterpolant(ts, (t * (p.(ts) + p.(-ts)) - abs.(ts) .* (p.(ts) + p.(-ts)) 
                                    + ts .* (p.(ts) - p.(-ts))), 
                                    f2p.(ts))
    pspline3 = CubicInterpolant(ts, t * (p.(ts) + p.(-ts)) + abs.(ts) .* (p.(ts) + p.(-ts)) 
                                - ts .* (p.(ts) - p.(-ts)),
                                f3p.(ts))
    constrain_interpolant_nonnegative!(model, pspline2 - pspline2err)
    constrain_interpolant_nonnegative!(model, pspline3 - pspline2err)
    set_silent(model)
    optimize!(model)
    return primal_status(model) == MathOptInterface.FEASIBLE_POINT
end

"""
    certify_newman_bound_lp(n)
Certify the order N newman bound for rational approximation to |x|
with a linear program
"""
function certify_newman_bound_lp(n, h)
    ts = -1:h:1
    @assert 0 in ts "We need 0 in ts for |x| to be differentiable on each subinterval"
    model = Model()
    set_optimizer(model, () -> Mosek.Optimizer(
        MSK_IPAR_INTPNT_SOLVE_FORM = MSK_SOLVE_DUAL))
    @variable(model, pspline1err)
    h = maximum(diff(ts))
    z = exp(-n^(-1/2))
    xvar = Polynomials.Polynomial([0,1])
    p = prod([xvar + z^k for k = 0:(n-1)])
    # Constraint that p(x) + p(-x) >= 0
    @constraint(model, pspline1err >= h^2 / 8 * (2 * n^2 * (2^(n-2))))
    @constraint(model, p.(ts) + p.(-ts) .>= pspline1err)

    # Now we constrain that ||x| - (p(x) - p(-x)) * x / (p(x) + p(-x))| <= t
    t = 3 * exp(-sqrt(n))
    # Need t * (p(x) + p(-x)) - |x| *(p(x) + p(-x)) + x * (p(x) - p(-x)) >= 0
    @variable(model, pspline2err)
    @constraint(model, pspline2err >= h^2 / 8 * ((2 * t + 4) * n^2 * (2^(n-2))+ 4*2 * n^1 * 2^(n-1)))
    @constraint(model, (t * (p.(ts) + p.(-ts)) - abs.(ts) .* (p.(ts) + p.(-ts)) 
    + ts .* (p.(ts) - p.(-ts))) .>= pspline2err)
    @constraint(model, (t * (p.(ts) + p.(-ts)) + abs.(ts) .* (p.(ts) + p.(-ts)) 
    - ts .* (p.(ts) - p.(-ts))) .>= pspline2err)
    set_silent(model)
    optimize!(model)
    return primal_status(model) == MathOptInterface.FEASIBLE_POINT
end


"""
get_num_required_point(n, eval_func)
Find the number of required points in the discretization
to certify the Newman bound
"""
function get_num_required_point(n, eval_func)
    # First find an upper bound on what's necessary
    estimated = 2
    certifies = false
    while ~certifies
        estimated = 10 * estimated - 1
        certifies = eval_func(n, 2/(estimated - 1))
    end
    # now binary search in between
    large_enough = estimated
    too_small = 1
    while large_enough - too_small > 2
        in_between = Int(floor((too_small + large_enough)/2))+1
        odd_in_between = 2 * Int(floor(in_between/2)) - 1
        certifies = eval_func(n, 2/(odd_in_between - 1))
        if certifies
            large_enough = odd_in_between
        else
            too_small = odd_in_between
        end
    end
    return large_enough
end

"""
    optimize_newman_bound(n, h)
Find the best bound with this n and h
"""
function optimize_newman_bound(n, h)
    ts = -1:h:1
    @assert 0 in ts "We need 0 in ts for |x| to be differentiable on each subinterval"
    model = Model()
    set_optimizer(model, () -> Mosek.Optimizer(
        MSK_IPAR_INTPNT_SOLVE_FORM = MSK_SOLVE_DUAL))
    @variable(model, pspline1err)
    h = maximum(diff(ts))
    z = exp(-n^(-1/2))
    xvar = Polynomials.Polynomial([0,1])
    p = prod([xvar + z^k for k = 0:(n-1)])
    pp = derivative(p, 1)
    # Constraint that p(x) + p(-x) >= 0
    pspline1 = CubicInterpolant(ts, p.(ts) + p.(-ts), pp.(ts) - pp.(-ts))
    @constraint(model, pspline1err >= h^4 / 384 * (2 * n^4 * (2^(n-4))))
    constrain_interpolant_nonnegative!(model, pspline1 - pspline1err)

    # Now we constrain that ||x| - (p(x) - p(-x)) * x / (p(x) + p(-x))| <= t
    @variable(model, t) #t = 3 * exp(-sqrt(n))
    # Need t * (p(x) + p(-x)) - |x| *(p(x) + p(-x)) + x * (p(x) - p(-x)) >= 0
    @variable(model, pspline2err)
    @constraint(model, pspline2err >= h^4 / 384 * ((2 * t + 4) * n^4 * (2^(n-4))+4*4*n^(4-1)*2^(n-3)))
    f2p = (x -> if x>= 0 
        -(t - 2* x) * pp(-x) + t * pp(x) - 2*p(-x)
    else 
        (-t * pp(-x) + (t + 2 * x) * pp(x) + 2 * p(x)) end)
    f3p = x -> (-abs(x) * pp(-x) + abs(x)* pp(x) + (sign(x) *p(-x))+
                (sign(x) *p(x)) - t* pp(-x) + t* pp(x) - x *pp(-x) 
                - x *pp(x) + p(-x) - p(x) )
    pspline2 = CubicInterpolant(ts, (t * (p.(ts) + p.(-ts)) - abs.(ts) .* (p.(ts) + p.(-ts)) 
                                    + ts .* (p.(ts) - p.(-ts))), 
                                    f2p.(ts))
    pspline3 = CubicInterpolant(ts, t * (p.(ts) + p.(-ts)) + abs.(ts) .* (p.(ts) + p.(-ts)) 
                                - ts .* (p.(ts) - p.(-ts)),
                                f3p.(ts))
    constrain_interpolant_nonnegative!(model, pspline2 - pspline2err)
    constrain_interpolant_nonnegative!(model, pspline3 - pspline2err)
    set_silent(model)
    @objective(model, Min, t)
    optimize!(model)
    if primal_status(model) != MathOptInterface.FEASIBLE_POINT
        return Inf
    end
    return value(t)
end

"""
    optimize_newman_bound_lp(n, h)
Find best newman bound with LP N newman bound for rational approximation to |x|
with a linear program
"""
function optimize_newman_bound_lp(n, h)
    ts = -1:h:1
    @assert 0 in ts "We need 0 in ts for |x| to be differentiable on each subinterval"
    model = Model()
    set_optimizer(model, () -> Mosek.Optimizer(
        MSK_IPAR_INTPNT_SOLVE_FORM = MSK_SOLVE_DUAL))
    @variable(model, pspline1err)
    h = maximum(diff(ts))
    z = exp(-n^(-1/2))
    xvar = Polynomials.Polynomial([0,1])
    p = prod([xvar + z^k for k = 0:(n-1)])
    # Constraint that p(x) + p(-x) >= 0
    @constraint(model, pspline1err >= h^2 / 8 * (2 * n^2 * (2^(n-2))))
    @constraint(model, p.(ts) + p.(-ts) .>= pspline1err)

    # Now we constrain that ||x| - (p(x) - p(-x)) * x / (p(x) + p(-x))| <= t
    @variable(model, t) #t = 3 * exp(-sqrt(n))
    # Need t * (p(x) + p(-x)) - |x| *(p(x) + p(-x)) + x * (p(x) - p(-x)) >= 0
    @variable(model, pspline2err)
    @constraint(model, pspline2err >= h^2 / 8 * ((2 * t + 4) * n^2 * (2^(n-2))+ 4*2 * n^1 * 2^(n-1)))
    @constraint(model, (t * (p.(ts) + p.(-ts)) - abs.(ts) .* (p.(ts) + p.(-ts)) 
    + ts .* (p.(ts) - p.(-ts))) .>= pspline2err)
    @constraint(model, (t * (p.(ts) + p.(-ts)) + abs.(ts) .* (p.(ts) + p.(-ts)) 
    - ts .* (p.(ts) - p.(-ts))) .>= pspline2err)
    set_silent(model)
    @objective(model, Min, t)
    optimize!(model)
    if primal_status(model) != MathOptInterface.FEASIBLE_POINT
        return Inf
    end
    return value(t)
end



"""
    find_rational_approx(n, h, t)
Find rational approx within t error (which we search with binary search) with this degree n and h
Use cheb basis for both numerator and denominator
If verbose > 1, print out all SOCP solve information
"""
function find_rational_approx(n, h, t; verbose = 0)
    ts = -1:h:1
    C42B = cheb2bernstein(n, 4)
    C32B = cheb2bernstein(n, 3)
    @assert 0 in ts "We need 0 in ts for |x| to be differentiable on each subinterval"
    model = Model()
    set_optimizer(model, () -> Mosek.Optimizer(
        MSK_IPAR_INTPNT_SOLVE_FORM = MSK_SOLVE_DUAL))
    @variable(model, pcoefs[1:n])
    @variable(model, qcoefs[1:n])
    # establish a scale:
    @constraint(model, qcoefs[1] == 1)
    @variable(model, qsplineerr)
    h = maximum(diff(ts))
    qspline = CubicInterpolant(ts, eval_cheb.(Ref(qcoefs), ts, Ref(0)),
        eval_cheb.(Ref(qcoefs), ts, Ref(1)))
    @variable(model, q4)
    @constraint(model, q4 .>= C42B * qcoefs)
    @constraint(model, q4 .>= -C42B * qcoefs)
    @constraint(model, qsplineerr >= h^4 / 384 * q4)
    constrain_interpolant_nonnegative!(model, qspline - qsplineerr)

    # Now we constrain that ||x| - p(x) / q(x)| <= t
    # Need tq - p + |x|q >= 0 and  tq+ p - |x|q >= 0
    @variable(model, pspline2err)
    pfunc = x -> eval_cheb(pcoefs, x, 0)
    ppfunc = x -> eval_cheb(pcoefs, x, 1)
    qfunc = x -> eval_cheb(qcoefs, x, 0)
    qpfunc = x -> eval_cheb(qcoefs, x, 1)
    # Let g =|x|*q, compute the derivative
    gp = x -> (if x >= 0 qfunc(x) + x*qpfunc(x) else -qfunc(x) - x*qpfunc(x) end)
    # |x| < 1, so we can bound |x| by 1. (x * q)^(4) = q^(4) + xq^(3)
    @variable(model, q3)
    @constraint(model, q3 .>= C32B * qcoefs)
    @constraint(model, q3 .>= -C32B * qcoefs)
    @variable(model, p4)
    @constraint(model, p4 .>= C42B * pcoefs)
    @constraint(model, p4 .>= -C42B * pcoefs)
    pspline2 = CubicInterpolant(ts, t * qfunc.(ts) - pfunc.(ts) + abs.(ts) .* qfunc.(ts), 
                                    t * qpfunc.(ts) - ppfunc.(ts) + gp.(ts))
    pspline3 = CubicInterpolant(ts, t * qfunc.(ts) + pfunc.(ts) - abs.(ts) .* qfunc.(ts), 
                                    t * qpfunc.(ts) + ppfunc.(ts) - gp.(ts))
    @constraint(model, pspline2err >= h^4 / 384 * (q4 + q3 + p4 + t * q4))
    constrain_interpolant_nonnegative!(model, pspline2 - pspline2err)
    constrain_interpolant_nonnegative!(model, pspline3 - pspline2err)
    if verbose <= 1
        set_silent(model)
    end
    optimize!(model)
    return (primal_status(model) == MathOptInterface.FEASIBLE_POINT, value.(pcoefs), value.(qcoefs))
end

"""
find_best_rational_approx(n, h)
Use binary search to find the best possible error t
verbose = 1 prints out the error at each step, = 2 prints all SOCP solve information
"""
function find_best_rational_approx(n, h; verbose = 0)
    t_best = t_curr = 1;
    solvable = true
    while solvable
        t_curr = t_best / 2;
        if verbose>0 print("testing t = $t_curr\n") end
        solvable, p, q = find_rational_approx(n, h, t_curr, verbose = verbose)
        if solvable
            t_best = t_curr
        end
    end
    # now search between t_curr and t_best
    ub = t_best
    lb = t_curr
    while ub - lb > 1e-8
        t_curr = (ub + lb) / 2
        if verbose>0 print("testing t = $t_curr\n") end
        solvable, p, q = find_rational_approx(n, h, t_curr, verbose = verbose)
        if solvable
            ub = t_curr
        else
            lb = t_curr
        end
    end
    return t_curr
end