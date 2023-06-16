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
    return value.(x)
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
    return value.(x)
end

# Now we repeat the above but for bump functions
"""
    get_bump_approx(ts, f, fp, fcnf4thdiv, bumps, bumpsderivs)
To approximate a function f with the best estimate in the
sup norm over [ts[1], ts[end]]. Derivative of function is fp. fcnfthdiv is a bound
on the fourth derivative of f. ts are the points used for the spline.
Bumps are a list of the bump functions we should use.
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
    # Each bumps' 4th derivative is bounded by 8315.9
    @constraint(model, [chebfrthdivbound; x .* 8315.9]
                        in MathOptInterface.NormOneCone(length(x)+1))
    @constraint(model, delta >= deltax^4 / 384 * (fcnf4thdiv + chebfrthdivbound))
    constrain_interpolant_nonnegative!(model, fspline - pspline + err - delta)
    constrain_interpolant_nonnegative!(model, pspline - fspline + err - delta)
    optimize!(model)
    return value.(x)
end


"""
    get_bump_approx_lp(ts, f, fcnsecdiv, bumps)
To approximate a function f with the best estimate in the
sup norm over [ts[1], ts[end]]. Derivative of function is fp. fcnsecdiv is a bound
on the second derivative of f. ts are the points used for linear program.
"""
function get_bump_approx_lp(ts, f, fcnsecdiv, bumps)
    model = Model()
    set_optimizer(model, () -> Mosek.Optimizer(
        MSK_IPAR_INTPNT_SOLVE_FORM = MSK_SOLVE_DUAL))
    @variable(model, x[1:length(bumps)])
    @variable(model, err)
    @objective(model, Min, err)
    @variable(model, chebsecdivbound)
    @variable(model, delta)
    deltax = maximum(diff(ts))
    #The below bound is more naive -- just the l1 norm of the coefficients.
    # Each bumps second derivative is bounded by 7.75
    @constraint(model, [chebsecdivbound; x .* 7.75]
                        in MathOptInterface.NormOneCone(length(x)+1))
    @constraint(model, delta >= deltax^2 / 8 * (fcnsecdiv + chebsecdivbound))
    fvals = f.(ts)
    pvals = hcat([b.(ts) for b in bumps]...)*x
    @constraint(model, fvals - pvals .+ err .>= delta)
    @constraint(model, pvals - fvals .+ err .>= delta)
    optimize!(model)
    return value.(x)
end