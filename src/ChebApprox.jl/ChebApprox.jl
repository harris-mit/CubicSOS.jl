# Implements a polynomial approximation to a given function
# in the Chebyshev or sup norm. We also use a Chebyshev basis
# for fun!


"""
    get_cheb_basis(K)
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
    @constraint(model, [chebfrthdivbound; x .* get_cheb_fourth_div_bound(K)]
                        in MathOptInterface.NormOneCone(length(x)+1))
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
    set_optimizer(model, () -> Mosek.Optimizer(
        MSK_IPAR_INTPNT_SOLVE_FORM = MSK_SOLVE_DUAL))
    @variable(model, x[1:K])
    @variable(model, err)
    @objective(model, Min, err)
    @variable(model, chebsecdivbound)
    @variable(model, delta)
    deltax = maximum(diff(ts))
    @constraint(model, [chebsecdivbound; x .* get_cheb_second_div_bound(K)]
                        in MathOptInterface.NormOneCone(length(x)+1))
    @constraint(model, delta >= deltax^2 / 8 * (fcnsecdiv + chebsecdivbound))
    fvals = f.(ts)
    pvals = eval_cheb.(Ref(x), ts, Ref(0))
    @constraint(model, fvals - pvals .+ err .>= delta)
    @constraint(model, pvals - fvals .+ err .>= delta)
    optimize!(model)
    return value.(x)
end