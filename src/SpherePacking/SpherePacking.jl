
export solve_delsarte_socp

"""
    solve_delsarte_socp(d, kmax, xs, fourth_div_bounds)
Compute the Delsarte bound on a spherical code.
f is required to be nonpositive on [xs[1], xs[end]]
We use 1:kmax Gegenbauer polynomials and return their coefficients so that the
spline interpolant of sum f_k G_k guarantees that the high degree polynomial
is nonpositive on the required interval.
fourth_div_bounds has a bound on |G_k^(4)| for the whole interval
(alternatively we could use a bound on each subinterval)

"""
function solve_delsarte_socp(d, kmax, xs)
    # Find Gegenbauer coefficients so that
    # spline of time domain is provably nonpositive on [-1, Amax]
    # !!! Here we assume xs only goes to Amax
    model = Model()
    set_optimizer(model, () -> Mosek.Optimizer(
        MSK_IPAR_INTPNT_SOLVE_FORM = MSK_SOLVE_DUAL))

    # Set up fk coefficients
    @variable(model, fk[1:kmax] >= 0)
    @constraint(model, fk[1] == 1) # constant coefficient is 1

    f = CubicInterpolant(xs,
            compute_from_expansion.(Ref(d), Ref(fk), xs),
            compute_deriv_from_expansion.(Ref(d), Ref(fk), xs)
            )
    num_xintervals = length(xs) - 1
    # @variable(model, delta[1:num_xintervals])
    # for xinterval = 1:num_xintervals
    #     deltax = xs[xinterval + 1] - xs[xinterval]
    #     @constraint(model, delta[xinterval] >= sum(fk .* fourth_div_bounds[:]) * deltax^4 / 384)
    #     constrain_interpolant_nonnegative!(model, -f - delta[xinterval], xinterval)
    # end
    # This bound computes the global optimum using high degree sum of squares programming
    # For this application, we want to avoid all SDP's for clarity, so leave this commented out.
    # L1 norm of the coefficients is a much weaker bound but works
    # Assuming a single spacing for a single delta
    @variable(model, delta)
    deltax = xs[2] - xs[1] # assuming even spacing
    @constraint(model, delta >= compute_gegenbauer_derivative_bound(d, fk, 4) * deltax^4 / 384)
    constrain_interpolant_nonnegative!(model, -f - delta)
    @objective(model, Min, compute_from_expansion(d, fk, 1))

    optimize!(model)
    @show termination_status(model)
    @show primal_status(model)
    return value.(fk)
end


"""
    solve_delsarte_lp(d, kmax, xs)
Compute the Delsarte bound on a spherical code.
f is required to be nonpositive on [xs[1], xs[end]]
We use 1:kmax Gegenbauer polynomials and return their coefficients so that the
sum f_k G_k is nonpositive on the required interval as
guaranteed by the linear program.
"""
function solve_delsarte_lp(d, kmax, xs)
    # Solves the same problem as solve_delsarte_socp but method is to
    model = Model()
    set_optimizer(model, () -> Mosek.Optimizer(
        MSK_IPAR_INTPNT_SOLVE_FORM = MSK_SOLVE_DUAL))

    # Set up fk coefficients
    @variable(model, fk[1:kmax] >= 0)
    @constraint(model, fk[1] == 1) # constant coefficient is 1

    fvals = compute_from_expansion.(Ref(d), Ref(fk), xs)

    # Assuming a single spacing for a single delta (that is computed more crudely...)
    @variable(model, delta)
    deltax = xs[2] - xs[1] # assuming even spacing
    @constraint(model, delta >= compute_gegenbauer_derivative_bound(d, fk, 2) * deltax^2 / 8)
    # certify just for half of interval -- the other half covered by next point
    @constraint(model, fvals .<= -delta)

    # This next bound computes the global optimum using high degree sum of squares programming
    # For this application, we want to avoid all SDP's for clarity, so leave this commented out.
    # L1 norm of the coefficients is a much weaker bound but works
    # Assuming a single spacing for a single delta
    # num_xintervals = length(xs) - 1
    # @variable(model, delta[1:num_xintervals])
    # for xinterval = 1:num_xintervals
    #     deltax = xs[xinterval + 1] - xs[xinterval]
    #     @constraint(model, delta[xinterval] >= sum(fk .* first_div_bounds[:]) * deltax / 2)
    #     @constraint(model, fvals .<= -delta[xinterval])
    # end
    @objective(model, Min, compute_from_expansion(d, fk, 1))

    optimize!(model)
    @show termination_status(model)
    @show primal_status(model)
    @show value.(delta)
    return value.(fk)
end
