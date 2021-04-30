# These are the main tools useful in running tests for the sphere packing application

function check_bound_feasibility(rad, Fcoefs, Fpcoefs, F4bnd, xs, ys)
    # Newest solver designed to use the CubicSpline object to represent splines.
    num_xs = length(xs)
    num_ys = length(ys)

    #model = Model(with_optimizer(Mosek.Optimizer, MSK_IPAR_INTPNT_SOLVE_FORM = MSK_SOLVE_DUAL))
    model = Model()
    set_optimizer(model, () -> Mosek.Optimizer(
        MSK_IPAR_INTPNT_SOLVE_FORM = MSK_SOLVE_DUAL))

    # Set up function f
    @variable(model, f_vals[1:num_xs])
    @variable(model, f_derivs[1:num_xs])
    f = CubicSpline(xs, f_vals, f_derivs)

    for i = 1:(num_xs-1)
        if xs[i] >= rad
            constrain_spline_nonnegative!(model, -f, i)
        end
    end

    C = get_hermite_basis_coefficients(f) # Arrange the spline coefficients in a nice way

    # Constrain delta < delta^4 / 384 * ||f^(4)||
    ydelta = ys[2] - ys[1]
    @variable(model, l1_F4bnd)
    @variable(model, delta)
    @constraint(model, [l1_F4bnd; vec(C .* F4bnd)] in MathOptInterface.NormOneCone(length(C)+1))
    @constraint(model, delta >= ydelta^4 / 384 * l1_F4bnd)

    # Create a Spline that is the Fourier transform of the spline of f
    fhat_vals = Array{GenericAffExpr{Float64,VariableRef}, 1}(undef, num_ys)
    fhat_derivs = Array{GenericAffExpr{Float64,VariableRef}, 1}(undef, num_ys)
    for i = 1:num_ys
        fhat_vals[i] = sum(Fcoefs[i, :, :] .* C)
        fhat_derivs[i] = sum(Fpcoefs[i, :, :] .* C)
    end
    fhat = CubicSpline(ys, fhat_vals, fhat_derivs)

    # We need that fhat - delta >= 0 everywhere.
    for i = 1:(num_ys-1)
        constrain_spline_nonnegative!(model, fhat - delta, i)
    end

    # Constrain that f(0) == fhat(0) and f(0) > 0
    @constraint(model, f.y_vals[1] == fhat.y_vals[1])
    @constraint(model, f.y_vals[1] >= 1)

    #MOI.set(model, MOI.Silent(), true)
    #set_optimizer_attribute(model, "MSK_IPAR_INTPNT_SOLVE_FORM", "MSK_SOLVE_DUAL")
    optimize!(model)
    #MOI.write_to_file(SumOfSquares.backend(model).optimizer.model, "n8r15.task.gz")
    @show termination_status(model)
    @show primal_status(model)
    @show value(l1_F4bnd)
    @show value(delta)
    return f, fhat
end

function check_bound_feasibility_cvx(rad, Fcoefs, Fpcoefs, F4bnd, num_xs, num_ys)
    # rad = the radius |x| > rad for which f(x) < 0
    # Other arguments assumed to be computed by populate_integrals!

    problem = satisfy()
    Q1s = []; Q2s = []; X1s = []; X2s = [];
    # Q1s, Q2s, for each interval we have a polynomial p which is written in terms of these matrices
    # X1s, X2s, the matrices defining the interpolants of \hat{f}(y) on each y interval
    num_xintervals = num_xs - 1
    num_yintervals = num_ys - 1
    # Always use lower triangular entry [2,1]
    for idx = 1:num_xintervals
        push!(Q1s, Variable(2,2))
        push!(Q2s, Variable(2,2))
    end
    for idx = 1:num_yintervals
        push!(X1s, Variable(2,2))
        push!(X2s, Variable(2,2))
    end

    # We replace the semidefinite constraint Q >= 0 with the equivalent SOCP
    for idx = 1:num_xintervals
        if xs[idx] >= rad
            constraint2x2sdptosocp!(problem, Q1s[idx])
            constraint2x2sdptosocp!(problem, Q2s[idx])
        end
    end

    # Need X is symmetric positive semidefinite
    for idx = 1:num_yintervals
        constraint2x2sdptosocp!(problem, X1s[idx])
        constraint2x2sdptosocp!(problem, X2s[idx])
    end

    # Now set p to interpolate the values of f. p_i is the polynomial on interval i
    # we need p_i(x_i) = p_{i+1}(x_i) and same for the derivative
    for xind = 2:num_xintervals
        xi = xs[xind]
        xip1 = xs[xind+1]
        xim1 = xs[xind-1]
        problem.constraints += Q1s[xind-1][1,1] == Q2s[xind][2,2]
        # Compute d/dx[p_i(x) exp(-pi*x^2)] via product rule and cancelling gives just
        problem.constraints += ((3 * Q1s[xind-1][1,1] - 2 * Q1s[xind-1][2,1] - Q2s[xind-1][1,1])/(xi - xim1) ==
                         (Q1s[xind][2,2] - 3 * Q2s[xind][2,2] + 2 * Q2s[xind][2,1])/(xip1 - xi))
    end

    # Write out the coefficients of p in the Lagrange-inspired basis
    # used in get_basis_factor of CubicSOS.FourierTransforms.jl
    C = Array{Convex.AbstractExpr, 2}(undef, num_xintervals, 4)
    for intnum = 1:num_xintervals
        C[intnum, 1] = -Q1s[intnum][1,1]
        C[intnum, 2] = -(-2 * Q1s[intnum][2,1] - Q2s[intnum][1,1])
        C[intnum, 3] = -(2 * Q2s[intnum][2,1] + Q1s[intnum][2,2])
        C[intnum, 4] = Q2s[intnum][2,2]
    end
    # We need to use -p to ensure this is less than 0, so all of these
    # coefficients are minus what they are in the paper.


    ygap = ys[2] - ys[1]
    # TODO: We could make this more general, but for the current theory
    # this is not a priority.

    # Finally, this polynomial interpolating the Fourier transform must be nonnegative
    # Constrain the interpolating polynomial to equal these values above

    # Bound the fourth derivative of the interpolant of the fourier transform
    delta = ygap^4 / 384 * sum(abs.(C .* F4bnd))

    for yintnum = 1:num_yintervals
        # Take a look at polynomial yintnum
        problem.constraints += sum(Fcoefs[yintnum, :, :] .* C) - delta >= X2s[yintnum][2,2] #left hand endpoint
        problem.constraints += sum(Fcoefs[yintnum+1, :, :] .* C) - delta >= X1s[yintnum][1,1]
        problem.constraints += (sum(Fpcoefs[yintnum, :, :] .* C) * (ys[yintnum+1]-ys[yintnum])
                == (X1s[yintnum][2,2] - 3 * X2s[yintnum][2,2] + 2 * X2s[yintnum][2,1]))
        problem.constraints += (sum(Fpcoefs[yintnum+1, :, :] .* C) * (ys[yintnum+1]-ys[yintnum])
                == (3 * X1s[yintnum][1,1] - 2 * X1s[yintnum][2,1] - X2s[yintnum][1,1]))
    end

    # We include the constraint \hat{f}(0) == f(0) > 0 noting that xs[1] == ys[1] == 0
    problem.constraints += sum(Fcoefs[1, :, :] .* C) == -Q2s[1][2,2];
    problem.constraints += -Q2s[1][2,2] >= .01;


    opt = () -> Mosek.Optimizer(MSK_IPAR_INTPNT_SOLVE_FORM=MSK_SOLVE_DUAL,
                                MSK_IPAR_PRESOLVE_USE = 0,
                                MSK_DPAR_INTPNT_CO_TOL_PFEAS = 1e-5,
                                MSK_DPAR_INTPNT_CO_TOL_NEAR_REL = 1)
    solve!(problem, opt)
    @show problem.status
    @show MOI.get(problem.model, MOI.PrimalStatus())
   return Q1s, Q2s, X1s, X2s
end


################
# Deprecated!!!
################
""" The following is old code that I'm just pasting here for reference. It should
be deleted once everything above is working properly.
It is more "obviously correct," but it is not clear how to rework it for use
in the monomial basis. Therefore, this should just serve as a reference.
The code below assumes the basis used for the Fourier integral calculations
used a monomial basis, unlike the Lagrange basis we are now using. DO NOT EXPECT
RELIABLE RESULTS FROM THIS CODE!!!!!.
"""
# re-pasting the above code for using the "efficient" integral computations"
function get_monom_basis_vector(xi, xip1)
    b = [1; x; x^2; x^3]
    return b
end
function check_bound_feasibility_with_sos(rad, Fcoefs, Fpcoefs, F4bnd, xs, ys)
    model = SOSModel(Mosek.Optimizer)
    @polyvar x
    num_xs = length(xs)
    num_ys = length(ys)
    num_xintervals = num_xs - 1
    num_yintervals = num_ys - 1
    @variable(model, poly_coefs[1:num_xintervals, 1:4]) # one coefficient for each basis
    @variable(model, delta[1:num_yintervals]) # The error in the 4th derivative estimate
    @variable(model, abscoefs[1:num_xintervals, 1:4]) # Bigger than absolute values of the coefficients
    @variable(model, hat_poly_coefs[1:num_yintervals, 1:4]) # The coefficients of interpolants of Fourier transform
    # enforce negativity
    for interval_num = 1:num_xintervals
        xi = xs[interval_num]
        xip1 = xs[interval_num+1]
        if xi >= rad
            S = @set x >= xi && x <= xip1
            b = get_monom_basis_vector(xi, xip1)
            p = sum(poly_coefs[interval_num, :] .* b)
            @constraint(model, p <= 0, domain = S)
        end
    end
    # enforce continuity of polynomials and derivatives
    for interval_num = 2:(num_xs-1)
        xi = xs[interval_num]
        blhs = get_monom_basis_vector(xs[interval_num-1], xi)
        brhs = get_monom_basis_vector(xi, xs[interval_num+1])
        plhs = sum(poly_coefs[interval_num-1,:] .* blhs)
        prhs = sum(poly_coefs[interval_num,:] .* brhs)
        @constraint(model, plhs(xi) == prhs(xi))
        @constraint(model, differentiate(plhs, x)(xi) == differentiate(prhs, x)(xi))
    end
    # bound delta
    @constraint(model, abscoefs .>= poly_coefs)
    @constraint(model, abscoefs .>= -poly_coefs)
    for interval_num = 1:num_yintervals
        deltay = ys[interval_num+1] - ys[interval_num]
        @constraint(model, delta[interval_num] >= deltay^4 / 384 * sum(abscoefs .* abs.(F4bndeff)))
    end
    # constrain the Fourier transform on each interval
    for interval_num = 1:num_yintervals
        yi = ys[interval_num]
        yip1 = ys[interval_num+1]
        b = get_monom_basis_vector(yi, yip1)
        p = sum(hat_poly_coefs[interval_num, :] .* b)
        @constraint(model, p(yi) == sum(Fcoefs[interval_num, :, :] .* poly_coefs))
        @constraint(model, p(yip1) == sum(Fcoefs[interval_num + 1, :, :] .* poly_coefs))
        @constraint(model, differentiate(p, x)(yi) == sum(Fpcoefs[interval_num, :, :] .* poly_coefs))
        @constraint(model, differentiate(p, x)(yip1) == sum(Fpcoefs[interval_num+1, :, :] .* poly_coefs))
        S = @set x >= yi && x <= yip1
        @constraint(model, p >= delta[interval_num], domain = S) #delta[interval_num]
    end

    # Set f(0) == \hat(f)(0)
    xbasis = get_monom_basis_vector(xs[1], xs[2])
    ybasis = get_monom_basis_vector(ys[1], ys[2])
    f0 = sum(poly_coefs[1,:] .* xbasis)(0)
    fhat0 = sum(hat_poly_coefs[1,:] .* ybasis)(0)
    @constraint(model, fhat0 == 1)
    @constraint(model, f0 == 1)
    MOI.set(model, MOI.Silent(), true)
    optimize!(model)
    #MOI.write_to_file(SumOfSquares.backend(model).optimizer.model, "n8r15.task.gz")
    @show termination_status(model)
    @show primal_status(model)
end


"""
Use the data from Cohn's paper to back out the optimal radius we should be able to find.
Warning!!! Only works for n we have recorded densities for below.
"""
function get_opt_radius(n)
    densities = Dict(
        1 => 1.0,
        2 => 0.906899683,
        3 => 0.779746762,
        4 => 0.647704966,
        5 => 0.524980022,
        6 => 0.417673416,
        7 => 0.327455611,
        8 => 0.253669508
    )
    return (densities[n] / pi^(n/2) * gamma(n/2 + 1))^(1/n) * 2
end
