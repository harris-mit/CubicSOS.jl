# These are the main tools useful in running tests for the sphere packing application

function solve_sphere_packing_socp(rad, Fcoefs, Fpcoefs, F4bnd, xs, ys)
    # Newest solver designed to use the CubicSpline object to represent splines.
    num_xs = length(xs)
    num_ys = length(ys)

    model = Model()
    set_optimizer(model, () -> Mosek.Optimizer(
        MSK_IPAR_INTPNT_SOLVE_FORM = MSK_SOLVE_DUAL,
        #MSK_DPAR_INTPNT_CO_TOL_DFEAS = 1e-10,
        #MSK_DPAR_INTPNT_CO_TOL_PFEAS = 1e-10,
        #MSK_IPAR_PRESOLVE_USE = 0, # MSK_PRESOLVE_MODE_FREE
        MSK_IPAR_NUM_THREADS = 2 ))

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
    constrain_spline_nonnegative!(model, fhat - delta)

    # Constrain that f(0) == fhat(0) and f(0) > 0
    @constraint(model, f.y_vals[1] == fhat.y_vals[1])
    @constraint(model, f.y_vals[1] >= 1)

    # MOI.set(model, MOI.Silent(), true) # This silences the solver.
    optimize!(model)
    # MOI.write_to_file(SumOfSquares.backend(model).optimizer.model, "n8r15.task.gz")
    # This line^ dumps the problem to a file.
    @show termination_status(model)
    @show primal_status(model)
    @show value(l1_F4bnd)
    @show value(delta)
    return f, fhat, delta
end


"""
    get_opt_radius(n)
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

"""
    check_feasibility(rad, f, fhat, xs, ys, delta)
Check the sign constraints on f and fhat.
Does not check the Fourier transform is correct.
"""
function check_feasibility(rad, f, fhat, xs, ys)
    num_checks = 1000
    x_check_pts = rad:(maximum(xs) - minimum(xs))/(num_checks-1):maximum(xs)
    y_check_pts = minimum(ys):(maximum(ys) - minimum(ys))/(num_checks-1):maximum(ys)
    fmax = maximum(value.(evaluate_cubic.(Ref(f),x_check_pts)) .* exp.(-pi .* x_check_pts.^2)) # want this to be < 0
    fhatmin = minimum(value.(evaluate_cubic.(Ref(fhat), y_check_pts)))
    return fmax, fhatmin
end


function solve_delsarte_socp(Amax, Gcoefs, xs)
    # Find a cubic spline function and find its Gegenbauer expansion
    num_xs = length(xs)

    model = Model()
    set_optimizer(model, () -> Mosek.Optimizer(
        MSK_IPAR_INTPNT_SOLVE_FORM = MSK_SOLVE_DUAL))

    # Set up function f
    @variable(model, f_vals[1:num_xs])
    @variable(model, f_derivs[1:num_xs])
    f = CubicSpline(xs, f_vals, f_derivs)

    for i = 1:num_xs - 1
        if xs[i + 1] <= Amax
            constrain_spline_nonnegative!(model, -f, i)
        end
    end

    C = get_hermite_basis_coefficients(f) # Arrange the spline coefficients in a nice way
    #f_0 = 1
    @constraint(model, sum(Gcoefs[1, :, :] .* C) == 1)
    # f_k >= 0 for all k (truncated at some point)
    for k = 2:size(Gcoefs, 1)
        @constraint(model, sum(Gcoefs[k, :, :] .* C) >= 0)
    end

    @objective(model, Min, f.y_vals[end])

    # MOI.set(model, MOI.Silent(), true) # This silences the solver.
    optimize!(model)
    @show termination_status(model)
    @show primal_status(model)
    return f
end

function solve_delsarte_socp2(d, kmax, xs)
    # Find Gegenbauer coefficients so that
    # spline of time domain is provably nonpositive on [-1, Amax]
    # !!! Here we assume xs only goes to Amax
    model = Model()
    set_optimizer(model, () -> Mosek.Optimizer(
        MSK_IPAR_INTPNT_SOLVE_FORM = MSK_SOLVE_DUAL))

    # Set up fk coefficients
    @variable(model, fk[1:kmax] >= 0)
    @constraint(model, fk[1] == 1) # constant coefficient is 1

    f = CubicSpline(xs,
            compute_from_expansion.(Ref(d), Ref(fk), xs),
            compute_deriv_from_expansion.(Ref(d), Ref(fk), xs)
            )
    @variable(model, delta)
    deltax = xs[2] - xs[1] # assuming even spacing
    @constraint(model, delta >= compute_4th_div_bound(d, fk) * deltax^4 / 384)
    #@constraint(model, delta >= deltax^4 / 384)
    constrain_spline_nonnegative!(model, -f - delta)
    @objective(model, Min, compute_from_expansion(d, fk, 1))

    # MOI.set(model, MOI.Silent(), true) # This silences the solver.
    optimize!(model)
    @show termination_status(model)
    @show primal_status(model)
    @show value(delta)
    return value.(fk)
end
