# These are the main tools useful in running tests for the sphere packing application

function check_bound_feasibility(rad, Fcoefs, Fpcoefs, F4bnd, xs, ys)
    # Newest solver designed to use the CubicSpline object to represent splines.
    num_xs = length(xs)
    num_ys = length(ys)

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
    # @show value(l1_F4bnd)
    # @show value(delta)
    return f, fhat
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
