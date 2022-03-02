# Designing a filter as in Chapter 8 of Goberna
# By Lars Abbe, citing Kortanek and Moulin in Reemtsen & Ruckman


"""
    compute_coding_gain(x, rk)
Returns the coding gain in dB from the a_k given by x for the correlations
defined by ak.
"""
function compute_coding_gain(x, rk)
    obj = rk[1]/2 + sum(x .* rk[2:2:length(rk)])
    # Since rk[1] = 1, the variance of the input signal is assumed to be 1.
    # (PCM band)...
    # This calculation matches the Kortanek and Abbe papers, but their prose description
    # is a little different... This seems to be what they meant.
    coding_gain = 0.5 / sqrt(obj * (1-obj))
    coding_gain_db = 10 * log10(coding_gain)
    return coding_gain_db
end


"""
    evaluate_abbe_g(x, omega)
Function value of Abbe's value of g in Section 6
"""
function evaluate_abbe_g(x, omega)
    val = -1
    for m = 0:length(x)-1
        val = val - 2 * x[m+1] * cos(2 * (2 * m + 1) * pi * omega)
    end
    return val
end

"""
    evaluate_abbe_g_deriv(x, omega)
Derivative value of Abbe's value of g in Section 6
"""
function evaluate_abbe_g_deriv(x, omega)
    val = 0
    for m = 0:length(x)-1
        val = val + (2 * (2 * m + 1) * pi) * 2 * x[m+1] * sin(2 * (2 * m + 1) * pi * omega)
    end
    return val
end


"""
    get_optimal_filter(rk, omegas, N)
Run the SOCP version of the problem described by Kortanek and Moulin
"""
function get_optimal_filter(rk, omegas, N)
    model = Model()
    set_optimizer(model, () -> Mosek.Optimizer(
        MSK_IPAR_INTPNT_SOLVE_FORM = MSK_SOLVE_DUAL))
    @variable(model, x[1:N])
    @objective(model, Min, -rk[1]/2 - sum(x .* rk[2:2:(2*N)]))
    p = CubicInterpolant(omegas, evaluate_abbe_g.(Ref(x), omegas), evaluate_abbe_g.(Ref(x), omegas))
    @variable(model, frthdiv)
    @variable(model, delta)
    @constraint(model, [frthdiv; 2^5 * pi^4 * (2 .* (0:N-1) .+1).^4 .* x] in MathOptInterface.NormOneCone(length(x)+1))
    @constraint(model, delta >= frthdiv * maximum(diff(omegas))^4 / 384)
    constrain_interpolant_nonnegative!(model, -p - delta)
    optimize!(model)
    @show value(frthdiv)
    @show value(delta)
    @show termination_status(model)
    @show primal_status(model)
    return value.(x)
end
