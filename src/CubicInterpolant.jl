# This file has the helper functions that constrain interpolating
# cubic polynomials to be SOS on an interval.

export CubicInterpolant
export constrain_interpolant_nonnegative!
export evaluate_cubic
export evaluate_cubic_derivative, get_hermite_basis, get_hermite_basis_coefficients # undocumented
export constrain_2x2_psd! # for testing only

# We want a type to describe any possible field of the spline variable.
VarOrExpressionOrConstant = Union{VariableRef, GenericAffExpr{T,VariableRef}, T} where T <: Real
VarOrExpressionOrConstantArray = Array{T, 1} where T <: VarOrExpressionOrConstant

struct CubicInterpolant
    x_vals::AbstractArray # The interpolation x points
    y_vals::VarOrExpressionOrConstantArray # Value at grid points
    deriv_vals::VarOrExpressionOrConstantArray # The derivative value at grid points
    function CubicInterpolant(x_vals::AbstractArray, y_vals::VarOrExpressionOrConstantArray,
        deriv_vals::VarOrExpressionOrConstantArray)
        if length(x_vals) != length(y_vals)
            error("Invalid number of function values")
        end
        if length(x_vals) != length(deriv_vals)
            error("Invalid number of derivative values")
        end
        new(x_vals, y_vals, deriv_vals)
    end
end

"""
    x2idx(cs, x)
Get the index that a point x belongs to.
"""
function x2idx(cs::CubicInterpolant, x)
    if x > maximum(cs.x_vals)
        error("x is too large -- out of bounds")
    end
    if x < minimum(cs.x_vals)
        error("x is too small -- out of bounds")
    end
    for j = 1:length(cs.x_vals)
        if x >= cs.x_vals[j] && j == (length(cs.x_vals) - 1) # in last interval
            return j
        elseif x >= cs.x_vals[j] && x < cs.x_vals[j+1]
            return j
        end
    end
end

"""
    get_hermite_basis(x, j, xi, xip1)
Evaluate the jth Hermite basis element for x.
"""
function get_hermite_basis(x, j, xi, xip1)
    val = 0
    delta = xip1 - xi
    t = (x - xi) / delta
    if j == 1
        val = 2 * t^3 - 3 * t^2 + 1
    elseif j == 2
        val = delta * (t^3 - 2 * t^2 + t)
    elseif j == 3
        val = -2 * t^3 + 3 * t^2
    elseif j == 4
        val = delta * (t^3 - t^2)
    else
        error("Invalid basis element requested")
    end
    return val
end

"""
    get_hermite_basis_derivatives(x, j, xi, xip1)
Evaluate the derivative of the jth Hermite basis element for x
with respect to x.
"""
function get_hermite_basis_derivatives(x, j, xi, xip1)
    val = 0
    delta = xip1 - xi
    t = (x - xi) / delta
    if j == 1
        val = 6 * t^2 - 6 * t
    elseif j == 2
        val = delta * (3 * t^2 - 4 * t + 1)
    elseif j == 3
        val = -6 * t^2 + 6 * t
    elseif j == 4
        val = delta * (3 * t^2 - 2 * t)
    else
        error("Invalid basis element requested")
    end
    return val / delta
end

"""
    get_hermite_basis_coefficients_interval(cs::CubicInterpolant, int_num)
Returns the vector of coefficient variables in the order corresponding
to the Hermite basis for a particular interval.
"""
function get_hermite_basis_coefficients_interval(cs::CubicInterpolant, int_num)
    x = [cs.y_vals[int_num], cs.deriv_vals[int_num], cs.y_vals[int_num + 1], cs.deriv_vals[int_num + 1]]
    return x
end


"""
    get_hermite_basis_coefficients(cs::CubicInterpolant)
Returns a matrix such that the ith row, jth column is the coefficient of the jth
polynomial in the ith interval.
"""
function get_hermite_basis_coefficients(cs::CubicInterpolant)
    coefs = Array{GenericAffExpr{Float64,VariableRef},2}(undef, length(cs) - 1, 4)
    for i = 1:length(cs) - 1
        x = get_hermite_basis_coefficients_interval(cs, i)
        for j = 1:4
            coefs[i, j] = x[j]
        end
    end
    return coefs
end


"""
    evaluate_cubic(cs::CubicInterpolant, x::Number)
Evaluates a cubic spline that has a value at a point
"""
function evaluate_cubic(cs::CubicInterpolant, x)
    idx = x2idx(cs, x)
    xi = cs.x_vals[idx]
    xip1 = cs.x_vals[idx + 1]
    get_one_base_element = j -> get_hermite_basis(x, j, xi, xip1)
    basis_elements = get_one_base_element.(1:4)
    hermite_basis_coefficients = get_hermite_basis_coefficients_interval(cs, idx)
    return sum(basis_elements .* hermite_basis_coefficients)
end
"""
    cs(x)
Shorthand for evaluation of a spline cs at point x.
"""
(cs::CubicInterpolant)(x) = evaluate_cubic(cs, x)

"""
    evaluate_cubic_derivative(cs::CubicInterpolant, x)
Evaluates the derivative of a cubic spline that has a value at a point
"""
function evaluate_cubic_derivative(cs::CubicInterpolant, x)
    idx = x2idx(cs, x)
    xi = cs.x_vals[idx]
    xip1 = cs.x_vals[idx + 1]
    get_one_base_element = j -> get_hermite_basis_derivatives(x, j, xi, xip1)
    basis_elements = get_one_base_element.(1:4)
    hermite_basis_coefficients = get_hermite_basis_coefficients_interval(cs, idx)
    return sum(basis_elements .* hermite_basis_coefficients)
end

"""
    constrain_interpolant_nonnegative!(model::AbstractModel, cs::CubicInterpolant, interval_number::Number)

If p is the spline, we add the constraint p >= 0 on the interval given by interval_number
to model.
"""
function constrain_interpolant_nonnegative!(model::AbstractModel, cs::CubicInterpolant, interval_number::Number)
    # For this interval, create the PSD matrices Q1 and Q2:
    if interval_number >= length(cs)
        error("There are at most n-1 intervals. This one is out of bounds")
    end
    Q1 = @variable(model, [1:2])
    Q2 = @variable(model, [1:2])

    # The matrix Q_1 is as follows:
    # p(x_{i+1})  Q1[1]
    # Q1[1]       Q1[2]
    # The matrix Q_2 is:
    # Q2[1]       Q2[2]
    # Q2[2]       p(x_i)
    #
    # We need to add the constraints that
    # p'(x_i) \Delta_i + 3 p(x_i) = Q_1^{22} + 2 * Q_2^{21}
    # p'(x_{i+1}) \Delta_i - 3 p(x_{i+1}) = -2*Q_1^{21} - Q_2^{11}
    # Along with the PSD constraints that Q_1 and Q_2 are PSD.

    delta = cs.x_vals[interval_number + 1] - cs.x_vals[interval_number]
    @constraint(model, (cs.deriv_vals[interval_number] * delta + 3 * cs.y_vals[interval_number]
            == Q1[2] + 2 * Q2[2])) # CHECK THIS CHANGE FROM Q1[1] to Q2[2]!
    @constraint(model, (cs.deriv_vals[interval_number + 1] * delta - 3 * cs.y_vals[interval_number + 1]
            == - 2 * Q1[1] - Q2[1]))

    constrain_2x2_psd!(model, [cs.y_vals[interval_number + 1], Q1[1], Q1[2]])
    constrain_2x2_psd!(model, [Q2[1], Q2[2], cs.y_vals[interval_number]])
end

"""
    constrain_interpolant_nonnegative!(model::AbstractModel, cs::CubicInterpolant, interval_numbers::AbstractArray)

If p is the spline, we add the constraint p >= 0 on the intervals given in interval_numbers
to model.
"""
function constrain_interpolant_nonnegative!(model::AbstractModel, cs::CubicInterpolant, interval_numbers::AbstractArray)
    for i in interval_numbers
        constrain_interpolant_nonnegative!(model, cs, interval_numbers[i])
    end
end

"""
    constrain_interpolant_nonnegative!(model::AbstractModel, cs::CubicInterpolant)

If p is the spline, we add the constraint p >= 0 for all intervals on which p is defined.
"""
function constrain_interpolant_nonnegative!(model::AbstractModel, cs::CubicInterpolant)
    constrain_interpolant_nonnegative!(model, cs, 1:length(cs)-1)
end

"""
    constrain_2x2_psd!(model::AbstractModel, x)
Adds an SOCP constraint that the matrix
x[1] x[2]
x[2] x[3]
is positive semidefinite to model.
"""
function constrain_2x2_psd!(model::AbstractModel, x)
    if length(x) != 3
        error("This constraint formulation only works for 2x2 PSD constraints specified in a particular way.")
    end
    @constraint(model, [(x[1] + x[3])/2, (x[1] - x[3])/2, x[2]] in SecondOrderCone())
end


# Overload some basic operations on cubic splines.

"""
    length(cs::CubicInterpolant)
Returns the number of points used in a spline.
"""
Base.length(cs::CubicInterpolant) = length(cs.x_vals)


"""
    cs1::CubicInterpolant + cs2::CubicInterpolant
Returns a new cubic spline whose values and derivatives are the sum
of the values and derivatives of the functions that cs1 and cs2 interpolate.
If the interpolation points are different, we take their union and evaluate
values and derivatives using the spline.
"""
function Base.:+(cs1::CubicInterpolant, cs2::CubicInterpolant)
    if all(cs1.x_vals == cs2.x_vals) # This is faster if interpolation points are same
        return CubicInterpolant(cs1.x_vals, cs1.y_vals .+ cs2.y_vals,
                            cs1.deriv_vals .+ cs2.deriv_vals)
    end
    x_vals = sort(union(cs1.x_vals, cs2.x_vals))
    return CubicInterpolant(x_vals,
            evaluate_cubic.(Ref(cs1), x_vals) .+ evaluate_cubic.(Ref(cs2), x_vals),
            evaluate_cubic_derivative.(Ref(cs1), x_vals) .+ evaluate_cubic_derivative.(Ref(cs2), x_vals))
end

"""
    cs1::CubicInterpolant - cs2::CubicInterpolant
Returns a new cubic spline whose values and derivatives are the difference
of the values and derivatives of the functions that cs1 and cs2 interpolate.
If the interpolation points are different, we take their union and evaluate
values and derivatives using the spline.
"""
function Base.:-(cs1::CubicInterpolant, cs2::CubicInterpolant)
    if all(cs1.x_vals == cs2.x_vals) # This is faster if interpolation points are same
        return CubicInterpolant(cs1.x_vals, cs1.y_vals .- cs2.y_vals,
                            cs1.deriv_vals .- cs2.deriv_vals)
    end
    x_vals = sort(union(cs1.x_vals, cs2.x_vals))
    return CubicInterpolant(x_vals,
            evaluate_cubic.(Ref(cs1), x_vals) .- evaluate_cubic.(Ref(cs2), x_vals),
            evaluate_cubic_derivative.(Ref(cs1), x_vals) .- evaluate_cubic_derivative.(Ref(cs2), x_vals))
end

"""
    cs::CubicInterpolant + t::VarOrExpressionOrConstant
Returns a new cubic spline whose values & derivs are the
values of the function cs + t. (Derivatives don't change.)
"""
function Base.:+(cs::CubicInterpolant, t::VarOrExpressionOrConstant)
    return CubicInterpolant(cs.x_vals, cs.y_vals .+ t, cs.deriv_vals)
end
function Base.:+(t::VarOrExpressionOrConstant, cs::CubicInterpolant)
    return t + cs
end

"""
    cs::CubicInterpolant - t::VarOrExpressionOrConstant
Returns a new cubic spline whose values & derivs are
the values of the function cs - t. (Derivatives don't change.)
"""
function Base.:-(cs::CubicInterpolant, t::VarOrExpressionOrConstant)
    return CubicInterpolant(cs.x_vals, cs.y_vals .- t, cs.deriv_vals)
end
function Base.:-(t::VarOrExpressionOrConstant, cs::CubicInterpolant)
    return CubicInterpolant(cs.x_vals, t .- cs.y_vals, cs.deriv_vals)
end

"""
    t::VarOrExpressionOrConstant * cs::CubicInterpolant
Returns a new cubic spline whose values and derivatives are
the values of the function t * cs.
A nonconstant should only be used if the entries of the spline are not variables
(otherwise there is a convexity issue).
"""
function Base.:*(cs::CubicInterpolant, t::VarOrExpressionOrConstant)
    return CubicInterpolant(cs.x_vals, cs.y_vals .* t, cs.deriv_vals .* t)
end
function Base.:*(t::VarOrExpressionOrConstant, cs::CubicInterpolant)
    return cs * t
end

"""
    cs1::CubicInterpolant * cs2::CubicInterpolant
Returns a new cubic spline whose values are pointwise products and derivatives
are computed from the product rule of two splines.
At least one of the spline should not be defined by variables or else there will
be convexity issues.
"""
function Base.:*(cs1::CubicInterpolant, cs2::CubicInterpolant)
    if all(cs1.x_vals == cs2.x_vals) # This is faster if interpolation points are same
        return CubicInterpolant(cs1.x_vals, cs1.y_vals .* cs2.y_vals,
                            cs1.deriv_vals .* cs2.y_vals + cs1.y_vals .* cs2.deriv_vals)
    else
        error("Must match domains")
    end
end


"""
    -cs::CubicInterpolant
Returns a new cubic spline whose values and derivatives are
the values of the function -cs.
"""
function Base.:-(cs::CubicInterpolant)
    return -1 * cs
end
