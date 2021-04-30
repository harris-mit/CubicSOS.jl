# This file has the helper functions that constrain interpolating
# cubic polynomials to be SOS on an interval.

struct CubicSpline5
    x_vals::AbstractArray # The interpolation x points
    y_vals::Array{VariableRef, 1} # Value at grid points
    deriv_vals::Array{VariableRef, 1} # The derivative value at grid points
    function CubicSpline5(x_vals::AbstractArray, y_vals::Array{VariableRef, 1},
        deriv_vals::Array{VariableRef, 1})
        if length(x_vals) != length(y_vals)
            error("Invalid number of function values")
        end
        if length(x_vals) != length(deriv_vals)
            error("Invalid number of derivative values")
        end
        new(x_vals, y_vals, deriv_vals)
    end
end
CubicSpline = CubicSpline5

"""
    x2idx(cs, x)
Get the index that a point x belongs to.
"""
function x2idx(cs::CubicSpline, x)
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
    get_hermite_basis_coefficients(cs::CubicSpline, int_num)
Returns the vector of coefficient variables in the order corresponding
to the Hermite basis.
"""
function get_hermite_basis_coefficients(cs::CubicSpline, int_num)
    x = [cs.y_vals[int_num], cs.deriv_vals[int_num], cs.y_vals[int_num + 1], cs.deriv_vals[int_num + 1]]
    return x
end

"""
    evaluate_cubic(cs::CubicSpline, x)
Evaluates a cubic spline that has a value at a point
"""
function evaluate_cubic(cs::CubicSpline, x)
    idx = x2idx(cs, x)
    xi = cs.x_vals[idx]
    xip1 = cs.x_vals[idx + 1]
    get_one_base_element = j -> get_hermite_basis(x, j, xi, xip1)
    basis_elements = get_one_base_element.(1:4)
    hermite_basis_coefficients = get_hermite_basis_coefficients(cs, idx)
    return sum(basis_elements .* hermite_basis_coefficients)
end

"""
    constrain_spline_nonnegative!(model::AbstractModel, cs::CubicSpline, interval_number::Number)

If p is the spline, we add the constraint p >= 0 on the interval given by interval_number
to model.
"""
function constrain_spline_nonnegative!(model::AbstractModel, cs::CubicSpline, interval_number::Number)
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
            == Q1[2] + 2 * Q1[1]))
    @constraint(model, (cs.deriv_vals[interval_number + 1] * delta - 3 * cs.y_vals[interval_number + 1]
            == - 2 * Q1[1] - Q2[1]))

    constrain_2x2_psd!(model, [cs.y_vals[interval_number + 1], Q1[1], Q1[2]])
    constrain_2x2_psd!(model, [Q2[1], Q2[2], cs.y_vals[interval_number]])
end


"""
    constrain_2x2_psd!(model::AbstractModel, x)
Adds an SOCP constraint that the matrix
x[1] x[2]
x[2] x[3]
is positive semidefinite.
"""
function constrain_2x2_psd!(model::AbstractModel, x)
    @constraint(model, [(x[1] + x[3])/2, (x[1] - x[3])/2, x[2]] in SecondOrderCone())
end


"""
    length(cs::CubicSpline)
Returns the number of points used in a spline.
"""
Base.length(cs::CubicSpline) = length(cs.x_vals) # Overload


#################################################

struct CubicSpline4
    x_vals::AbstractArray # interpolation x values
    Q1s::Convex.Variable # stores information from the first matrix in spline
    Q2s::Convex.Variable # stores information from second matrix
    function CubicSpline4(x_vals::AbstractArray, Q1s::Convex.Variable,
        Q2s::Convex.Variable)
        if length(Q1s) != 3 * length(x_vals) - 3
            error("Invalid length of Q1")
        end
        if length(Q2s) != 1 + length(x_vals)
            error("Invalid length of Q2")
        end
        # Q1s = Convex.Variable(3 * (length(x_vals) - 1))
        # Q2s = Convex.Variable(1 + length(x_vals))
        # Aside from first matrix,
        # two parameters from Q2 are redundant with the continuity conditions
        return new(x_vals, Q1s, Q2s)
    end
end

CubicSpline = CubicSpline4

"""
    get_matrix(cs::CubicSpline, interval_index, matrix_index)

matrix_index: either the first or the second index
Returns the three variables in the vector x where the matrix is
x[1] x[2]
x[2] x[3]
"""
function get_matrix(cs::CubicSpline, interval_index, matrix_index)
    if interval_index >= length(x_vals)
        error("Index out of bounds")
    end
    if !(matrix_index in [1,2])
        error("Matrix index out of bounds")
    end
    if matrix_index == 1
        return cs.Q1s[3 * interval_index - 2: 3 * interval_index]
    else
        # matrix index == 2, the first element is stored directly when
        # it's the first matrix.
        if interval_index == 1
            return cs.Q2s[1:3]
        end
        # otherwise, we only store Q2s[1,1] in this array, and we have that
        # Q2s[i+1][2,2] = Q1s[i][1,1] by matching the value of the polynomial and
        # for the derivative
        # 1/(x_i - x_{i-1}) * (3 * Q1[i-1][1,1] - 2 * Q1[i-1][1,2] - Q2[i-1][1,1] ==
        # 1/(x_{i+1} - x_i) * (Q1[i][2,2] - 3 * Q2[i][2,2]  + 2 * Q2[i][1,2])
        x = Array{Convex.AbstractExprOrValue}(undef, 3)
        x[1] = cs.Q2s[interval_index + 2] # This is Q2s[interval_index][1,1] stored exactly
        prev_Q1 = get_matrix(cs, interval_index - 1, 1)
        prev_Q211 = cs.Q2s[interval_index + 1]
        fp = (1/(cs.x_vals[interval_index] - cs.x_vals[interval_index - 1]) *
                (3 * prev_Q1[1] - 2 * prev_Q1[2] - prev_Q211))
        # = F'(x_vals[interval_index])
        x[2] = .5 * ((cs.x_vals[interval_index + 1] - cs.x_vals[interval_index]) * fp
                        - cs.Q1s[3 * interval_index] + 3 * prev_Q1[1])
        #TODO CHECK SCALINGGGG!!!! AND EQUATIONS
        x[3] = prev_Q1[1]
        println("running again")
        @show x[1]
        @show x[2]
        @show x[3]
        @show interval_index
        return x
    end
end

"""
Gets the values and derivatives at the interpolation points.
"""
function get_vals_and_derivs(cs::CubicSpline)
    vals = Array{Convex.AbstractExprOrValue, 2}(undef, length(cs), 2)
    # First column is value, second is derivative
    # For the first point we use the left hand point formula.
    vals[1,1] = get_matrix(cs, 1, 2)[3] #p(x_1) = Q_2^{22} for the first interval
    deli = cs.x_vals[2] - cs.x_vals[1]
    vals[1,2] = 1/deli * (get_matrix(cs, 1, 1)[3] - 3 * get_matrix(cs, 1, 2)[3]
                            + 2 * get_matrix(cs, 1, 2)[2])

    # It is more efficient from our data storage method to access the right hand endpoint
    for int_num = 1:(length(cs) - 1) # compute p(x_{int_num +1})
        # The value at the point is p(x_{int_num+1}) = Q_1^{11}
        vals[int_num+1, 1] = get_matrix(cs, int_num, 1)[1]
        deli = cs.x_vals[int_num + 1] - cs.x_vals[int_num]
        vals[int_num+1, 2] = 1/deli*(3 * get_matrix(cs, int_num, 1)[1]
                                    - 2 * get_matrix(cs, int_num, 1)[2] +
                                    - get_matrix(cs, int_num, 2)[1])
    end
    return vals
end

"""
Gets the Lagrange-inspired basis coefficients.
"""
function get_poly_coefficients(cs::CubicSpline)
    coefs = Array{Convex.AbstractExprOrValue, 2}(undef, length(cs) - 1, 4)
    for int_num = 1:(length(cs) - 1)
        Q1 = get_matrix(cs, int_num, 1)
        Q2 = get_matrix(cs, int_num, 2)
        coefs[int_num, 1] = Q1[1]
        coefs[int_num, 2] = -2 * Q1[2] - Q2[1]
        coefs[int_num, 3] = 2 * Q2[2] + Q1[3]
        coefs[int_num, 4] = -Q2[3]
    end
    return coefs
end

function get_lagrange_basis(x, basis_index, xip1, xi)
    delta3 = (xip1 - xi)^3
    val = 0
    if basis_index == 1
        val = (x - xi)^3 / delta3
    elseif basis_index == 2
        val = (x - xi)^2 * (x - xip1)/delta3
    elseif basis_index == 3
        val = (x - xi) * (x - xip1)^2 / delta3
    elseif basis_index == 4
        val = (x - xip1)^3 / delta3
    else
        error("Basis index for evaluation is invalid. Must be integer 1 to 4.")
    end
    return val
end




"""
    evaluate_cubic(cs::cubicSpline, xx)
Evaluate the cubic spline at the points xx
"""
function evaluate_cubic(cs::CubicSpline, xx)
    coefs = get_poly_coefficients(cs)
    coefs = evaluate.(coefs)
    coef_vals = zeros(size(coefs))
    for i = 1:size(coefs,1)
        for j = 1:size(coefs,2)
            coef_vals[i,j] = coefs[i,j][1]
        end
    end
    evaluation_vals = zeros(length(xx))
    for i = 1:length(xx)
        int_num = x2idx(cs, xx[i])
        xi = cs.x_vals[int_num]
        xip1 = cs.x_vals[int_num + 1]
        get_basis_elements = j -> get_lagrange_basis(xx[i], j, xip1, xi)
        basis_vec = get_basis_elements.(1:4)
        evaluation_vals[i] = sum(coef_vals[int_num, :] .* basis_vec)
    end
    return evaluation_vals
end

"""
Add a constraint to a Convex.jl problem that the cubics on an interval
are nonnegative by constraining these matrices to be PSD.
"""
function add_nonnegative_constraint!(problem, interval_index, cs::CubicSpline)
    Q1 = get_matrix(cs::CubicSpline, interval_index, 1)
    Q2 = get_matrix(cs::CubicSpline, interval_index, 2)
    add_psd_constraint!(problem, Q1)
    add_psd_constraint!(problem, Q2)
end

"""
Add a PSD constraint when a symmetric matrix variable is given as
X1 X2
XX X3
"""
function add_psd_constraint!(problem, X)
    problem.constraints += (X[1] + X[3])/2 >= norm([X[2];(X[1] - X[3])/2],2)
end


struct CubicSplineInterval1 # A single interval's cubic spline
    x1::Real
    x2::Real
    Q1::Convex.ConvexVariable
    Q2::Convex.ConvexVariable


function getCubicSplineVals(cs::CubicSpline)


struct CubicSplineVals1 <: CubicSpline
    x_vals::AbstractArray # interpolation points
    function_vals::AbstractArray # evaluation values
    deriv_vals::AbstractArray # derivative evaluation values

    function CubicSplineVals1(x_vals::AbstractArray, function_vals::AbstractArray,
        deriv_vals::AbstractArray)
        if length(x_vals) != length(function_vals)
            error("Invalid number of values for splining")
        elseif length(x_vals) != length(deriv_vals)
            error("Invalid number of derivative values for splining")
        end
        return new(x_vals, function_vals, deriv_vals)
    end
end

# The matrix values, these are Convex.jl Variables
struct CubicSplineMats1 <: CubicSpline
    xs::AbstractArray
    Q1s::Array{Convex.AbstractExprOrValue}
    Q2s::Array{Convex.AbstractExprOrValue}
end


CubicSpline = CubicSpline6

# For each spline, write p as a quadratic form of Q1 and Q2
# Q1s[i] is first matrix for the ith interval

# A final representation is in the basis described in the reference



# overload basic operators: Finding lengths, adding, constraining positivity
# Evaluating and plotting.


function get_matrix_representation(cs::CubicSpline)
    Q1s = Array{Convex.AbstractExprOrValue}(undef, length(cs) - 1)
    Q2s = Array{Convex.AbstractExprOrValue}(undef, length(cs) - 1)
    for i = 1:length(cs) - 1
        # These matrices are structured as
        # Q1 Q2
        # XX Q3
        Q1s[i] = Variable(3)
        Q2s[i] = Variable(3)
        Q2s[i][3] = cs.function_vals[i] # Left hand side value
        Q1s[i][1] = cs.function_vals[i+1] # Right hand side value

end

"""
    cs1::CubicSpline + cs2::CubicSpline
Returns a new cubic spline whose values and derivatives are the sum
of the values and derivatives of the functions that cs1 and cs2 interpolate.
"""
function Base.:+(cs1::CubicSpline, cs2::CubicSpline)
    if ~all(cs1.x_vals == cs2.x_vals)
        error("Adding splines require they have same interpolation points")
    end
    return CubicSpline(cs1.x_vals, cs1.function_vals .+ cs2.function_vals,
                        cs1.deriv_vals .+ cs2.deriv_vals)
end

"""
    cs1::CubicSpline - cs2::CubicSpline
Returns a new cubic spline whose values and derivatives are the difference
of the values and derivatives of the functions that cs1 and cs2 interpolate.
"""
function Base.:-(cs1::CubicSpline, cs2::CubicSpline)
    if ~all(cs1.x_vals == cs2.x_vals)
        error("Subtracting splines require they have same interpolation points")
    end
    return CubicSpline(cs1.x_vals, cs1.function_vals .- cs2.function_vals,
                        cs1.deriv_vals .- cs2.deriv_vals)
end


"""
    evaluate(cs::CubicSpline, x)
Evaluates the cubic spline interpolant at a point x
"""
function evaluate(cs::CubicSpline, x)
    if x < minimum(cs.x_vals) || x > maximum(cs.x_vals)
        error("x is out of range of the spline")
    end



function add_positive_constraint!(problem::Convex.Problem, cs::CubicSpline)
# On each subdomain constrain p >= 0

end

"""
Evaluation at arbitrary points
"""
function evaluate_cubic(cs::CubicSpline, xs::Vector{AbstractFloat})

end


"""
    make2x2sdp(problem, X)
For a Convex.jl problem and a 2x2 matrix Variable X, add SOCP constraints to
problem that ensure that X is PSD
"""
function constraint2x2sdptosocp!(problem::Convex.Problem, X)
    #TODO: Check that X is 2x2
    problem.constraints += (X[1,1] + X[2,2])/2 >= norm([X[2,1];(X[1,1] - X[2,2])/2],2)
end



# TODO: Factor out more of the code inside of SpherePacking.jl so we can
# use easily for other purposes...
