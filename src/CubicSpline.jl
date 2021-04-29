# This file has the helper functions that constrain interpolating
# cubic polynomials to be SOS on an interval.

"""
    make2x2sdp(problem, X)
For a Convex.jl problem and a 2x2 matrix Variable X, add SOCP constraints to
problem that ensure that X is PSD
"""
function constraint2x2sdptosocp!(problem, X)
    #TODO: Check that X is 2x2
    problem.constraints += X[1,2] == X[2,1]
    problem.constraints += (X[1,1] + X[2,2])/2 >= norm([X[1,2];(X[1,1] - X[2,2])/2],2)
end

# TODO: Factor out more of the code inside of SpherePacking.jl so we can
# use easily for other purposes...
