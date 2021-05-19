module CubicSOS

using Convex, Mosek, MosekTools
using JuMP
using SumOfSquares, DynamicPolynomials, MathOptInterface
using HCubature, SpecialFunctions

include("CubicSpline.jl")
include("SpherePacking/SpherePacking.jl")
include("SpherePacking/PlottingUtils.jl")
include("SpherePacking/FourierTransforms.jl")

end # module
