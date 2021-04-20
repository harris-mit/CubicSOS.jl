#__precompile__()

module CubicSOS

using Convex, Mosek, MosekTools
using SumOfSquares, DynamicPolynomials, MathOptInterface
using HCubature, SpecialFunctions

include("SpherePacking/SpherePacking.jl")
include("SpherePacking/PlottingUtils.jl")
include("SpherePacking/FourierTransforms.jl")

end # module
