module CubicSOS

using Convex, Mosek, MosekTools
using SumOfSquares, DynamicPolynomials, MathOptInterface
using HCubature, SpecialFunctions

include("interpolatingCubic.jl")
include("SpherePacking/SpherePacking.jl")
include("SpherePacking/PlottingUtils.jl")
include("SpherePacking/FourierTransforms.jl")

end # module
