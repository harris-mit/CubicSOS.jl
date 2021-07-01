module CubicSOS

using Convex, Mosek, MosekTools
using JuMP
using MathOptInterface
using HCubature, SpecialFunctions
using Polynomials, SpecialPolynomials

include("CubicSpline.jl")
include("SpherePacking/SpherePacking.jl")
include("SpherePacking/PlottingUtils.jl")
include("SpherePacking/FourierTransforms.jl")
include("SpherePacking/GegenbauerTransforms.jl")
end # module
