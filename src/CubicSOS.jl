module CubicSOS

using Convex, Mosek, MosekTools
using SCS
using JuMP
using MathOptInterface, LinearAlgebra
using HCubature, SpecialFunctions
using Polynomials, SpecialPolynomials, SumOfSquares, DynamicPolynomials

include("CubicInterpolant.jl")
include("SpherePacking/SpherePacking.jl")
include("SpherePacking/GegenbauerTransforms.jl")
include("FilterDesign/FilterDesign.jl")
include("ChebApprox.jl/ChebApprox.jl")
end # module
