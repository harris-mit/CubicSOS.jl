using CubicSOS, Test

@time begin
    @time @testset "Cubic Spline unit tests..." begin include("runCubicSpline.jl") end
    # These transforms are unused, and the tests are somewhat slow, so they're not going to run. Uncomment if desired.
    #@time @testset "Transforms unit tests..." begin include("runTransforms.jl") end
    @time @testset "Delsarte sphere packing bound tests..." begin include("runDelsarte.jl") end
end
