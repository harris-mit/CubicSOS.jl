using CubicSOS, Test

@time begin
    @time @testset "Cubic Spline unit tests..." begin include("runCubicSpline.jl") end
end
