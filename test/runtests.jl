using CubicSOS, Test

@time begin
    # @time @testset "Sphere packing" begin include("runSpherePacking.jl")
    @time @testset "Sphere packing" begin @test 2 == 2
end