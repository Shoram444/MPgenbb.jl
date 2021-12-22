using MPgenbb
using Test

@testset "MPgenbb.jl" begin
    # Write your tests here.
    @test getpMag(2, 1) == 2.8284271247461903
    @test getParticleString(3, 2.0, 1.0, 0.2, 0.6) == "3 1.565215413930794 0.31728487213663276 2.3344016402296104 0"
end
