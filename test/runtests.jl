using Migmatites
using Test
using CairoMakie
using JPerpleX

@testset "Migmatites.jl" begin
    # Write your tests here.

    melt, host = equilibrate_open_system("23SD20A_melt-test1/MeltSource","23SD20A_melt-test1/Host",875,10000,800,9000)
    @show melt
    @show host
end
