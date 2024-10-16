using Migmatites
using Test
using CairoMakie
using JPerpleX


@testset "Migmatites.jl" begin
    # Write your tests here.

    source,melt, host = equilibrate_open_system("23SD20A_melt-test1/MeltSource","23SD20A_melt-test1/Host",850,10000,800,9000)
    @show source.composition
    @show name.(source.phases)
    @show melt.composition
    @show name.(melt.phases)
    @show host.composition
    @show name.(host.phases)
end
