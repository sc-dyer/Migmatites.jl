using Migmatites
using Test
using CairoMakie
using JPerpleX


@testset "Migmatites.jl" begin
    # Write your tests here.
    source_compo = init_meemum("23SD20A_melt-test1/MeltSource")
    h2ostart = 0.0001
    h2oend = 20.0

    
    source_compo1 = change_list_component(source_compo,h2ostart,"H2O")
    source_compo2 = change_list_component(source_compo,h2oend,"H2O")

    sources, melts, hosts = equilibrate_open_system("23SD20A_melt-test1/MeltSource","23SD20A_melt-test1/Host",875,10000,800,9000,source_compo1,source_compo2)

  

    host_melt_percent = Float64[]
    # for compo in source_compo_range
    #     source,melt, host = equilibrate_open_system.("23SD20A_melt-test1/MeltSource","23SD20A_melt-test1/Host",875,10000,800,9000,source_compo = compo,suppresswarn=true)
    #     host_melt = Phase(name="melt",composition = Component[])
    #     try
    #         host_melt = get_melt(host)
    #     catch
    #         push!(host_melt_percent,0.0)
    #     else
    #         push!(host_melt_percent,host_melt.vol/host.vol*100)
    #     end
    # end
    for host in hosts 
        host_melt = get_melt(host)
        push!(host_melt_percent,host_melt.vol/host.vol*100)
    end
    @show range(h2ostart,h2oend,length=100)
    @show host_melt_percent
end
