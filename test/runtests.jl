using Migmatites
using Test
using CairoMakie
using JPerpleX


@testset "Migmatites.jl" begin
    #Test1
    hostlib = init_meemum("23SD20A_melt-test1/Host")
    host = minimizepoint(hostlib,800,9000,μ1 = -316240)

    hosth2o = getchemical(host.composition,"H2O")

    @test round(hosth2o.mol,digits=3) ≈ 0.414
    @test round(hosth2o.μ,sigdigits = 6) ≈ -316240
    @test round(host.phases[1].vol/host.vol*100,digits=2) ≈ 5.15
    close_meemum!(hostlib)

    source, melt, host = equilibrate_open_system("23SD20A_melt-test1/MeltSource","23SD20A_melt-test1/Melt","23SD20A_melt-test1/Host",875,10000,800,9000)

    sourcemelt = get_melt(source)
    sourcemelth2o= getchemical(sourcemelt.composition,"H2O")
    @test round(sourcemelth2o.mol,digits=5) ≈ 0.40998
    @test round(sourcemelt.vol/source.vol * 100,digits=2) ≈ 6.89
    # @test sourcemelt.composition .≈ melt.composition

    melth2o = getchemical(melt.composition,"H2O")
    @test round(melth2o.μ,sigdigits=6) ≈ -316241
    @test round(melth2o.mol,digits=5) ≈ 0.40998
    
    hosth2o = getchemical(host.composition,"H2O")

    @test round(hosth2o.mol,digits=3) ≈ 0.414
    @test round(hosth2o.μ,sigdigits = 6) ≈ -316241
    @test round(host.phases[1].vol/host.vol*100,digits=2) ≈ 5.15

    #Test2
    source, melt, host = equilibrate_open_system("23SD20A_melt-test2/MeltSource","23SD20A_melt-test1/Melt","23SD20A_melt-test2/Host",875,10000,800,9000)

    sourcemelt = get_melt(source)
    sourcemelth2o= getchemical(sourcemelt.composition,"H2O")
    @test round(sourcemelth2o.mol,digits=5) ≈ 0.67574
    @test round(sourcemelt.vol/source.vol * 100,digits=2) ≈ 97.78
    # @test sourcemelt.composition .≈ melt.composition

    melth2o = getchemical(melt.composition,"H2O")
    @test round(melth2o.μ,sigdigits=6) ≈ -311522
    @test round(melth2o.mol,digits=5) ≈ 0.67574
    
    hosth2o = getchemical(host.composition,"H2O")

    @test round(hosth2o.mol,digits=3) ≈ 47.143
    @test round(hosth2o.μ,sigdigits=6) ≈ -311522
    hostmelt = get_melt(host)
    @test round(hostmelt.vol/host.vol*100,digits=2) ≈ 96.07
    hostmelth2o = getchemical(hostmelt.composition,"H2O")
    @test round(hostmelth2o.mol,digits=2) ≈ 0.67



    # source_compo = init_meemum("23SD20A_melt-test1/MeltSource")
    # h2ostart = 0.0
    # h2oend = 1.0

    # h2oindex = findchemical[host.composition,"H2O"]
    
    # source_compo1 = change_list_component(source_compo,h2ostart,"H2O")
    # source_compo2 = change_list_component(source_compo,h2oend,"H2O")
    # @show source_compo1
    # sources, melts, hosts = equilibrate_open_system("23SD20A_melt-test1/MeltSource","23SD20A_melt-test1/Host",875,10000,800,9000,source_compo1,source_compo2, steps =10)

   
    # # T_start = 850
    # # T_end = 950
    # # sources, melts, hosts = equilibrate_open_system("23SD20A_melt-test1/MeltSource","23SD20A_melt-test1/Host",T_start,T_end,10000,800,9000)
    
    # host_melt_percent = Float64[]
    # hostH2O = Float64[]
    # for host in hosts 
    #     host_melt = get_melt(host)
    #     push!(host_melt_percent,host_melt.vol/host.vol*100)
    #     h2oindex = findchemical(host.composition,"H2O")
    #     if h2oindex > 0
    #         push!(hostH2O,concentration(host.composition[h2oindex]))
    #     else
    #         push!(hostH2O,0)
    #     end
    # end

    # @show range(h2ostart,h2oend,length=100)
    # # @show range(T_start,T_end,length=100)
    # @show host_melt_percent
    # @show hostH2O

    # meltH2O = Float64[]
    # meltuh2o = Float64[]
    
    # for melt in melts
    #     h2oindex = findchemical(melt.composition,"H2O")
    #     if h2oindex > 0
    #         push!(meltH2O,concentration(melt.composition[h2oindex]))
    #         push!(meltuh2o,melt.composition[h2oindex].μ)
    #     else
    #         push!(meltH2O,0)
    #         push!(meltuh2o,0)
    #     end
    # end

    # @show meltH2O
    # @show meltuh2o

    # source_melt_percent = Float64[]
    # sourceH2O = Float64[]
    # for source in sources
    #     source_melt = get_melt(source)
    #     push!(source_melt_percent, source_melt.vol/source.vol*100)
    #     h2oindex = findchemical(source.composition,"H2O")
    #     if h2oindex > 0
    #         push!(sourceH2O,concentration(source.composition[h2oindex]))
    #     else
    #         push!(sourceH2O,0)
    #     end
    # end

   
    # @show source_melt_percent
    # @show sourceH2O
end
