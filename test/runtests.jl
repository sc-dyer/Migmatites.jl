using Migmatites
using Test
using CairoMakie
using JPerpleX


@testset "Migmatites.jl" begin
    #Test1
    @testset "Open system tests" begin
        hostlib = init_meemum("23SD20A_melt-test1/Host")
        host = minimizepoint(hostlib,800,9000,μ1 = -316240)

        hosth2o = getchemical(host.composition,"H2O")

        @test round(hosth2o.mol,digits=3) ≈ 0.414
        @test round(hosth2o.μ,sigdigits = 6) ≈ -316240
        @test round(host.phases[1].vol/host.vol*100,digits=2) ≈ 5.15
        close_meemum!(hostlib)

        source, melt, host = equilibrate_open_system("23SD20A_melt-test1/MeltSource","23SD20A_melt-test1/Host",875,10000,800,9000)

        sourcemelt = getmelt(source)
        sourcemelth2o= getchemical(sourcemelt.composition,"H2O")
        @test round(sourcemelth2o.mol,digits=5) ≈ 0.40998
        @test round(sourcemelt.vol/source.vol * 100,digits=2) ≈ 6.89
        # @test sourcemelt.composition .≈ melt.composition

        melth2o = getchemical(melt.composition,"H2O")
        @test round(melth2o.μ,sigdigits=6) ≈ -316239
        @test round(melth2o.mol,digits=3) ≈ 40.998
        
        hosth2o = getchemical(host.composition,"H2O")

        @test round(hosth2o.mol,digits=3) ≈ 0.414
        @test round(hosth2o.μ,sigdigits = 6) ≈ -316239
        @test round(host.phases[1].vol/host.vol*100,digits=2) ≈ 5.15

        #Test2
        source, melt, host = equilibrate_open_system("23SD20A_melt-test2/MeltSource","23SD20A_melt-test2/Host",875,10000,800,9000)

        sourcemelt = getmelt(source)
        sourcemelth2o= getchemical(sourcemelt.composition,"H2O")
        @test round(sourcemelth2o.mol,digits=5) ≈ 0.67574
        @test round(sourcemelt.vol/source.vol * 100,digits=2) ≈ 97.78
        # @test sourcemelt.composition .≈ melt.composition

        melth2o = getchemical(melt.composition,"H2O")
        @test round(melth2o.μ,sigdigits=6) ≈ -311522
        @test round(melth2o.mol,digits=3) ≈ 67.574
        
        hosth2o = getchemical(host.composition,"H2O")

        @test round(hosth2o.mol,digits=3) ≈ 47.143
        @test round(hosth2o.μ,sigdigits=6) ≈ -311522
        hostmelt = getmelt(host)
        @test round(hostmelt.vol/host.vol*100,digits=2) ≈ 96.07
        hostmelth2o = getchemical(hostmelt.composition,"H2O")
        @test round(hostmelth2o.mol,digits=2) ≈ 0.67



        sourcelib = init_meemum("23SD20A_melt-test1/MeltSource")
        source_compo = getcompo(sourcelib)
        close_meemum!(sourcelib)

        h2ostart = 1.0
        h2oend = 50.0
        
        source_compo1 = change_list_component(source_compo,h2ostart,"H2O")
        source_compo2 = change_list_component(source_compo,h2oend,"H2O")
    
        sources, melts, hosts = equilibrate_open_system("23SD20A_melt-test1/MeltSource","23SD20A_melt-test1/Host",875,10000,800,9000,source_compo1,source_compo2, steps =10)

        source = sources[1]
        melt = melts[1]
        host = hosts[1]
        sourcemelt = getmelt(source)
        sourcemelth2o= getchemical(sourcemelt.composition,"H2O")
        @test round(sourcemelth2o.mol,digits=5) ≈ 0.40998
        @test round(sourcemelt.vol/source.vol * 100,digits=2) ≈ 6.89
        # @test sourcemelt.composition .≈ melt.composition

        melth2o = getchemical(melt.composition,"H2O")
        @test round(melth2o.μ,sigdigits=6) ≈ -316239
        @test round(melth2o.mol,digits=3) ≈ 40.998
        
        hosth2o = getchemical(host.composition,"H2O")

        @test round(hosth2o.mol,digits=3) ≈ 0.414
        @test round(hosth2o.μ,sigdigits = 6) ≈ -316239
        @test round(host.phases[1].vol/host.vol*100,digits=2) ≈ 5.15

        source = sources[10]
        melt = melts[10]
        host = hosts[10]
        sourcemelt = getmelt(source)
        sourcemelth2o= getchemical(sourcemelt.composition,"H2O")
        @test round(sourcemelth2o.mol,digits=5) ≈ 0.67574
        @test round(sourcemelt.vol/source.vol * 100,digits=2) ≈ 97.78
        # @test sourcemelt.composition .≈ melt.composition

        melth2o = getchemical(melt.composition,"H2O")
        @test round(melth2o.μ,sigdigits=6) ≈ -311522
        @test round(melth2o.mol,digits=3) ≈ 67.574
        
        hosth2o = getchemical(host.composition,"H2O")

        @test round(hosth2o.mol,digits=3) ≈ 47.143
        @test round(hosth2o.μ,sigdigits=6) ≈ -311522
        hostmelt = getmelt(host)
        @test round(hostmelt.vol/host.vol*100,digits=2) ≈ 96.07
        hostmelth2o = getchemical(hostmelt.composition,"H2O")
        @test round(hostmelth2o.mol,digits=2) ≈ 0.67

        fig = Figure(size = (600,450))
        ax = Axis(fig[1,1])
        phasemode!(ax,range(h2ostart,h2oend,10),hosts)
        fig[1,2] = Legend(fig,ax)
        save("23SD20A_melt-test1/Host.svg",fig)

    end

    @testset "Closed system tests" begin

        hostlib = init_meemum("23SD20A_melt-test3/Host")
        hostcompo = hostlib.composition
        close_meemum!(hostlib)

        sourcelib = init_meemum("23SD20A_melt-test3/MeltSource")
        sourcecompo = sourcelib.composition
        

        compo1, compo2 = balance_component(sourcecompo,hostcompo,"H2O",0.01,2.0)

        @test round(molfrac(compo1,"H2O"),sigdigits = 3) ≈ 0.0833
        @test round(molfrac(compo2,"H2O"),sigdigits = 3) ≈ 0.0256

        source, melt, host = equilibrate_closed_system(sourcelib,sourcecompo,hostcompo,875,10000,800,9000,2.0)
        
        sourcemelt = getmelt(source)
        sourcemelth2o= getchemical(sourcemelt.composition,"H2O")
        @test round(sourcemelth2o.mol,digits=5) ≈ 0.40998
        @test round(sourcemelt.vol/source.vol * 100,digits=2) ≈ 6.89
        # @test sourcemelt.composition .≈ melt.composition

        melth2o = getchemical(melt.composition,"H2O")
        @test round(melth2o.μ,sigdigits=6) ≈ -316239
        @test round(melth2o.mol,digits=3) ≈ 40.998
        
        hosth2o = getchemical(host.composition,"H2O")

        @test round(hosth2o.mol,digits=3) ≈ 0.414
        @test round(hosth2o.μ,sigdigits = 6) ≈ -316239
        @test round(host.phases[1].vol/host.vol*100,digits=2) ≈ 5.15

        h2ostart = 1.0
        h2oend = 50.0
        
        source_compo1 = change_list_component(sourcecompo,h2ostart,"H2O")
        source_compo2 = change_list_component(sourcecompo,h2oend,"H2O")

        sourcecomporange = range(source_compo1,source_compo2,10)
        @show sourcecomporange
        systems = equilibrate_closed_system.(fill(sourcelib,10),sourcecomporange,fill(hostcompo,10),875,10000,800,9000,2.0)
       
        source = systems[10][1]
        melt = systems[10][2]
        host = systems[10][3]
        sourcemelt = getmelt(source)
        sourcemelth2o= getchemical(sourcemelt.composition,"H2O")
        @test round(sourcemelth2o.mol,digits=5) ≈ 0.67574
        @test round(sourcemelt.vol/source.vol * 100,digits=2) ≈ 97.78
        # @test sourcemelt.composition .≈ melt.composition

        melth2o = getchemical(melt.composition,"H2O")
        @test round(melth2o.μ,sigdigits=6) ≈ -311522
        @test round(melth2o.mol,digits=3) ≈ 67.574
        
        hosth2o = getchemical(host.composition,"H2O")

        @test round(hosth2o.mol,digits=3) ≈ 47.143
        @test round(hosth2o.μ,sigdigits=6) ≈ -311522
        hostmelt = getmelt(host)
        @test round(hostmelt.vol/host.vol*100,digits=2) ≈ 96.07
        hostmelth2o = getchemical(hostmelt.composition,"H2O")
        @test round(hostmelth2o.mol,digits=2) ≈ 0.67

        hosts = [sys[3] for sys in systems]
        fig = Figure(size = (600,450))
        ax = Axis(fig[1,1])
        phasemode!(ax,range(h2ostart,h2oend,10),hosts)
        fig[1,2] = Legend(fig,ax)
        save("23SD20A_melt-test3/Host.svg",fig)

        close_meemum!(sourcelib)
    end

    
end
