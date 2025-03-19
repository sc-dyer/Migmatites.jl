using Migmatites
using Test
using CairoMakie
using JPerpleX
using DataFrames
# using CSV

@testset "Migmatites.jl" begin
    #Test1
    # @testset "Open system tests" begin
    #     hostlib = init_meemum("23SD20A_melt-test1/Host")
    #     host = minimizepoint(hostlib,800,9000,μ1 = -316240)

    #     hosth2o = getchemical(host.composition,"H2O")

    #     @test round(hosth2o.mol,digits=3) ≈ 0.414
    #     @test round(hosth2o.μ,sigdigits = 6) ≈ -316240
    #     @test round(host.phases[1].vol/host.vol*100,digits=2) ≈ 5.15
    #     close_meemum!(hostlib)

    #     source, melt, host = equilibrate_open_system("23SD20A_melt-test1/MeltSource","23SD20A_melt-test1/Host",875,10000,800,9000)

    #     sourcemelt = getmelt(source)
    #     sourcemelth2o= getchemical(sourcemelt.composition,"H2O")
    #     @test round(sourcemelth2o.mol,digits=5) ≈ 0.40998
    #     @test round(sourcemelt.vol/source.vol * 100,digits=2) ≈ 6.89
    #     # @test sourcemelt.composition .≈ melt.composition

    #     melth2o = getchemical(melt.composition,"H2O")
    #     @test round(melth2o.μ,sigdigits=6) ≈ -316239
    #     @test round(melth2o.mol,digits=3) ≈ 40.998
        
    #     hosth2o = getchemical(host.composition,"H2O")

    #     @test round(hosth2o.mol,digits=3) ≈ 0.414
    #     @test round(hosth2o.μ,sigdigits = 6) ≈ -316239
    #     @test round(host.phases[1].vol/host.vol*100,digits=2) ≈ 5.15

    #     #Test2
    #     source, melt, host = equilibrate_open_system("23SD20A_melt-test2/MeltSource","23SD20A_melt-test2/Host",875,10000,800,9000)

    #     sourcemelt = getmelt(source)
    #     sourcemelth2o= getchemical(sourcemelt.composition,"H2O")
    #     @test round(sourcemelth2o.mol,digits=5) ≈ 0.67574
    #     @test round(sourcemelt.vol/source.vol * 100,digits=2) ≈ 97.78
    #     # @test sourcemelt.composition .≈ melt.composition

    #     melth2o = getchemical(melt.composition,"H2O")
    #     @test round(melth2o.μ,sigdigits=6) ≈ -311522
    #     @test round(melth2o.mol,digits=3) ≈ 67.574
        
    #     hosth2o = getchemical(host.composition,"H2O")

    #     @test round(hosth2o.mol,digits=3) ≈ 47.143
    #     @test round(hosth2o.μ,sigdigits=6) ≈ -311522
    #     hostmelt = getmelt(host)
    #     @test round(hostmelt.vol/host.vol*100,digits=2) ≈ 96.07
    #     hostmelth2o = getchemical(hostmelt.composition,"H2O")
    #     @test round(hostmelth2o.mol,digits=2) ≈ 0.67



    #     sourcelib = init_meemum("23SD20A_melt-test1/MeltSource")
    #     source_compo = getcompo(sourcelib)
    #     close_meemum!(sourcelib)

    #     h2ostart = 1.0
    #     h2oend = 50.0
        
    #     source_compo1 = change_list_component(source_compo,h2ostart,"H2O")
    #     source_compo2 = change_list_component(source_compo,h2oend,"H2O")
    
    #     sources, melts, hosts = equilibrate_open_system("23SD20A_melt-test1/MeltSource","23SD20A_melt-test1/Host",875,10000,800,9000,source_compo1,source_compo2, steps =10)

    #     source = sources[1]
    #     melt = melts[1]
    #     host = hosts[1]
    #     sourcemelt = getmelt(source)
    #     sourcemelth2o= getchemical(sourcemelt.composition,"H2O")
    #     @test round(sourcemelth2o.mol,digits=5) ≈ 0.40998
    #     @test round(sourcemelt.vol/source.vol * 100,digits=2) ≈ 6.89
    #     # @test sourcemelt.composition .≈ melt.composition

    #     melth2o = getchemical(melt.composition,"H2O")
    #     @test round(melth2o.μ,sigdigits=6) ≈ -316239
    #     @test round(melth2o.mol,digits=3) ≈ 40.998
        
    #     hosth2o = getchemical(host.composition,"H2O")

    #     @test round(hosth2o.mol,digits=3) ≈ 0.414
    #     @test round(hosth2o.μ,sigdigits = 6) ≈ -316239
    #     @test round(host.phases[1].vol/host.vol*100,digits=2) ≈ 5.15

    #     source = sources[10]
    #     melt = melts[10]
    #     host = hosts[10]
    #     sourcemelt = getmelt(source)
    #     sourcemelth2o= getchemical(sourcemelt.composition,"H2O")
    #     @test round(sourcemelth2o.mol,digits=5) ≈ 0.67574
    #     @test round(sourcemelt.vol/source.vol * 100,digits=2) ≈ 97.78
    #     # @test sourcemelt.composition .≈ melt.composition

    #     melth2o = getchemical(melt.composition,"H2O")
    #     @test round(melth2o.μ,sigdigits=6) ≈ -311522
    #     @test round(melth2o.mol,digits=3) ≈ 67.574
        
    #     hosth2o = getchemical(host.composition,"H2O")

    #     @test round(hosth2o.mol,digits=3) ≈ 47.143
    #     @test round(hosth2o.μ,sigdigits=6) ≈ -311522
    #     hostmelt = getmelt(host)
    #     @test round(hostmelt.vol/host.vol*100,digits=2) ≈ 96.07
    #     hostmelth2o = getchemical(hostmelt.composition,"H2O")
    #     @test round(hostmelth2o.mol,digits=2) ≈ 0.67

    #     fig = Figure(size = (600,450))
    #     ax = Axis(fig[1,1])
    #     modebox!(ax,range(h2ostart,h2oend,10),hosts)
    #     fig[1,2] = Legend(fig,ax)
    #     save("23SD20A_melt-test1/Host.svg",fig)

    # end

    # @testset "Closed system tests" begin

    #     hostlib = init_meemum("23SD20A_melt-test3/Host")
    #     hostcompo = hostlib.composition
    #     close_meemum!(hostlib)

    #     sourcelib = init_meemum("23SD20A_melt-test3/MeltSource")
    #     sourcecompo = sourcelib.composition
        

    #     compo1, compo2 = balance_component(sourcecompo,hostcompo,"H2O",0.01,2.0)

    #     @test round(molfrac(compo1,"H2O"),sigdigits = 3) ≈ 0.161 
    #     @test round(molfrac(compo2,"H2O"),sigdigits = 3) ≈ 0.0206

    #     source, melt, host = equilibrate_closed_system(sourcelib,sourcecompo,hostcompo,875,10000,800,9000,2.0)
        
    #     sourcemelt = getmelt(source)
    #     sourcemelth2o= getchemical(sourcemelt.composition,"H2O")
    #     @test round(sourcemelth2o.mol,digits=5) ≈ 0.46122
    #     @test round(sourcemelt.vol/source.vol * 100,digits=2) ≈ 96.65
    #     # @test sourcemelt.composition .≈ melt.composition
        
    #     melth2o = getchemical(melt.composition,"H2O")
    #     @test round(melth2o.μ,sigdigits=6) ≈ -316142
    #     @test round(melth2o.mol,digits=3) ≈ 45.111
        
    #     hosth2o = getchemical(host.composition,"H2O")

    #     @test round(hosth2o.mol,digits=3) ≈ 0.87
    #     @test round(hosth2o.μ,sigdigits = 6) ≈ -316143
    #     hostmelt = getmelt(host)
    #     @test round(hostmelt.vol/host.vol*100,digits=2) ≈ 2.24

    #     h2ostart = 1.0
    #     h2oend = 50.0
        
    #     source_compo1 = change_list_component(sourcecompo,h2ostart,"H2O")
    #     source_compo2 = change_list_component(sourcecompo,h2oend,"H2O")

    #     sourcecomporange = range(source_compo1,source_compo2,10)
    #     # sources = PetroSystem[]
    #     # melts = PetroSystem[]
    #     # hosts = PetroSystem[]
    #     # for i in 1:lastindex(sourcecomporange)
    #     #     sourcei, melti, hosti = equilibrate_closed_system(sourcelib,sourcecomporange[i],hostcompo,875,10000,800.9000,2.0)
    #     # end
    #     systems = equilibrate_closed_system.((sourcelib,),sourcecomporange,(hostcompo,),875,10000,800,9000,2.0)
    #     # for i in 1:lastindex(systems)
    #     #     @show molfrac(getcompo(systems[i][2]),"H2O")
    #     #     @show molfrac(getcompo(systems[i][3]),"H2O")
    #     # end
    #     source = systems[10][1]
    #     melt = systems[10][2]
    #     host = systems[10][3]
    #     sourcemelt = getmelt(source)
    #     sourcemelth2o= getchemical(sourcemelt.composition,"H2O")
    #     @test round(sourcemelth2o.mol,digits=5) ≈ 0.67574
    #     @test round(sourcemelt.vol/source.vol * 100,digits=2) ≈ 97.78
    #     # @test sourcemelt.composition .≈ melt.composition

    #     melth2o = getchemical(melt.composition,"H2O")
    #     @test round(melth2o.μ,sigdigits=6) ≈ -314457
    #     @test round(melth2o.mol,digits=3) ≈ 39.122
        
    #     hosth2o = getchemical(host.composition,"H2O")

    #     @test round(hosth2o.mol,digits=3) ≈ 28.919
    #     @test round(hosth2o.μ,sigdigits=6) ≈ -314458
    #     hostmelt = getmelt(host)
    #     @test round(hostmelt.vol/host.vol*100,digits=2) ≈ 95.26
    #     hostmelth2o = getchemical(hostmelt.composition,"H2O")
    #     @test round(hostmelth2o.mol,digits=2) ≈ 0.55

    #     hosts = [sys[3] for sys in systems]
    #     fig = Figure(size = (600,450))
    #     ax = Axis(fig[1,1])
    #     modebox!(ax,range(h2ostart,h2oend,10),hosts)
    #     fig[1,2] = Legend(fig,ax)
    #     save("23SD20A_melt-test3/Host.svg",fig)

    #     close_meemum!(sourcelib)
    # end

    # @testset "Melting path tests" begin
    #     sourcelib = init_meemum("23SD20A_melt-test1/MeltSource")
    #     source_compo = getcompo(sourcelib)
    #     t_path = [600,1100]
    #     p_path = [10000,10000]

    #     melts, restite_compos, system_steps, extract_T, extract_P = opensystem_melting_path(sourcelib, source_compo,t_path,p_path)
    #     h2ostart = 0.001
        
        
    #     # source_compo = change_list_component(source_compo,h2ostart,"H2O")
    #     # t_path2 = [800,800]
    #     # x_range = [h2ostart,0.2]
        
    #     # melts, restite_compos, system_steps, extract_T, extract_X = opensystem_melting_path(sourcelib, source_compo,t_path2,x_range,10000,numsteps=1000)
    #     # @show melts
    #     # @show extract_X

    #     close_meemum!(sourcelib)
        
       
    #     sourcefile = "23SD20A_melt-test4/MeltSource"
    #     hostfile = "23SD20A_melt-test4/Host"
    #     melts_f, hosts, restite_compos, system_steps, extract_T, extract_P = melt_cycle_model_open(hostfile,sourcefile,800,9000,t_path,p_path,numsteps = 1000)

        
    # end

    @testset "Partitioning tests" begin
        
        function feldspar_conditions(phase)
            if lowercase(phase.name) == "fsp"
                k2o = getchemical(phase.composition,"K2O")
                na2o = getchemical(phase.composition,"Na2O")
                cao = getchemical(phase.composition,"CaO")
            
                if mol(k2o*2) >0.1
                    return changename(phase,"Afs")
                else
                    return changename(phase,"Pl")
                end
            else
                return phase
            end
        end

        testdf, xvar, yvar = read_werami_output("20SD06_R2c_1_trimmed.phm",iscelsius=true,iskbar = true)

        petrodf = werami_to_petrosys(testdf, xvar, yvar,phasefunc =[feldspar_conditions])
        
        kd_vals = import_Kd_values("zr_kd.csv")
        add_te!.(petrodf[!,:system],(TraceElement("Zr",91.224,1001),))
        calc_Zr!.(petrodf[!,:system],petrodf[!,xvar],(kd_vals,))
        petrodf[!,"Zircon (mass frac)"] = calc_zircon!.(petrodf[!,:system])
        petrorows = eachrow(petrodf)

        @show petrorows[126]
        
        testitem = petrorows[126][:system]
        melt = getmelt(testitem)
        meltcompo = melt.composition
        na = mol(getchemical(meltcompo,"Na2O"))*2
        k = mol(getchemical(meltcompo,"K2O"))*2
        ca = mol(getchemical(meltcompo,"CaO"))
        al = mol(getchemical(meltcompo,"Al2O3"))*2
        si = mol(getchemical(meltcompo,"SiO2"))
        
        @show massfrac(melt,testitem)
        @show zircon_saturation(getmelt(testitem),799)
        zirconweight = petrorows[126]["Zircon (mass frac)"]
        @test round(concentration(getchemical(getphase(testitem,"melt")[1],"Zr")),sigdigits=3) ≈ 205
        @test round(concentration(getchemical(getphase(testitem,"Cpx")[1],"Zr")),sigdigits=3) ≈ 25.6
        @test round(concentration(getchemical(getphase(testitem,"Ilm")[1],"Zr")),sigdigits=3) ≈ 472
        @test round(concentration(getchemical(getphase(testitem,"Sp")[1],"Zr")),sigdigits=3) ≈ 24.6
        @test round(zirconweight,sigdigits=3) ≈ 0.00198
    end
    
end
