module Migmatites

export
    getmelt,
    equilibrate_open_system,
    equilibrate_closed_system,
    balance_component

using
    DocStringExtensions,
    JPerpleX,
    Optim

#Function that takes a PetroSystem and calculates u_H2O for given P-T-X_H2O conditions

#(Open system calculation) Function that calculates X_H2O required for given u_H2O -> Investigate if this is possible with Perple_X

#Function that iteratively calculates X_H2O at equilibrium for melt and host (PetroSystems) 
#of given proportions (options for volume, mol, mass?) at given P-T

#Function that calculates and extracts melt phase of parent rock at given P-T conditions
#Will need to run init_meemum before running this
"""
$(TYPEDSIGNATURES)

Calculates the phase equilibria at specified T and P conditions for given composition and 
returns the melt phase in the calculated system. Must run init_meemum before running this function.
"""
function getmelt(composition,T,P; suppresswarn = false)
    system = minimizepoint(composition, T, P, suppresswarn = suppresswarn)
    return getmelt(system)
end

"""
$(TYPEDSIGNATURES)

Returns the melt phase in a 'PetroSystem' if one exists, returns an empty phase with 0 mols of every component 
if no melt exists
"""
function getmelt(system)
    return getphase(system,r"melt")[1]
end

"""
$(TYPEDSIGNATURES)

This is a two step function which calculates the system and phase properties of the melt source rock at 'source_T' and 
'source_P', then use the chemical potential of `equilib_component` to calculate the system and phase properties of 
the host rock at `host_T` and `host_P`. If `source_compo` or `host_compo` are not provided, it will default to 
the compositions defined in the given dat files: `sourcefile` and `hostfile`. For this function to work as expected 
your sourcefile must have only P and T as independent variables, and your hostfile must have P, T and μ of 
one component as independent variables.
"""
function equilibrate_open_system(sourcefile,hostfile,source_T,source_P, host_T,host_P;source_compo = nothing, host_compo = nothing,equilib_component = "H2O",suppresswarn = false, phasefunc = [])
    #Need to run these in sequence,first calc melt system at T and P
    #Then get uH2O from that melt and calc host system for given uH2O
    #Complicated by how datafiles are defined
    sourcelib = init_meemum(sourcefile)
    scomp = getcompo(sourcelib)
    if !isnothing(source_compo)
        scomp = source_compo
    end
    source_system = minimizepoint(sourcelib,source_T,source_P,suppresswarn=suppresswarn, composition = scomp, phasefunc = phasefunc)
   
    melt = getmelt(source_system)
    if mol(melt) ≈ 0
        println("No melt generated at ",source_T," °C and ", source_P, " bar" )
        return
    end

    # close_meemum!(sourcelib)
    #Have to calculate melt conditions at host T and P
    # meltlib = init_meemum(meltfile)
    melt_system = minimizepoint(sourcelib,host_T,host_P,suppresswarn=suppresswarn,composition = melt.composition.*100, phasefunc = phasefunc)

    μ = melt_system.composition[findchemical(melt_system.composition,equilib_component)].μ
    
    close_meemum!(sourcelib)

    hostlib = init_meemum(hostfile)
    hcomp = getcompo(hostlib)
    if !isnothing(host_compo)
        hcomp = host_compo
    end
    
    host_system = minimizepoint(hostlib,host_T,host_P,μ1 = μ,suppresswarn=suppresswarn,composition=hcomp, phasefunc = phasefunc)

    close_meemum!(hostlib)
    return source_system,melt_system, host_system

end
"""
$(TYPEDSIGNATURES)

This does the same as the "single point" version of this function but provides a range of output for 
variable source compositions on a linear range between 'source_compo_start' and 'source_compo_end'.
"""
function equilibrate_open_system(sourcefile, hostfile,source_T,source_P, host_T,host_P, source_compo_start, source_compo_end; host_compo = nothing,equilib_component = "H2O",suppresswarn = false, steps = 100, phasefunc = [])
    #Should define another function for multiple compositions at source
    #So that init_meemum does not need to be run a billion times
    
    #Need to run these in sequence,first calc melt system at T and P
    #Then get uH2O from that melt and calc host system for given uH2O
    #Complicated by how datafiles are defined
    sourcelib = init_meemum(sourcefile)   

    source_compo_range = range(source_compo_start,source_compo_end,steps)
    source_systems = PetroSystem[]
   
    for compo in source_compo_range
        
        ssys= minimizepoint(sourcelib,source_T,source_P,suppresswarn=suppresswarn, composition = compo, phasefunc = phasefunc)
        push!(source_systems, ssys)
    end

  
    melt_systems = PetroSystem[]
    for ssys in source_systems
        melt = getmelt(ssys)
        if mol(melt) ≈ 0
            println("No melt generated at ",source_T," °C and ", source_P, " bar with composition: ", ssys.composition)
            #Done to maintain array shapes
            push!(melt_systems,PetroSystem())
        else
            msys = minimizepoint(sourcelib,host_T,host_P,suppresswarn=suppresswarn,composition = melt.composition.*100, phasefunc = phasefunc)
            push!(melt_systems,msys)
        end
    end
    
    close_meemum!(sourcelib)


    #Have to calculate melt conditions at host T and P
    hostlib = init_meemum(hostfile)
    hcomp = getcompo(hostlib)
    if !isnothing(host_compo)
        hcomp = host_compo
    end
    host_systems = PetroSystem[]
    for msys in melt_systems 
       
        if length(msys.composition) > 0
            index = findchemical(msys.composition,equilib_component)
            if index > 0
                μ = msys.composition[index].μ
                hsys = minimizepoint(hostlib,host_T,host_P,μ1 = μ,suppresswarn=suppresswarn,composition=hcomp, phasefunc = phasefunc)
                push!(host_systems,hsys)
            else
                push!(host_systems,PetroSystem())
            end
        else
            push!(host_systems,PetroSystem())
        end

    end
    close_meemum!(hostlib)
    return source_systems,melt_systems, host_systems

end

"""
$(TYPEDSIGNATURES)

This does the same as the "single point" version of this function but provides a range of output for 
variable source compositions on a linear range between 'source_T1' and 'source_T2'.
"""
function equilibrate_open_system(sourcefile, hostfile,source_T1, source_T2, source_P, host_T,host_P;source_compo = nothing, host_compo = nothing,equilib_component = "H2O",suppresswarn = false, steps = 100, phasefunc = [])
    #Should define another function for multiple compositions at source
    #So that init_meemum does not need to be run a billion times
    
    #Need to run these in sequence,first calc melt system at T and P
    #Then get uH2O from that melt and calc host system for given uH2O
    #Complicated by how datafiles are defined

    sourcelib = init_meemum(sourcefile)   
    scomp = getcompo(sourcelib)
    if !isnothing(source_compo)
        scomp = source_compo
    end

    Trange = range(source_T1, source_T2,steps)
    source_systems = PetroSystem[]
   
    for T in Trange  
        ssys= minimizepoint(sourcelib,T,source_P,suppresswarn=suppresswarn, composition = scomp, phasefunc = phasefunc)
        push!(source_systems, ssys)
    end

    
    melt_systems = PetroSystem[]
    for ssys in source_systems
        melt = getmelt(ssys)
        if mol(melt) ≈ 0
            println("No melt generated at ",source_T," °C and ", source_P, " bar with composition: ", ssys.composition)
            #Done to maintain array shapes
            push!(melt_systems,PetroSystem())
        else
            msys = minimizepoint(sourcelib,host_T,host_P,suppresswarn=suppresswarn,composition = melt.composition.*100, phasefunc = phasefunc)
            push!(melt_systems,msys)
        end
    end
    
    close_meemum!(sourcelib)


    #Have to calculate melt conditions at host T and P
    hostlib = init_meemum(hostfile)
    hcomp = getcompo(hostlib)
    if !isnothing(host_compo)
        hcomp = host_compo
    end
    host_systems = PetroSystem[]
    for msys in melt_systems 
       
        if length(msys.composition) > 0
            index = findchemical(msys.composition,equilib_component)
            if index > 0
                μ = msys.composition[index].μ
                hsys = minimizepoint(hostlib,host_T,host_P,μ1 = μ,suppresswarn=suppresswarn,composition=hcomp, phasefunc = phasefunc)
                push!(host_systems,hsys)
            else
                push!(host_systems,PetroSystem())
            end
        else
            push!(host_systems,PetroSystem())
        end

    end
    close_meemum!(hostlib)
    return source_systems,melt_systems, host_systems

end

"""
$(TYPEDSIGNATURES)

This is a two step function which calculates the system and phase properties of the melt source rock at 'source_T' and
'source_P', then use the chemical potential of `equilib_component` and `melt_proportion` to calculate the system and 
phase properties of the host rock at `host_T` and `host_P` attempts to converge on a proportion of the component of 
interest where both melt and host have the same chemical potential in that component. `melt_proportion` is the molar ratio 
of melt to host rock. If `source_compo` or `host_compo` are not provided, it will default to the compositions defined 
in the given dat files: `sourcefile` and `hostfile`. For this function to work as expected your sourcefile and 
hostfile must have only P and T as independent variables and they must contain the same components.
"""
function equilibrate_closed_system(meemumlib,sourcecompo,hostcompo,source_T,source_P, host_T,host_P, melt_proportion;equilib_component = "H2O",suppresswarn = false, phasefunc = [],μthreshold = 2.0,xthreshold = 1E-5)
    #Set up: get the host rock composition
   
    source_system = minimizepoint(meemumlib,source_T,source_P,suppresswarn=suppresswarn, composition = sourcecompo, phasefunc = phasefunc)
   
    melt = getmelt(source_system)
    if mol(melt) ≈ 0
        println("No melt generated at ",source_T," °C and ", source_P, " bar" )
        return
    end

    # melt_system = minimizepoint(meemumlib,host_T,host_P,composition = melt.composition.*100, suppresswarn = suppresswarn, phasefunc = phasefunc)
    # host_system = minimizepoint(meemumlib,host_T,host_P,composition = hostcompo, suppresswarn = suppresswarn, phasefunc = phasefunc)

    # μ_melt = melt_system.composition[findchemical(melt_system.composition,equilib_component)].μ
    # μ_host = host_system.composition[findchemical(host_system.composition,equilib_component)].μ

    # if μ_melt < μ_host
    #     println("Host μ is higher than melt μ")
    #     return source_system,melt_system, host_system
    # end
    
    
    # count = 1
    
    # x_change = 0.001
    # if molfrac(melt_system.composition,equilib_component) > 0.15
    #     x_change = 0.01
    # end
    
    # # melt_system, host_system = recursive_equilibrate(sourcelib,melt.composition,hcomp,melt_proportion,host_T, host_P,equilib_component=equilib_component,suppresswarn=suppresswarn,phasefunc=phasefunc)
    # while abs(μ_melt - μ_host) >μthreshold && molfrac(melt_system.composition,equilib_component) > xthreshold && abs(x_change) > 1e-9 && count < 100
    #     println("Loop ", count)
    #     # x_change = molfrac(melt_system.composition,equilib_component)/divisor
    #     @show x_change
    #     newmeltcompo, newhostcompo = balance_component(melt_system.composition,host_system.composition,equilib_component,x_change,melt_proportion)
    #     chemi = findchemical(newmeltcompo,equilib_component)
    #     chemj = findchemical(newhostcompo,equilib_component)
    #     if chemi > 0 && chemj > 0
    #         newmelt_system = minimizepoint(meemumlib,host_T,host_P,composition = newmeltcompo, suppresswarn = suppresswarn, phasefunc = phasefunc)
    #         newhost_system = minimizepoint(meemumlib,host_T,host_P,composition = newhostcompo, suppresswarn = suppresswarn, phasefunc = phasefunc)
    #         nμ_melt = getchemical(newmelt_system.composition,equilib_component).μ
    #         nμ_host = getchemical(newhost_system.composition,equilib_component).μ
            
    #         if abs(nμ_melt-nμ_host) > abs(μ_melt-μ_host)
    #             x_change /= 10
    #         else
    #             if (nμ_melt - nμ_host)/x_change < 0
    #                 x_change /= -10
    #             end
    #             melt_system = newmelt_system
    #             host_system = newhost_system
    #             μ_melt = nμ_melt
    #             μ_host = nμ_host
            
    #             Xh2o_melt = molfrac(melt_system.composition,equilib_component)
    #             Xh2o_host = molfrac(host_system.composition,equilib_component)
                
    #             @show μ_host, μ_melt
    #             @show Xh2o_host, Xh2o_melt
    #         end
            
    #     else
    #         x_change /= 10
    #     end
    #     count +=1
    # end

    
    meltcompo = melt.composition.*100
    res = optimize(x -> closed_system_optim(meemumlib,melt.composition,hostcompo,host_T,host_P,melt_proportion,x[1],equilib_component=equilib_component,suppresswarn=suppresswarn,phasefunc=phasefunc),[0.0])
    
    x_change = Optim.minimizer(res)[1]
   
    newmeltcompo, newhostcompo = balance_component(meltcompo,hostcompo,equilib_component,x_change,melt_proportion)
    

    melt_system = minimizepoint(meemumlib,host_T,host_P,composition = newmeltcompo, suppresswarn = suppresswarn, phasefunc = phasefunc)
    host_system = minimizepoint(meemumlib,host_T,host_P,composition = newhostcompo, suppresswarn = suppresswarn, phasefunc = phasefunc)

    return source_system,melt_system, host_system

end

function closed_system_optim(meemumlib,meltcompo,hostcompo,host_T,host_P,melt_proportion,x_change;equilib_component = "H2O",suppresswarn=false,phasefunc = [])
    

    newmeltcompo, newhostcompo = balance_component(meltcompo,hostcompo,equilib_component,x_change,melt_proportion)

    melt_system = minimizepoint(meemumlib,host_T,host_P,composition = newmeltcompo, suppresswarn = suppresswarn, phasefunc = phasefunc)
    host_system = minimizepoint(meemumlib,host_T,host_P,composition = newhostcompo, suppresswarn = suppresswarn, phasefunc = phasefunc)

    chemi = findchemical(melt_system.composition,equilib_component)
    chemj = findchemical(host_system.composition,equilib_component)
    if chemi == 0 || chemj == 0
        return Inf
    end
    if melt_system.composition[chemi].mol < 0.0 || host_system.composition[chemj].mol < 0.0
        return Inf
    end
    μ_melt = melt_system.composition[chemi].μ
    μ_host = host_system.composition[chemj].μ

    return abs(μ_melt-μ_host)
end

function balance_component(compo1, compo2, component,diff, proportion)
    #Takes from system1 and adds to system2 according to proportion of system1 to system2
    #So if there is 2 times as much of system1 as system2, and system1 is losing 1% of its H2O 
    #to system2, system2 will gain twice as much of the H2O
    
    

    #Normalize to 1 mol
    x1 = molfrac(compo1,component)
    x2 = molfrac(compo2,component)

    x1 = x1 - diff
    x2 = x2 + diff*proportion
    #Reverse molfrac normalization
    mol1 = x1*sum_mols(compo1)
    mol2 = x2*sum_mols(compo2)
    compo1 = change_list_component(compo1,mol1,component)
    compo2 = change_list_component(compo2,mol2,component)

    return compo1, compo2
    #Get proportion of component of interestt
end


end