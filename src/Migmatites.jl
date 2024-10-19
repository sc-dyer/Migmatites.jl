module Migmatites

export
    get_melt,
    equilibrate_open_system

using
    DocStringExtensions,
    JPerpleX

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
function get_melt(composition,T,P; suppresswarn = false)
    system = minimizepoint(composition, T, P, suppresswarn = suppresswarn)
    return get_melt(system)
end

"""
$(TYPEDSIGNATURES)

Returns the melt phase in a 'PetroSystem' if one exists, returns an empty phase with 0 mols of every component
if no melt exists
"""
function get_melt(system)
    for phase in system.phases
        if contains(lowercase(phase.name),"melt")
            return phase
        end
    end
    return Phase(name="Melt",composition = zero.(system.composition))
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
function equilibrate_open_system(sourcefile,meltfile,hostfile,source_T,source_P, host_T,host_P;source_compo = nothing, host_compo = nothing,equilib_component = "H2O",suppresswarn = false)
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
    source_system = minimizepoint(sourcelib,source_T,source_P,suppresswarn=suppresswarn, composition = scomp)
   
    melt = get_melt(source_system)
    if mol(melt) ≈ 0
        println("No melt generated at ",source_T," °C and ", source_P, " bar" )
        return
    end

    close_meemum!(sourcelib)
    #Have to calculate melt conditions at host T and P
    meltlib = init_meemum(meltfile)
    melt_system = minimizepoint(meltlib,host_T,host_P,suppresswarn=suppresswarn,composition = melt.composition)

    μ = melt_system.composition[findchemical(melt_system.composition,equilib_component)].μ
    
    close_meemum!(meltlib)

    hostlib = init_meemum(hostfile)
    hcomp = getcompo(hostlib)
    if !isnothing(host_compo)
        hcomp = host_compo
    end
    
    host_system = minimizepoint(hostlib,host_T,host_P,μ1 = μ,suppresswarn=suppresswarn,composition=hcomp)

    close_meemum!(hostlib)
    return source_system,melt_system, host_system

end

function equilibrate_open_system(sourcefile,meltfile, hostfile,source_T,source_P, host_T,host_P, source_compo_start, source_compo_end; host_compo = nothing,equilib_component = "H2O",suppresswarn = false, steps = 100)
    #Should define another function for multiple compositions at source
    #So that init_meemum does not need to be run a billion times
    
    #Need to run these in sequence,first calc melt system at T and P
    #Then get uH2O from that melt and calc host system for given uH2O
    #Complicated by how datafiles are defined
    sourcelib = init_meemum(sourcefile)   

    source_compo_range = range(source_compo_start,source_compo_end,steps)
    source_systems = PetroSystem[]
   
    for compo in source_compo_range
        
        ssys= minimizepoint(sourcelib,source_T,source_P,suppresswarn=suppresswarn, composition = compo)
        push!(source_systems, ssys)
    end

    close_meemum!(sourcelib)

    meltlib = init_meemum(meltfile)
    melt_systems = PetroSystem[]
    for ssys in source_systems
        melt = get_melt(ssys)
        if mol(melt) ≈ 0
            println("No melt generated at ",source_T," °C and ", source_P, " bar with composition: ", ssys.composition)
            #Done to maintain array shapes
            push!(melt_systems,PetroSystem())
        else
            msys = minimizepoint(meltlib,host_T,host_P,suppresswarn=suppresswarn,composition = melt.composition)
            push!(melt_systems,msys)
        end
    end
    
    close_meemum!(meltlib)


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
                hsys = minimizepoint(hostlib,host_T,host_P,μ1 = μ,suppresswarn=suppresswarn,composition=hcomp)
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



# function equilibrate_open_system(sourcefile,meltfile, hostfile,source_T1, source_T2, source_P, host_T,host_P;source_compo = nothing, host_compo = nothing,equilib_component = "H2O",suppresswarn = false, steps = 100)
#     #Should define another function for multiple compositions at source
#     #So that init_meemum does not need to be run a billion times
    
#     #Need to run these in sequence,first calc melt system at T and P
#     #Then get uH2O from that melt and calc host system for given uH2O
#     #Complicated by how datafiles are defined
#     scomp = init_meemum(sourcefile)
#     if !isnothing(source_compo)
#         scomp = source_compo
#     end

#     Trange = range(source_T1, source_T2,steps)
#     source_systems = PetroSystem[]
#     melt_systems = PetroSystem[]
#     for T in Trange
#         ssys= minimizepoint(scomp,T,source_P,suppresswarn=suppresswarn)
#         push!(source_systems, ssys)
#         melt = get_melt(ssys)
#         if mol(melt) ≈ 0
#             println("No melt generated at ",source_T," °C and ", source_P, " bar with composition: ", compo)
#             #Done to maintain array shapes
#             push!(melt_systems,PetroSystem())

#         else
            
           
#             msys = minimizepoint(meltcompo,host_T,host_P,suppresswarn=suppresswarn)
#             push!(melt_systems,msys)
#         end
#     end
    
#     #Have to calculate melt conditions at host T and P
#     hcomp = init_meemum(hostfile)
#     if !isnothing(host_compo)
#         hcomp = host_compo
#     end
#     host_systems = PetroSystem[]
#     for msys in melt_systems
       
#         if length(msys.composition) > 0
#             μ = getchemical(msys,equilib_component).μ
#             hsys = minimizepoint(hcomp,host_T,host_P,μ1 = μ,suppresswarn=suppresswarn)
#             push!(host_systems,hsys)
#         else
#             push!(host_systems,PetroSystem())
#         end

#     end

#     return source_systems,melt_systems, host_systems

# end

end