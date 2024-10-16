module Migmatites

export
    get_melt,
    equilibrate_open_system

using
    DocStringExtensions,
    JPerpleX,
    PetroBase

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

Returns the melt phase in a 'PetroSystem' if one exists, throws an error if it does not
"""
function get_melt(system)
    for phase in system.phases
        if contains(lowercase(phase.name),"melt")
            return phase
        end
    end
    throw(ErrorException("System does not contain melt"))
end

"""
$(TYPEDSIGNATURES)

This is a two step function which calculates the system and phase properties of the melt source rock at 'source_T' and
'source_P', then use the chemical potential of `equilib_component` to calculate the system and phase properties of
the host rock at `host_T` and `host_P`. If `meltsource_compo` or `host_compo` are not provided, it will default to
the compositions defined in the given dat files: `meltsourcefile` and `hostfile`. For this function to work as expected
your meltsourcefile must have only P and T as independent variables, and your hostfile must have P, T and μ of 
one component as independent variables.
"""
function equilibrate_open_system(meltsourcefile, hostfile,source_T,source_P, host_T,host_P;meltsource_compo = nothing, host_compo = nothing,equilib_component = "H2O",suppresswarn = false)

    #Need to run these in sequence,first calc melt system at T and P
    #Then get uH2O from that melt and calc host system for given uH2O
    #Complicated by how datafiles are defined
    mscomp = init_meemum(meltsourcefile)
    if !isnothing(meltsource_compo)
        mscomp = meltsource_compo
    end
    melt_system = minimizepoint(mscomp,source_T,source_P,suppresswarn=suppresswarn)
    melt = Phase(name="melt",composition = Component[])
    try
        melt = get_melt(melt_system)
    catch
        println("No melt generated at ",source_T," °C and ", source_P, " bar" )
        return
    end
    μ = melt.composition[findchemical(melt.composition,equilib_component)].μ
    @show μ
    hcomp = init_meemum(hostfile)

    if !isnothing(host_compo)
        hcomp = host_compo
    end

    host_system = minimizepoint(hcomp,host_T,host_P,μ1 = μ,suppresswarn=suppresswarn)

    return melt_system, host_system

end



end
