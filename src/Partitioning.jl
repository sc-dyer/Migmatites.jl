#Sub-module for partition coefficients
const ZIRCON_ZR = 497657
const ZR = TraceElement("Zr",91.224,NaN)
const ZRN_ρ = 4650
"""
$(TYPEDSIGNATURES)
Calculates zircon saturation of a melt (i.e. maximum concentration of Zr in melt) following Watson & Harrison 1983.
Temperature is expected in celsius
"""
function zircon_saturation(melt,temp)
    temp += 273.15    
    #Calculate M parameter based on weight percentages of 
    function calc_m(melt)
        meltcompo = melt.composition
        na = mol(getchemical(meltcompo,"Na2O"))*2
        k = mol(getchemical(meltcompo,"K2O"))*2
        ca = mol(getchemical(meltcompo,"CaO"))
        al = mol(getchemical(meltcompo,"Al2O3"))*2
        si = mol(getchemical(meltcompo,"SiO2"))
        total = na + k + ca + al + si
        na = na/total
        k = k/total
        ca = ca/total
        al = al/total
        si = si/total
        return (na + k + 2*ca)/(al*si)
    end
    
    m = calc_m(melt)
    lnD = -3.8-(0.85*(m-1))+12900/temp
    
    meltZr = ZIRCON_ZR/exp(lnD)

    return meltZr
end

"""
$(TYPEDSIGNATURES)
Import Kd values to a dictionary, following a csv file template with two columns with the first column being the phase
name and the second column being the Kd value (mineral/melt). Mineral names must match the names of phases in the PetroSystem, 
and every phase present with melt must have an associated Kd value or an error will be thrown.
"""
function import_Kd_values(filename)
    return Dict(CSV.File(filename))
end

"""
$(TYPEDSIGNATURES)
Calculates the bulk partition coefficient for a given system and partition coefficients
"""
function bulk_D(system,kd_dict)
    
    bulk_d = 0

    for phase in system.phases
        if !contains(lowercase(phase.name),"melt")
            bulk_d += massfrac(phase,system)*kd_dict[phase.name]
        end
    end

    return bulk_d

end

"""
$(TYPEDSIGNATURES)
Calculates the concentration of Zr in each phase of a system with temp in degreese Celsius
"""
function calc_Zr!(system,temp,kd_dict)
    
    melt = getmelt(system)
    if melt.mol == 0
       
        for phase in system.phases
            push!(phase.traceelements,ZR)
        end
    else
        systemZr = concentration(getchemical(system,"Zr"))
        bulkD = bulk_D(system,kd_dict)
        melt_massfrac = massfrac(melt,system)
        cl_co = 1/(bulkD*(1-melt_massfrac)+melt_massfrac)
        modelledZr = cl_co*systemZr
        meltZr = modelledZr
        zr_sat = zircon_saturation(melt,temp)
        if modelledZr > zr_sat
            meltZr = zr_sat
        end

        push!(melt.traceelements,TraceElement(ZR,meltZr))

        for phase in system.phases
            if !contains(lowercase(phase.name),"melt")
                phaseZr = kd_dict[phase.name]*meltZr
                push!(phase.traceelements,TraceElement(ZR,phaseZr))
            end
        end
    end

end

"""
$(TYPEDSIGNATURES)
Calculates the mass fraction of zircon in a given system that has already had Zr concentration calculated
"""
function calc_zircon!(system)
    totalZr = 0.0

    for phase in system.phases
        totalZr += massfrac(phase,system)*concentration(getchemical(phase,ZR))
    end

    excessZr = concentration(getchemical(system,ZR)) - totalZr

    zircon_compo = [Component("ZrO2",MOLAR_MASSES["ZrO2"],1.0),Component("SiO2",MOLAR_MASSES["SiO2"],1.0)]
    
    zrn_mass_frac =excessZr/ZIRCON_ZR
    zrn_mass = zrn_mass_frac*system.mass
    zrn_molarmass = sum_mass(zircon_compo)
    zrn_vol = ((zrn_mass/10^3)/ZRN_ρ)/10^5
    zrn_mol = zrn_mass/zrn_molarmass
    zrn_phase = Phase(name = "Zrn",
                    composition = zircon_compo,
                    molarmass = zrn_molarmass,
                    mol = zrn_mol,
                    ρ = ZRN_ρ,
                    vol = zrn_vol,
                    Vmol = zrn_vol/zrn_mol,
                    traceelements = [TraceElement(ZR,ZIRCON_ZR)],
                    mass = zrn_mass)

    push!(system.phases,zrn_phase)
    return zrn_mass_frac
end
