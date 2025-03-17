module Migmatites

export
    getmelt,
    equilibrate_open_system,
    equilibrate_closed_system,
    balance_component,
    gridcalc,
    opensystem_melting_path,
    melt_cycle_model_open

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
    return getphase(system,"melt")[1]
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
function equilibrate_closed_system(meemumlib,sourcecompo,hostcompo,source_T,source_P, host_T,host_P, melt_proportion;equilib_component = "H2O",suppresswarn = false, phasefunc = [],μthreshold = 5.0,xthreshold = 1E-5)
    #Set up: get the host rock composition
   
    source_system = minimizepoint(meemumlib,source_T,source_P,suppresswarn=suppresswarn, composition = sourcecompo, phasefunc = phasefunc)
   
    melt = getmelt(source_system)
    if mol(melt) ≈ 0
        println("No melt generated at ",source_T," °C and ", source_P, " bar" )
        return
    end

    
    
    meltcompo = melt.composition*100
    lowerlimit = -molfrac(hostcompo,equilib_component)/melt_proportion
    
    upperlimit = molfrac(meltcompo,equilib_component)
    @show getchemical(sourcecompo,equilib_component).mol
    x_change = dumbsolve(x -> closed_system_optim(meemumlib,meltcompo,hostcompo,host_T,host_P,
                                            melt_proportion,x,equilib_component=equilib_component,
                                            suppresswarn=suppresswarn,phasefunc=phasefunc),lowerlimit,upperlimit,ftol = μthreshold, xtol = xthreshold)
    
    
    
    # @show x_change
    newmeltcompo, newhostcompo = balance_component(meltcompo,hostcompo,equilib_component,x_change,melt_proportion)
  
   
    melt_system = minimizepoint(meemumlib,host_T,host_P,composition = newmeltcompo, suppresswarn = suppresswarn, phasefunc = phasefunc)
    host_system = minimizepoint(meemumlib,host_T,host_P,composition = newhostcompo, suppresswarn = suppresswarn, phasefunc = phasefunc)
    # chemi = findchemical(melt_system.composition,equilib_component)
    # chemj = findchemical(host_system.composition,equilib_component)
    # μ_melt = melt_system.composition[chemi].μ
    # μ_host = host_system.composition[chemj].μ


    return source_system,melt_system, host_system

end

function closed_system_optim(meemumlib,meltcompo,hostcompo,host_T,host_P,melt_proportion,x_change;equilib_component = "H2O",suppresswarn=false,phasefunc = [])
    

    newmeltcompo, newhostcompo = balance_component(meltcompo,hostcompo,equilib_component,x_change,melt_proportion)

    melt_system = minimizepoint(meemumlib,host_T,host_P,composition = newmeltcompo, suppresswarn = suppresswarn, phasefunc = phasefunc)
    host_system = minimizepoint(meemumlib,host_T,host_P,composition = newhostcompo, suppresswarn = suppresswarn, phasefunc = phasefunc)

    chemi = findchemical(melt_system.composition,equilib_component)
    chemj = findchemical(host_system.composition,equilib_component)
    if chemi == 0 || chemj == 0
        return NaN
    end
    if melt_system.composition[chemi].mol <= 0.0
        return -Inf
    end
    if  host_system.composition[chemj].mol <= 0.0
        return Inf
    end
   
    μ_melt = melt_system.composition[chemi].μ
    μ_host = host_system.composition[chemj].μ

    return μ_melt-μ_host
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


function dumbsolve(f, intervalstart,intervalend;numchecks = 10, xtol = 1e-8, ftol = 1e-8, maxiter = 1000)

    checkrange = Vector(range(intervalstart,intervalend,numchecks))

    res = f.(checkrange)

    deleteat!(checkrange,isnan.(res))
    deleteat!(res,isnan.(res))
    minpos = Inf
    posi = 0
    maxneg = -Inf
    negi = 0
    for i in 1:lastindex(res)

        if res[i] < minpos && res[i] >= 0
            minpos = res[i]
            posi = i
        end
        if res[i] > maxneg && res[i] <= 0
            maxneg = res[i]
            negi = i
        end
    end
    
    if (posi == 0 && negi == 0) || (posi == 1 && negi == numchecks)
        return dumbsolve(f,intervalstart,intervalend, numchecks = numchecks*2, xtol = xtol, ftol = ftol, maxiter = maxiter-1)
    end

    if posi == 0
        newstart = intervalstart
        newend = checkrange[1]
        if negi == lastindex(res)
            newstart = checkrange[end]
            newend = intervalend
        end
        return dumbsolve(f,newstart,newend, numchecks = numchecks, xtol = xtol, ftol = ftol, maxiter = maxiter-1)
    end

    if negi == 0
        newstart = intervalstart
        newend = checkrange[1]
        if posi == lastindex(res)
            newstart = checkrange[end]
            newend = intervalend
        end
        return dumbsolve(f,newstart,newend, numchecks = numchecks, xtol = xtol, ftol = ftol, maxiter = maxiter-1)
    end
   
    if minpos <= ftol
        return checkrange[posi]
    end
    if maxneg >= -ftol
        return checkrange[negi]
    end
    if abs(checkrange[posi] - checkrange[negi]) <= xtol
        return checkrange[posi]
    end
    if maxiter == 0
        @warn "Did not converge, returning best result"
        if abs(minpos) > abs(maxneg)
            println("Minimum = ",abs(minpos))
            return checkrange[posi]
        else
            println("Minimum = ",abs(maxneg))
            return checkrange[negi]
        end
    end
    if numchecks >= 2000
        @warn "Too many value checks, returning best result"
        if abs(minpos) > abs(maxneg)
            println("Minimum = ",abs(minpos))
            return checkrange[posi]
        else
            println("Minimum = ",abs(maxneg))
            return checkrange[negi]
        end
    end
    newstart = checkrange[negi]
    newend = checkrange[posi]
    if checkrange[negi] > checkrange[posi]
        newstart = checkrange[posi]
        newend = checkrange[negi]
    end
  
    return dumbsolve(f,newstart,newend, numchecks = numchecks, xtol = xtol, ftol = ftol, maxiter = maxiter-1)

end

function gridcalc(f, xrange, yrange; xnodes = 10, ynodes = 10)

    xs = LinRange(xrange[1],xrange[2],xnodes)
    ys = LinRange(yrange[1],yrange[2],ynodes)
    results = []
    for x in xs
        for y in ys
            results = push!(results,f(x,y))
        end
    end

    return results
end

function opensystem_melting_path(meemumlib,composition,t_path,p_path;vol_threshold = 0.07,vol_residual = 0.01, numsteps = 100,suppresswarn=false,phasefunc = [])

    if length(t_path) != length(p_path)
        throw(ErrorException("Temperature and pressure path inputs must be the same length"))
    end
    t_steps = []
    p_steps = []
    interval_steps = trunc(Int,numsteps/(length(t_path)-1))
    for i in 2:lastindex(t_path)
        t_steps = [t_steps;range(t_path[i-1],t_path[i],interval_steps)]
        p_steps = [p_steps;range(p_path[i-1],p_path[i],interval_steps)]
    end

    current_compo = composition
    melts = []
    restite_compos = []
    system_steps = []
    extract_T = []
    extract_P = []
    for i in 1:lastindex(t_steps)
        system = minimizepoint(meemumlib,t_steps[i],p_steps[i],composition = current_compo, suppresswarn = suppresswarn, phasefunc = phasefunc)
        push!(system_steps,system)
        melt_volprop = get_volprop(system,"melt")
        if melt_volprop >= vol_threshold
            meltphase = getmelt(system)
            percentremoved = (melt_volprop - vol_residual)/melt_volprop
            vol_removed = percentremoved*meltphase.vol
            mol_removed = meltphase.mol/meltphase.vol * vol_removed
            newsyscompo = current_compo .- meltphase.composition*mol_removed
            current_compo = newsyscompo
            push!(melts,meltphase)
            push!(restite_compos, current_compo)
            push!(extract_T,t_steps[i])
            push!(extract_P,p_steps[i])
        end
    end

    return melts, restite_compos, system_steps, extract_T, extract_P

end

# function opensystem_melting_path(meemumlib,composition,t_path,x_range, pressure;vol_threshold = 0.07,vol_residual = 0.01, numsteps = 100,suppresswarn=false,phasefunc = [],x_component = "H2O")

#     t_steps = []
#     x_step = (x_range[2] - x_range[1])/numsteps
#     interval_steps = trunc(Int,numsteps/length(t_path))
#     for i in 2:lastindex(t_path)
#         t_steps = [t_steps;range(t_path[i-1],t_path[i],interval_steps)]
#     end

#     current_compo = composition
#     melts = []
#     restite_compos = []
#     system_steps = []
#     extract_T = []
#     extract_X = []
#     for i in 1:lastindex(t_steps)
#         system = minimizepoint(meemumlib,t_steps[i],pressure,composition = current_compo, suppresswarn = suppresswarn, phasefunc = phasefunc)
#         push!(system_steps,system)
#         melt_volprop = get_volprop(system,"melt")
#         if melt_volprop >= vol_threshold
#             meltphase = getmelt(system)
#             percentremoved = (melt_volprop - vol_residual)/melt_volprop
#             vol_removed = percentremoved*meltphase.vol
#             mol_removed = meltphase.mol/meltphase.vol * vol_removed
#             newsyscompo = current_compo .- meltphase.composition*mol_removed
#             current_compo = newsyscompo
#             push!(melts,meltphase)
#             push!(restite_compos, current_compo)
#             push!(extract_T,t_steps[i])
#             push!(extract_X,molfrac(current_compo,x_component))
#         end

#         x = molfrac(current_compo,x_component)
       
#         x_new = x+x_step
        
#         #Reverse molfrac normalization
#         mol_new = x_new*sum_mols(current_compo)
        
#         current_compo = change_list_component(current_compo,mol_new,x_component)
#     end

#     return melts, restite_compos, system_steps, extract_T, extract_X

# end
function equilibrate_melt_opensys(hostfile,sourcefile,melt_compo,host_T,host_P;host_compo = nothing,equilib_component = "H2O",suppresswarn = false, phasefunc = [])
    sourcelib = init_meemum(sourcefile)
    melt_system = minimizepoint(sourcelib,host_T,host_P,suppresswarn=suppresswarn,composition = melt_compo, phasefunc = phasefunc)

    μ = melt_system.composition[findchemical(melt_system.composition,equilib_component)].μ
    
    close_meemum!(sourcelib)

    hostlib = init_meemum(hostfile)
    hcomp = getcompo(hostlib)
    if !isnothing(host_compo)
        hcomp = host_compo
    end
    
    host_system = minimizepoint(hostlib,host_T,host_P,μ1 = μ,suppresswarn=suppresswarn,composition=hcomp, phasefunc = phasefunc)

    close_meemum!(hostlib)
    return melt_system, host_system
end

function equilibrate_melts_opensys(hostfile,sourcefile,melt_compos,host_T,host_P;host_compo = nothing,equilib_component = "H2O",suppresswarn = false, phasefunc = [])
    sourcelib = init_meemum(sourcefile)
    melt_systems = []
    μs = []
    for compo in melt_compos
        melt = minimizepoint(sourcelib,host_T,host_P,suppresswarn=suppresswarn,composition = compo, phasefunc = phasefunc)
        push!(melt_systems,melt)
        μ = melt.composition[findchemical(melt.composition,equilib_component)].μ
        push!(μs, μ)
    end
    
    
    close_meemum!(sourcelib)

    hostlib = init_meemum(hostfile)
    hcomp = getcompo(hostlib)
    if !isnothing(host_compo)
        hcomp = host_compo
    end
    host_systems = []
    for μ in μs
        host = minimizepoint(hostlib,host_T,host_P,μ1 = μ,suppresswarn=suppresswarn,composition=hcomp, phasefunc = phasefunc)
        push!(host_systems,host)
    end
    close_meemum!(hostlib)
    return melt_systems, host_systems
end

function melt_cycle_model_open(hostfile, sourcefile, host_T, host_P, t_path, p_path;host_compo = nothing,equilib_component = "H2O",suppresswarn = false, phasefunc = [],vol_threshold = 0.07,vol_residual = 0.01, numsteps = 100,source_compo=nothing)
    sourcelib = init_meemum(sourcefile)
    scomp = getcompo(sourcelib)
    if !isnothing(source_compo)
        scomp = source_compo
    end
    melts_i, restite_compos, system_steps, extract_T, extract_P = opensystem_melting_path(sourcelib,scomp,t_path,p_path,
                                                    vol_threshold = vol_threshold, vol_residual = vol_residual, 
                                                    numsteps = numsteps, suppresswarn=suppresswarn,phasefunc=phasefunc)
    close_meemum!(sourcelib)
    if length(melts_i) > 0
       
        meltcompos = getcompo.(melts_i)
  
        melts_f, hosts = equilibrate_melts_opensys(hostfile,sourcefile,meltcompos,host_T,host_P;
                                                    host_compo = host_compo, equilib_component = equilib_component,
                                                    suppresswarn = suppresswarn, phasefunc = phasefunc)

        return melts_f, hosts, restite_compos, system_steps, extract_T, extract_P
    else
        return [], [], [], system_steps, [], []
    end
end

end