# recalculation of diple operator
function recalculate_dipole_operators!(lab::LabSystem; new_objects::Bool=false)
    # check if new objects necessary
    if new_objects
        # new site dipole operators
        lab.dipoles_hor = Vector{DipoleOperator}(undef, length(lab.sites))
        lab.dipoles_ver = Vector{DipoleOperator}(undef, length(lab.sites))
        for s in 1:length(lab.sites)
            # transform the global beam parameters to the coordinate frame of site s
            site_q_in      = get_in_inner_coordinates(lab.sample, lab.q_in)
            site_q_out     = get_in_inner_coordinates(lab.sample, lab.q_out)
            site_eps_in    = get_in_inner_coordinates(lab.sites[s], get_in_inner_coordinates(lab.sample, lab.epsilon_in))
            site_eps_out_h = get_in_inner_coordinates(lab.sites[s], get_in_inner_coordinates(lab.sample, lab.epsilon_out_hor))
            site_eps_out_v = get_in_inner_coordinates(lab.sites[s], get_in_inner_coordinates(lab.sample, lab.epsilon_out_ver))
            site_position  = lab.sites[s].position
            # create new dipole operators
            lab.dipoles_hor[s] = DipoleOperator(
                        edge = lab.edge,
                        site = s,
                        site_position = site_position,
                        eps_in  = site_eps_in,
                        eps_out = site_eps_out_h,
                        q_in    = site_q_in,
                        q_out   = site_q_out
                    )
            lab.dipoles_ver[s] = DipoleOperator(
                        edge = lab.edge,
                        site = s,
                        site_position = site_position,
                        eps_in  = site_eps_in,
                        eps_out = site_eps_out_v,
                        q_in    = site_q_in,
                        q_out   = site_q_out
                    )
            if lab.spin_quantization == :local
                lab.dipoles_hor[s].spin_quantization = CoordinateFrame()
                lab.dipoles_ver[s].spin_quantization = CoordinateFrame()
            elseif lab.spin_quantization == :sample
                lab.dipoles_hor[s].spin_quantization = lab.sites[s]
                lab.dipoles_ver[s].spin_quantization = lab.sites[s]
            else
                error("Unknown spin quantization axis "*string(lab.spin_quantization))
            end
            if lab.dipoles_ver[s].edge != lab.edge
                lab.dipoles_ver[s].edge = lab.edge
                recalculate!(lab.dipoles_ver[s], true, true)
            end
            if lab.dipoles_hor[s].edge != lab.edge
                lab.dipoles_hor[s].edge = lab.edge
                recalculate!(lab.dipoles_hor[s], true, true)
            end
        end
        # build sums of diple operators
        N = length(basis(lab.hamiltonian)[1].occupation)
        basis_dipole = getMultiParticleBasis(getMultiSiteBasis(getT2GBasisXYZ(), length(lab.sites)), N)
        lab.dipole_hor = ProjectorOperator( sum([
            MPGeneralizedSPOperator(
                basis_dipole,
                ProjectorOperator(
                    dh, basis_dipole.single_particle_basis
                )
            )
            for dh in lab.dipoles_hor  ]), basis(lab.hamiltonian))
        lab.dipole_ver = ProjectorOperator( sum([
            MPGeneralizedSPOperator(
                basis_dipole,
                ProjectorOperator(
                    dv, basis_dipole.single_particle_basis
                )
            )
            for dv in lab.dipoles_ver  ]), basis(lab.hamiltonian))
    else
        # update the dipole operators
        for s in 1:length(lab.sites)
            # transform the global beam parameters to the coordinate frame of site s
            site_q_in      = get_in_inner_coordinates(lab.sample, lab.q_in)
            site_q_out     = get_in_inner_coordinates(lab.sample, lab.q_out)
            site_eps_in    = get_in_inner_coordinates(lab.sites[s], get_in_inner_coordinates(lab.sample, lab.epsilon_in))
            site_eps_out_h = get_in_inner_coordinates(lab.sites[s], get_in_inner_coordinates(lab.sample, lab.epsilon_out_hor))
            site_eps_out_v = get_in_inner_coordinates(lab.sites[s], get_in_inner_coordinates(lab.sample, lab.epsilon_out_ver))
            site_position  = lab.sites[s].position
            # update diple operators
            lab.dipoles_hor[s].position = site_position
            lab.dipoles_hor[s].eps_in   = site_eps_in
            lab.dipoles_hor[s].eps_out  = site_eps_out_h
            lab.dipoles_hor[s].q_in     = site_q_in
            lab.dipoles_hor[s].q_out    = site_q_out
            lab.dipoles_ver[s].position = site_position
            lab.dipoles_ver[s].eps_in   = site_eps_in
            lab.dipoles_ver[s].eps_out  = site_eps_out_v
            lab.dipoles_ver[s].q_in     = site_q_in
            lab.dipoles_ver[s].q_out    = site_q_out
            if lab.spin_quantization == :local
                lab.dipoles_hor[s].spin_quantization = CoordinateFrame()
                lab.dipoles_ver[s].spin_quantization = CoordinateFrame()
            elseif lab.spin_quantization == :sample
                lab.dipoles_hor[s].spin_quantization = lab.sites[s]
                lab.dipoles_ver[s].spin_quantization = lab.sites[s]
            else
                error("Unknown spin quantization axis "*string(lab.spin_quantization))
            end
            if lab.dipoles_ver[s].edge != lab.edge
                lab.dipoles_ver[s].edge = lab.edge
                recalculate!(lab.dipoles_ver[s], true, true)
            end
            if lab.dipoles_hor[s].edge != lab.edge
                lab.dipoles_hor[s].edge = lab.edge
                recalculate!(lab.dipoles_hor[s], true, true)
            end
        end
    end
    # recalculate the dipole operator top down
    recalculate!(lab.dipole_hor, true, false)
    recalculate!(lab.dipole_ver, true, false)
end
# recalculate the hamiltonian
function recalculate_hamiltonian!(ls :: LabSystem, basis_change::Bool=true, rediagonalize::Bool=true)
    # recalculate hamiltonian
    recalculate!(ls.hamiltonian, true, basis_change)
    # maybe rediagonalize
    if rediagonalize
        ls.eigensys = eigensystem(ls.hamiltonian)
    end
end

# possibly recalculate the matrix representation
function recalculate!(ls :: LabSystem, basis_change::Bool=true, rediagonalize::Bool=true)
    # maybe recalculate dipole operators
    recalculate_dipole_operators!(ls, new_objects=basis_change)
    # recalculate hamiltonian
    recalculate_hamiltonian!(ls, basis_change, rediagonalize)
end


function basis(ls::LabSystem)
    return basis(ls.hamiltonian)
end


################################################################################################################################


# obtaining a spectrum
function get_spectrum(
            ls :: LabSystem,
            args...
            ;
            kwargs...
        ) :: Spectrum where {B <: AbstractBasis}

    # construct spectrum
    return get_spectrum(ls.eigensys, ls.dipole_hor, args...; kwargs...) + get_spectrum(ls.hamiltonian, ls.dipole_ver, args...; kwargs...)
end





function get_dq_eps(line :: AbstractString; factor :: Real = 1.0)
    parts = split(line, "\t")
    return (
        [
            Meta.eval(Meta.parse(parts[1])),
            Meta.eval(Meta.parse(parts[2])),
            Meta.eval(Meta.parse(parts[3]))
        ] .* factor, [
            Meta.eval(Meta.parse(parts[4])),
            Meta.eval(Meta.parse(parts[5])),
            Meta.eval(Meta.parse(parts[6]))
        ]
    )
end
function get_dq_eps(fn :: AbstractString, linenumber :: Integer; factor :: Real = 1.0)
    # open file and extract values
    f = open(fn, "r")
    lines = strip.(readlines(f))
    lines = [l for l in lines if !startswith(l, "#") && isdigit(l[1])]
    close(f)
    parts = split(lines[linenumber], "\t")
    return (
        [
            Meta.eval(Meta.parse(parts[1])),
            Meta.eval(Meta.parse(parts[2])),
            Meta.eval(Meta.parse(parts[3]))
        ] .* factor, [
            Meta.eval(Meta.parse(parts[4])),
            Meta.eval(Meta.parse(parts[5])),
            Meta.eval(Meta.parse(parts[6]))
        ]
    )
end




################################################################################################################################



# set a parameter (returns (found parameter?, changed matrix?))
function set_parameter!(ls :: LabSystem, parameter :: Symbol, value; print_result::Bool=false, recalculate::Bool=true, site=:all, relative_to::Symbol=:sample, rediagonalize::Bool=true, kwargs...)
    # check if it can be set in 1
    if site == :all
        # define dfault parameters
        found_param, changed_matrix = false, false
        # go through all sites
        for s in 1:length(ls.sites)
            # set the value on site s
            found_param_s, changed_matrix_s = set_parameter!(ls, parameter, value; print_result=print_result, recalculate=false, site=s, relative_to=relative_to, kwargs...)
            # process parameters
            found_param_s, changed_matrix_s = (found_param || found_param_s), (changed_matrix || changed_matrix_s)
        end
        # recalculate everything
        recalculate_hamiltonian!(ls, false, found_param && rediagonalize)
        # return
        return (found_param, changed_matrix)
    elseif typeof(site) <: Integer
        # transform maybe to local axis
        if typeof(value) <: Vector{<:Any} && length(value) == 3
            if relative_to == :site
                # no transformation needed
            elseif relative_to == :sample
                # transform to site coordinates
                value = get_in_inner_coordinates(ls.sites[site], value)
            elseif relative_to == :global || relative_to == :lab
                # transform to site coordinates
                value = get_in_inner_coordinates(ls.sites[site], get_in_inner_coordinates(ls.sample, value))
            else
                error("relative orientation unknown: $(relative_to)")
            end
        end
        # set parameter
        found_param, changed_matrix = set_parameter!(ls.hamiltonian, parameter, value, print_result=print_result, recalculate=true; site=site, kwargs...)
        return (found_param, changed_matrix)
    else
        error("specify a corret site: $(site)")
    end
end



