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