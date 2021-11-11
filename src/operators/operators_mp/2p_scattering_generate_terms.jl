# generate a term of the spin flip exchange form for site i (ELECTRON)
function generate2PScatteringElectronInteractionSpinFlip(basis :: MPB, site :: Int64, prefactor :: Real) where {
        N, SPSSBS<:BasisStateXYZ, SPBS <: SPMSBasisState{SPSSBS},
        MPB <: MPBasis{N,SPBS}
    }

    # select all single particle states that are on site site
    sp_state_indices = Int64[i for i in 1:length(basis.single_particle_basis) if basis.single_particle_basis[i].site == site]
    sp_states = SPBS[basis.single_particle_basis[i] for i in sp_state_indices]

    # generate a new operator
    op = MPElectron2PScatteringOperator(basis, Tuple{Tuple{Int64,Int64},Tuple{Int64,Int64}}[], prefactor)

    # try generating all pairs ((a,b),(c,d))
    # s.t. orb(a)=orb(b) != orb(c)=orb(d) , spin(a) < spin(b), spin(c) > spin(d), spin(c) = spin(b)
    # meaning: b->a, d->c scattering

    # find scattering event 1
    for a in 1:length(sp_states)
    for b in 1:length(sp_states)
        # check if it is orb(a) = orb(b) but spin(a) != spin(b)
        if sp_states[a].state.orbital == sp_states[b].state.orbital && sp_states[a].state.ms > sp_states[b].state.ms
            # find scattering event 2
            for c in 1:length(sp_states)
            for d in 1:length(sp_states)
                if sp_states[c].state.orbital == sp_states[d].state.orbital && sp_states[c].state.orbital != sp_states[b].state.orbital &&
                    sp_states[c].state.ms < sp_states[d].state.ms && sp_states[c].state.ms == sp_states[b].state.ms
                    push!(op.interacting_orbitals, ((sp_state_indices[a],sp_state_indices[b]),(sp_state_indices[c],sp_state_indices[d])))
                end
            end
            end
        end
    end
    end

    # recalculate the operator
    recalculate!(op)

    # return the finished operator
    return op
end
export generate2PScatteringElectronInteractionSpinFlip

# generate a term of the spin conserving exchange form for site i (ELECTRON)
function generate2PScatteringElectronInteractionSpinConserve(basis :: MPB, site :: Int64, prefactor :: Real) where {
        N, SPSSBS<:BasisStateXYZ, SPBS <: SPMSBasisState{SPSSBS},
        MPB <: MPBasis{N,SPBS}
    }

    # select all single particle states that are on site site
    sp_state_indices = Int64[i for i in 1:length(basis.single_particle_basis) if basis.single_particle_basis[i].site == site]
    sp_states = SPBS[basis.single_particle_basis[i] for i in sp_state_indices]

    # generate a new operator
    op = MPElectron2PScatteringOperator(basis, Tuple{Tuple{Int64,Int64},Tuple{Int64,Int64}}[], prefactor)

    # try generating all pairs ((a,b),(c,d))
    # s.t. orb(a)=orb(c)!=orb(b)=orb(d) , spin(a)=spin(b)<spin(c)=spin(d)
    # meaning: b->a, d->c scattering

    # find scattering event 1
    for a in 1:length(sp_states)
    for b in 1:length(sp_states)
        # check if it is orb(a)!= orb(b) but spin(a)=spin(b)
        if sp_states[a].state.orbital != sp_states[b].state.orbital && sp_states[a].state.ms == sp_states[b].state.ms
            # find scattering event 2
            for c in 1:length(sp_states)
            for d in 1:length(sp_states)
                if sp_states[c].state.orbital == sp_states[a].state.orbital && sp_states[d].state.orbital == sp_states[b].state.orbital &&
                    sp_states[c].state.ms > sp_states[a].state.ms && sp_states[c].state.ms == sp_states[d].state.ms
                    push!(op.interacting_orbitals, ((sp_state_indices[a],sp_state_indices[b]),(sp_state_indices[c],sp_state_indices[d])))
                end
            end
            end
        end
    end
    end

    # recalculate the operator
    recalculate!(op)

    # return the finished operator
    return op
end
export generate2PScatteringElectronInteractionSpinConserve

# generate a term of the spin flip exchange form for site i (ELECTRON)
function generate2PScatteringHoleInteractionSpinFlip(basis :: MPB, site :: Int64, prefactor :: Real) where {
        N, SPSSBS<:BasisStateXYZ, SPBS <: SPMSBasisState{SPSSBS},
        MPB <: MPBasis{N,SPBS}
    }

    # select all single particle states that are on site site
    sp_state_indices = Int64[i for i in 1:length(basis.single_particle_basis) if basis.single_particle_basis[i].site == site]
    sp_states = SPBS[basis.single_particle_basis[i] for i in sp_state_indices]

    # generate a new operator
    op = MPHole2PScatteringOperator(basis, Tuple{Tuple{Int64,Int64},Tuple{Int64,Int64}}[], prefactor)

    # try generating all pairs ((a,b),(c,d))
    # s.t. orb(a)=orb(b) != orb(c)=orb(d) , spin(a) < spin(b), spin(c) > spin(d), spin(c) = spin(b)
    # meaning: b->a, d->c scattering

    # find scattering event 1
    for a in 1:length(sp_states)
    for b in 1:length(sp_states)
        # check if it is orb(a) = orb(b) but spin(a) != spin(b)
        if sp_states[a].state.orbital == sp_states[b].state.orbital && sp_states[a].state.ms > sp_states[b].state.ms
            # find scattering event 2
            for c in 1:length(sp_states)
            for d in 1:length(sp_states)
                if sp_states[c].state.orbital == sp_states[d].state.orbital && sp_states[c].state.orbital != sp_states[b].state.orbital &&
                    sp_states[c].state.ms < sp_states[d].state.ms && sp_states[c].state.ms == sp_states[b].state.ms
                    push!(op.interacting_orbitals, ((sp_state_indices[a],sp_state_indices[b]),(sp_state_indices[c],sp_state_indices[d])))
                end
            end
            end
        end
    end
    end

    # recalculate the operator
    recalculate!(op)

    # return the finished operator
    return op
end
export generate2PScatteringHoleInteractionSpinFlip


# generate a term of the spin conserving exchange form for site i (ELECTRON)
function generate2PScatteringHoleInteractionSpinConserve(basis :: MPB, site :: Int64, prefactor :: Real) where {
        N, SPSSBS<:BasisStateXYZ, SPBS <: SPMSBasisState{SPSSBS},
        MPB <: MPBasis{N,SPBS}
    }

    # select all single particle states that are on site site
    sp_state_indices = Int64[i for i in 1:length(basis.single_particle_basis) if basis.single_particle_basis[i].site == site]
    sp_states = SPBS[basis.single_particle_basis[i] for i in sp_state_indices]

    # generate a new operator
    op = MPHole2PScatteringOperator(basis, Tuple{Tuple{Int64,Int64},Tuple{Int64,Int64}}[], prefactor)

    # try generating all pairs ((a,b),(c,d))
    # s.t. orb(a)=orb(c)!=orb(b)=orb(d) , spin(a)=spin(b)<spin(c)=spin(d)
    # meaning: b->a, d->c scattering

    # find scattering event 1
    for a in 1:length(sp_states)
    for b in 1:length(sp_states)
        # check if it is orb(a)!= orb(b) but spin(a)=spin(b)
        if sp_states[a].state.orbital != sp_states[b].state.orbital && sp_states[a].state.ms == sp_states[b].state.ms
            # find scattering event 2
            for c in 1:length(sp_states)
            for d in 1:length(sp_states)
                if sp_states[c].state.orbital == sp_states[a].state.orbital && sp_states[d].state.orbital == sp_states[b].state.orbital &&
                    sp_states[c].state.ms > sp_states[a].state.ms && sp_states[c].state.ms == sp_states[d].state.ms
                    push!(op.interacting_orbitals, ((sp_state_indices[a],sp_state_indices[b]),(sp_state_indices[c],sp_state_indices[d])))
                end
            end
            end
        end
    end
    end

    # recalculate the operator
    recalculate!(op)

    # return the finished operator
    return op
end
export generate2PScatteringHoleInteractionSpinConserve
