# generate a term of the form n_aup n_adown for site i (ELECTRON)
function generateDensityDensityElectronInteractionSameOrbital(basis :: MPB, site :: Int64, prefactor :: Real) where {
        N, SPSSBS<:BasisStateXYZ, SPBS <: SPMSBasisState{SPSSBS},
        MPB <: MPBasis{N,SPBS}
    }

    # select all single particle states that are on site site
    sp_state_indices = Int64[i for i in 1:length(basis.single_particle_basis) if basis.single_particle_basis[i].site == site]
    sp_states = SPBS[basis.single_particle_basis[i] for i in sp_state_indices]

    # generate a new operator
    op = MPElectronDensityDensityOperator(basis, Tuple{Int64,Int64}[], prefactor)

    # try generating all pairs (a,b) s.t. orb(a) = orb(b) but spin(a) != spin(b)
    for a in 1:length(sp_states)
    for b in 1:length(sp_states)
        # check if it is orb(a) = orb(b) but spin(a) != spin(b)
        if sp_states[a].state.orbital == sp_states[b].state.orbital && sp_states[a].state.ms > sp_states[b].state.ms
            push!(op.interacting_orbitals, (sp_state_indices[a],sp_state_indices[b]))
        end
    end
    end

    # recalculate the operator
    recalculate!(op)

    # return the finished operator
    return op
end
export generateDensityDensityElectronInteractionSameOrbital

# generate a term of the form n_aup n_apup + n_adown n_apdown for site i (ELECTRON)
function generateDensityDensityElectronInteractionSameSpin(basis :: MPB, site :: Int64, prefactor :: Real) where {
        N, SPSSBS<:BasisStateXYZ, SPBS <: SPMSBasisState{SPSSBS},
        MPB <: MPBasis{N,SPBS}
    }

    # select all single particle states that are on site site
    sp_state_indices = Int64[i for i in 1:length(basis.single_particle_basis) if basis.single_particle_basis[i].site == site]
    sp_states = SPBS[basis.single_particle_basis[i] for i in sp_state_indices]

    # generate a new operator
    op = MPElectronDensityDensityOperator(basis, Tuple{Int64,Int64}[], prefactor)

    # try generating all pairs (a,b) s.t. orb(a) != orb(b) but spin(a) = spin(b)
    for a in 1:length(sp_states)
    for b in 1:length(sp_states)
        # check if it is orb(a) != orb(b) but spin(a) = spin(b)
        if sp_states[a].state.orbital != sp_states[b].state.orbital && sp_states[a].state.ms == sp_states[b].state.ms
            push!(op.interacting_orbitals, (sp_state_indices[a],sp_state_indices[b]))
        end
    end
    end

    # recalculate the operator
    recalculate!(op)

    # return the finished operator
    return op
end
export generateDensityDensityElectronInteractionSameSpin

# generate a term of the form n_aup n_apdown  for site i (ELECTRON)
function generateDensityDensityElectronInteractionRemainder(basis :: MPB, site :: Int64, prefactor :: Real) where {
        N, SPSSBS<:BasisStateXYZ, SPBS <: SPMSBasisState{SPSSBS},
        MPB <: MPBasis{N,SPBS}
    }

    # select all single particle states that are on site site
    sp_state_indices = Int64[i for i in 1:length(basis.single_particle_basis) if basis.single_particle_basis[i].site == site]
    sp_states = SPBS[basis.single_particle_basis[i] for i in sp_state_indices]

    # generate a new operator
    op = MPElectronDensityDensityOperator(basis, Tuple{Int64,Int64}[], prefactor)

    # try generating all pairs (a,b) s.t. orb(a) != orb(b) but spin(a) != spin(b)
    for a in 1:length(sp_states)
    for b in 1:length(sp_states)
        # check if it is orb(a) != orb(b) but spin(a) != spin(b)
        if sp_states[a].state.orbital != sp_states[b].state.orbital && sp_states[a].state.ms > sp_states[b].state.ms
            push!(op.interacting_orbitals, (sp_state_indices[a],sp_state_indices[b]))
        end
    end
    end

    # recalculate the operator
    recalculate!(op)

    # return the finished operator
    return op
end
export generateDensityDensityElectronInteractionRemainder



# generate a term of the form n_aup n_adown for site i (HOLE)
function generateDensityDensityHoleInteractionSameOrbital(basis :: MPB, site :: Int64, prefactor :: Real) where {
        N, SPSSBS<:BasisStateXYZ, SPBS <: SPMSBasisState{SPSSBS},
        MPB <: MPBasis{N,SPBS}
    }

    # select all single particle states that are on site site
    sp_state_indices = Int64[i for i in 1:length(basis.single_particle_basis) if basis.single_particle_basis[i].site == site]
    sp_states = SPBS[basis.single_particle_basis[i] for i in sp_state_indices]

    # generate a new operator
    op = MPHoleDensityDensityOperator(basis, Tuple{Int64,Int64}[], prefactor)

    # try generating all pairs (a,b) s.t. orb(a) = orb(b) but spin(a) != spin(b)
    for a in 1:length(sp_states)
    for b in 1:length(sp_states)
        # check if it is orb(a) = orb(b) but spin(a) != spin(b)
        if sp_states[a].state.orbital == sp_states[b].state.orbital && sp_states[a].state.ms > sp_states[b].state.ms
            push!(op.interacting_orbitals, (sp_state_indices[a],sp_state_indices[b]))
        end
    end
    end

    # recalculate the operator
    recalculate!(op)

    # return the finished operator
    return op
end
export generateDensityDensityHoleInteractionSameOrbital

# generate a term of the form n_aup n_apup + n_adown n_apdown for site i (HOLE)
function generateDensityDensityHoleInteractionSameSpin(basis :: MPB, site :: Int64, prefactor :: Real) where {
        N, SPSSBS<:BasisStateXYZ, SPBS <: SPMSBasisState{SPSSBS},
        MPB <: MPBasis{N,SPBS}
    }

    # select all single particle states that are on site site
    sp_state_indices = Int64[i for i in 1:length(basis.single_particle_basis) if basis.single_particle_basis[i].site == site]
    sp_states = SPBS[basis.single_particle_basis[i] for i in sp_state_indices]

    # generate a new operator
    op = MPHoleDensityDensityOperator(basis, Tuple{Int64,Int64}[], prefactor)

    # try generating all pairs (a,b) s.t. orb(a) != orb(b) but spin(a) = spin(b)
    for a in 1:length(sp_states)
    for b in 1:length(sp_states)
        # check if it is orb(a) != orb(b) but spin(a) = spin(b)
        if sp_states[a].state.orbital != sp_states[b].state.orbital && sp_states[a].state.ms == sp_states[b].state.ms
            push!(op.interacting_orbitals, (sp_state_indices[a],sp_state_indices[b]))
        end
    end
    end

    # recalculate the operator
    recalculate!(op)

    # return the finished operator
    return op
end
export generateDensityDensityHoleInteractionSameSpin

# generate a term of the form n_aup n_apdown  for site i (HOLE)
function generateDensityDensityHoleInteractionRemainder(basis :: MPB, site :: Int64, prefactor :: Real) where {
        N, SPSSBS<:BasisStateXYZ, SPBS <: SPMSBasisState{SPSSBS},
        MPB <: MPBasis{N,SPBS}
    }

    # select all single particle states that are on site site
    sp_state_indices = Int64[i for i in 1:length(basis.single_particle_basis) if basis.single_particle_basis[i].site == site]
    sp_states = SPBS[basis.single_particle_basis[i] for i in sp_state_indices]

    # generate a new operator
    op = MPHoleDensityDensityOperator(basis, Tuple{Int64,Int64}[], prefactor)

    # try generating all pairs (a,b) s.t. orb(a) != orb(b) but spin(a) != spin(b)
    for a in 1:length(sp_states)
    for b in 1:length(sp_states)
        # check if it is orb(a) != orb(b) but spin(a) != spin(b)
        if sp_states[a].state.orbital != sp_states[b].state.orbital && sp_states[a].state.ms > sp_states[b].state.ms
            push!(op.interacting_orbitals, (sp_state_indices[a],sp_state_indices[b]))
        end
    end
    end

    # recalculate the operator
    recalculate!(op)

    # return the finished operator
    return op
end
export generateDensityDensityHoleInteractionRemainder
