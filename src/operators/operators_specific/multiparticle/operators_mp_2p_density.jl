##############################################################
#
#   MP operators for density density interactions
#   - type definition
#   - interface functions
#   - convenience functions
#
##############################################################





##############################################################
#   Type definition
##############################################################

# define some abstract density density operator
abstract type AbstractMPDensityDensityOperator{MPB} <: AbstractMP2POperator{MPB} end

# define a MP operator for density density interactions of ELECTRONS (n*n)
mutable struct MPElectronDensityDensityOperator{
        MPB <: MPBasis{N,SPBS} where {N,SPBS}
    } <: AbstractMPDensityDensityOperator{MPB}

    # the MP basis
    basis :: MPB
    # the matrix representation
    matrix_rep :: Matrix{Complex{Float64}}
    # the list of interacting orbitals
    interacting_orbitals :: Vector{Tuple{Int64,Int64}}
    # the prefactor
    prefactor :: Float64

    # custom constructor
    function MPElectronDensityDensityOperator(basis :: MPB, interacting_orbitals :: Vector{Tuple{Int64,Int64}}, prefactor :: Real) where {
                N, SPBS <: AbstractSPBasisState,
                MPB <: MPBasis{N,SPBS}
            }
        # create a new operator
        op = new{MPB}(basis, zeros(length(basis), length(basis)), interacting_orbitals, prefactor)
        # recalculate the matrix representation
        recalculate!(op)
        # return the operator
        return op
    end
end
export MPElectronDensityDensityOperator

# define a MP operator for density density interactions of HOLES ((1-n)*(1-n))
mutable struct MPHoleDensityDensityOperator{
        MPB <: MPBasis{N,SPBS} where {N,SPBS}
    } <: AbstractMPDensityDensityOperator{MPB}

    # the MP basis
    basis :: MPB
    # the matrix representation
    matrix_rep :: Matrix{Complex{Float64}}
    # the list of interacting orbitals
    interacting_orbitals :: Vector{Tuple{Int64,Int64}}
    # the prefactor
    prefactor :: Float64

    # custom constructor
    function MPHoleDensityDensityOperator(basis :: MPB, interacting_orbitals :: Vector{Tuple{Int64,Int64}}, prefactor :: Real) where {
                N, SPBS <: AbstractSPBasisState,
                MPB <: MPBasis{N,SPBS}
            }
        # create a new operator
        op = new{MPB}(basis, zeros(length(basis), length(basis)), interacting_orbitals, prefactor)
        # recalculate the matrix representation
        recalculate!(op)
        # return the operator
        return op
    end
end
export MPHoleDensityDensityOperator


import Base.show
function Base.show(io::IO, op::MPHoleDensityDensityOperator{MPB}) where {
            N,
            SPBS <: AbstractSPBasisState,
            MPB <: MPBasis{N,SPBS}
        }
    if haskey(io, :compact)
        print(io, "2-particle (hole) density-density interaction with "*string(length(op.interacting_orbitals))*" n*n terms")
    else
        print(io, "2-particle (hole) density-density interaction\n")
        print(io, "Interactions include "*string(length(op.interacting_orbitals))*" terms:\n")
        for orbital_pair in op.interacting_orbitals
            orb_1 = summary(basis(op).single_particle_basis[orbital_pair[1]], "{}")
            orb_2 = summary(basis(op).single_particle_basis[orbital_pair[2]], "{}")
            print(io, " + n_"*orb_1*"*n_"*orb_2*"\n")
        end
        print(io, "Overall prefactor is "*string(op.prefactor)*"\n")
        print(io, "Multi-particle basis contains "*string(length(basis(op)))*" states in total, with "*string(N)*" particles per state\n")
    end
end
function Base.show(io::IO, op::MPElectronDensityDensityOperator{MPB}) where {
            N,
            SPBS <: AbstractSPBasisState,
            MPB <: MPBasis{N,SPBS}
        }
    if haskey(io, :compact)
        print(io, "2-particle (electron) density-density interaction with "*string(length(op.interacting_orbitals))*" n*n terms")
    else
        print(io, "2-particle (electron) density-density interaction\n")
        print(io, "Interactions include "*string(length(op.interacting_orbitals))*" terms:\n")
        for orbital_pair in op.interacting_orbitals
            orb_1 = summary(basis(op).single_particle_basis[orbital_pair[1]], "{}")
            orb_2 = summary(basis(op).single_particle_basis[orbital_pair[2]], "{}")
            print(io, " + n_"*orb_1*"*n_"*orb_2*"\n")
        end
        print(io, "Overall prefactor is "*string(op.prefactor)*"\n")
        print(io, "Multi-particle basis contains "*string(length(basis(op)))*" states in total, with "*string(N)*" particles per state\n")
    end
end




##############################################################
#   Interface functions
##############################################################

# obtain the current basis (ELECTRON & HOLE)
function basis(operator :: MPDDOP) :: MPB where {
            N, SPBS <: AbstractSPBasisState,
            MPB <: MPBasis{N,SPBS},
            MPDDOP <: AbstractMPDensityDensityOperator{MPB}
        }
    return operator.basis
end

# obtain the matrix representation (ELECTRON & HOLE)
function matrix_representation(operator :: MPDDOP) :: Matrix{Complex{Float64}} where {
            N, SPBS <: AbstractSPBasisState,
            MPB <: MPBasis{N,SPBS},
            MPDDOP <: AbstractMPDensityDensityOperator{MPB}
        }
    return operator.matrix_rep .* operator.prefactor
end


# possibly recalculate the matrix representation (ELECTRON & HOLE) (Fallback for non XYZ)
function recalculate!(operator :: MPDDOP, recursive::Bool=true, basis_change::Bool=true) where {
            N, SPBS <: AbstractSPBasisState,
            MPB <: MPBasis{N,SPBS},
            MPDDOP <: AbstractMPDensityDensityOperator{MPB}
        }
    @error "currently only recalculation of density-density operator implemented for XYZ basis states, not basis states of type $(SPBS)" stacktrace()
end

# possibly recalculate the matrix representation (ELECTRON & HOLE)
function recalculate!(operator :: MPDDOP, recursive::Bool=true, basis_change::Bool=true) where {
            N, SPBS <: Union{SPMSBasisState{BasisStateXYZ}, BasisStateXYZ},
            MPB <: MPBasis{N,SPBS},
            MPDDOP <: AbstractMPDensityDensityOperator{MPB}
        }
    # create new matrix
    operator.matrix_rep = zeros(Complex{Float64}, length(basis(operator)), length(basis(operator)))
    # generate all diagonal element contributions to the many body hamiltonian
    for alpha in 1:length(basis(operator))
        # add the expectation with ab to the matrix
        operator.matrix_rep[alpha, alpha] = orbital_contribution(operator, basis(operator)[alpha])
    end
end


# contribution of one orbital (ELEKTRON)
function orbital_contribution(operator :: MPElectronDensityDensityOperator{MPB}, state :: MPBasisState{N}) :: Complex{Float64} where {
            N, SPBS <: AbstractSPBasisState,
            MPB <: MPBasis{N,SPBS}
        }
    # interaction energy
    interaction = 0.0
    # check all orbital interactions
    for orbital in operator.interacting_orbitals
        # check if they are both occupied in the current state
        n_a = orbital[1] in state.occupation ? 1 : 0
        n_b = orbital[2] in state.occupation ? 1 : 0
        interaction += n_a * n_b
    end
    # return the interaction contribution
    return interaction
end

# contribution of one orbital (HOLE)
function orbital_contribution(operator :: MPHoleDensityDensityOperator{MPB}, state :: MPBasisState{N}) :: Complex{Float64} where {
            N, SPBS <: AbstractSPBasisState,
            MPB <: MPBasis{N,SPBS}
        }
    # interaction energy
    interaction = 0.0
    # check all orbital interactions
    for orbital in operator.interacting_orbitals
        # check if they are both occupied in the current state
        n_a = orbital[1] in state.occupation ? 1 : 0
        n_b = orbital[2] in state.occupation ? 1 : 0
        interaction += (1-n_a) * (1-n_b)
    end
    # return the interaction contribution
    return interaction
end

export orbital_contribution

##############################################################
#   Convenience functions for generating suitable terms
##############################################################

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
