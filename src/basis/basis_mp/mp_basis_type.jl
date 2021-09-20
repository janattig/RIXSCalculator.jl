##########################
# Concrete MP Basis type #
##########################

# Type definition
mutable struct MPBasis{N,SPBS<:AbstractSPBasisState} <: AbstractBasis{MPBasisState{N}}
    # states of the many particle basis
    states :: Vector{MPBasisState{N}}
    # single particle basis
    single_particle_basis :: SPBasis{SPBS}

    # lookup: occupation(MPBS) --> (sign of permutation, index)
    lookup_index :: Dict{Vector{Int64}, Tuple{Int64,Int64}}
    # lookup: sp state index i --> Vector{state index of states containing sp state i}
    lookup_sp_states :: Vector{Vector{Int64}}
end
export MPBasis

# custom constructors
function MPBasis{N,SPBS}(states :: Vector{MPBasisState{N}}, single_particle_basis :: SPBasis{SPBS}) where {N, SPBS<:AbstractSPBasisState}
    # use default constructor
    basis = MPBasis{N,SPBS}(states, single_particle_basis, Dict{Vector{Int64}, Tuple{Int64,Int64}}(), Vector{Int64}[])
    # recalculate lookups
    resetLookupIndex!(basis)
    resetLookupSPStates!(basis)
    # return basis
    return basis
end


# custom showing of MP basis state (wtth basis reference)
import Base.summary
function Base.summary(io::IO, basis::MPBasis{N,SPBS}) where {N,SPBS}
    print(io, string(length(basis))*"-element $(N)-particle basis.\n")
    print(io, "SP basis is ")
    display(TextDisplay(io),basis.single_particle_basis)
    print(io, "\nmulti-particle states are")
end



# interface funtion for accessing states
function states(basis :: MPBasis{N,SPBS}) where {N,SPBS<:AbstractSPBasisState}
    return basis.states
end
function get_sites(basis :: MPBasis{N,SPBS}) where {N,SPBS<:AbstractSPBasisState}
    return get_sites(basis.single_particle_basis)
end
export states, get_sites