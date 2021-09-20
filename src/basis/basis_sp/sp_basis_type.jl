# Type definition of SINGLE PARTICLE basis
"""
    SPBasis{BS <: AbstractSPBasisState} <: AbstractBasis{BS}

This object contains a list of single particle basis elements.
"""
mutable struct SPBasis{BS <: AbstractSPBasisState} <: AbstractBasis{BS}
    # only contains a list of basis elements
    states :: Vector{BS}
end
export  SPBasis



# interface function
function states(basis :: SPBasis{BS}) where {BS}
    return basis.states
end
function get_sites(basis::SPBasis{BS}) where {BS <: AbstractSPSSBasisState}
    return Int64[-1, ]
end
export states, get_sites




# custom print function
import Base.summary
function Base.summary(io::IO, basis::SPBasis{BS}) where {BS}
    print(io, string(length(basis))*"-element SP basis for states of type "*string(BS))
end







# index function for finding a basis state in a basis
function index(basis :: SPBasis{SPBS}, state::SPBS) where {SPBS <: AbstractSPBasisState}
    for i in 1:length(basis)
        if state == basis[i]
            return i
        end
    end
    return -1
end
