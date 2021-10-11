# Abstract super type for all bases
"""
    AbstractBasis{BS <: AbstractBasisState} <: AbstractArray{BS,1}

Abstract supertype for all bases. 
"""
abstract type AbstractBasis{BS <: AbstractBasisState} <: AbstractArray{BS,1} end
export  AbstractBasis




# Interface function for accessing states that has to be overwritten by all bases
"""
    states(basis::B) where {BS, B<:AbstractBasis{BS}}

Interface function for accessing states.
"""
function states(basis::B) where {BS, B<:AbstractBasis{BS}}
    error("function 'states' not implemented for basis of type "*string(B))
end
export states

# Interface function for obtaining a list of sites on which the states are defined
"""
    get_sites(basis::B) where {BS, B<:AbstractBasis{BS}}

Interface function for obtaining a list of sites on which the states are defined.
"""
function get_sites(basis::B) where {BS, B<:AbstractBasis{BS}}
    error("function 'get_sites' not implemented for basis of type "*string(B))
end
export get_sites







# Indexing
import Base.IndexStyle
function Base.IndexStyle(::Type{<:AbstractBasis})
    return IndexLinear()
end

# size
import Base.size
function Base.size(basis::B) where {BS, B<:AbstractBasis{BS}}
    return size(states(basis))
end

# getindex
import Base.getindex
function Base.getindex(basis::B, i::Int64) where {BS, B<:AbstractBasis{BS}}
    return getindex(states(basis), i)
end

# setindex
import Base.setindex!
function Base.setindex!(basis::B, v, i::Int64) where {BS, B<:AbstractBasis{BS}}
    return setindex!(states(basis), v, i)
end
