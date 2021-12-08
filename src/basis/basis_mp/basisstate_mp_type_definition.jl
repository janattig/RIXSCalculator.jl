#############################################
# Definition of MULTI PARTICLE Basis states #
#############################################

# (defined on any number of particles N)
"""
    MPBasisState{N} <: AbstractBasisState

Mutable struct defining the multi particle basis state on any number of particles `N`.
"""
mutable struct MPBasisState{N} <: AbstractBasisState
    # only characterized by the occupation of different single particle states
    occupation :: Vector{Int64}

    # information on basis: index and sign
    basis_index :: Int64
    basis_sign  :: Int64

    # Custom constructor: only use occupation
    function MPBasisState{N}(occupation :: Vector{<:Integer}) where {N}
        return new{N}(occupation, -1, 0)
    end
end

# export
export MPBasisState


# custom showing of MP basis state (without basis reference) ⊗
import Base.show
function Base.show(io::IO, bs::MPBasisState{N}) where {N}
    if haskey(io, :compact)
        if bs.basis_index < 0
            bs.basis_sign = permutation_sign(bs.occupation)
        end
        operator_string = bs.basis_sign<0 ? "-|" : (bs.basis_sign==0 ? "0|" : "+|")
        for o in bs.occupation
            operator_string *= "($(o))"
        end
        operator_string *= "⟩"
        print(io, operator_string)
    else
        # calculate the sign of the permutation
        if bs.basis_index < 0
            bs.basis_sign = permutation_sign(bs.occupation)
        end
        operator_string = bs.basis_sign<0 ? "(-)" : (bs.basis_sign==0 ? "(0)" : "(+)")
        for o in bs.occupation
            operator_string *= " a^d_{$(o)}"
        end
        print(io, operator_string * " |vac⟩" * (bs.basis_index > 0 ? "   (i=$(bs.basis_index))" : ""))
    end
end




# general overlap function as defined in abstract way
function overlap(state_1 :: MPBasisState{N1}, state_2 :: MPBasisState{N2}) :: Complex{Float64} where {N1,N2}
    @error "overlap cannot be called for multi particle states without basis!" stacktrace()
    return NaN+ NaN*im
end
