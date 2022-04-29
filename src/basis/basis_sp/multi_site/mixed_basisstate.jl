struct MixedBasisState{B1<:AbstractSPSSBasisState,B2<:AbstractSPSSBasisState} <: AbstractSPBasisState
    # the SPSS basis state itself
    state1 :: B1
    state2 :: B2
    # bonding type
    bonding_type :: Symbol
    # coefficients
#     coeff :: Vector{Complex{Float64}}
    # sites
    site1 :: Int64
    site2 :: Int64
end

# custom print function

import Base.show
function Base.show(io::IO, state::MixedBasisState{BS} where BS) 
    
    # note: state.state must be of a type whose name is "BasisState"*s, where s is a string of arbitrary legth
    
    bsstr1 =  haskey(io, :compact) ? "" : string(typeof(state.state1))[11:length(string(typeof(state.state1)))]
    bsstr2 =  haskey(io, :compact) ? "" : string(typeof(state.state2))[11:length(string(typeof(state.state2)))]
    
    # print(io, bsstr1*summary(state.state1, ["|#"*string(state.site1)*",","⟩ "])*( state.bonding_type == :b ? '+' : '-' )*" "*bsstr2*summary(state.state2, ["|#"*string(state.site2)*",","⟩ "])
    print(io, bsstr1*bsstr2*summary(state.state1, ["|"*(state.bonding_type == :b ? "b " : "ab")*",","⟩ "])
)
end

# custom summary function
function summary(bs::MixedBasisState, brackets="()")
    # state1 and state2 are supposed to have the same quantum numbers
    return brackets[1]*(bs.bonding_type == :b ? "b " : "ab")*summary(bs.state1,[',', brackets[2]])
end


# export the basis type
export MixedBasisState


"""
    getMixedBasis(basis1 :: SPBasis{BS1},basis2 :: SPBasis{BS2},site1::Int64, site2::Int64) :: SPBasis{MixedBasisState{BS1,BS2}} where {BS1<:AbstractSPSSBasisState, BS2<:AbstractSPSSBasisState}

This function provides the pre-implemented single particle - two site bonding and antibonding basis for the face-sharing configuration. Note that this basis has been tested
only for the BasisStateJ1 and BasisStateJ2 types.
"""
function getMixedBasis(basis1 :: SPBasis{BS1},basis2 :: SPBasis{BS2},site1::Int64, site2::Int64) :: SPBasis{MixedBasisState{BS1,BS2}} where {BS1<:AbstractSPSSBasisState, BS2<:AbstractSPSSBasisState}
    # make a list of multisite states
    mixed_states = Union{MixedBasisState{BS1}, MixedBasisState{BS2}}[]
    # push all states for all sites
    for i in 1:length(basis1)
        push!(mixed_states, MixedBasisState{BS1,BS2}(basis1[i], basis2[i], :b, site1, site2))
        push!(mixed_states, MixedBasisState{BS1,BS2}(basis1[i], basis2[i], :ab, site1, site2))
    end
    # return the multisite basis
    return SPBasis{MixedBasisState{BS1,BS2}}(mixed_states)
end

export getMixedBasis