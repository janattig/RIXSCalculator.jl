struct DelocalizedBasisState{B<:AbstractSPSSBasisState} <: AbstractSPBasisState
    # the SPSS basis state itself
    state :: B
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
function Base.show(io::IO, state::DelocalizedBasisState{BS} where BS) 
    
    # note: state.state must be of a type whose name is "BasisState"*s, where s is a string of arbitrary legth
    
    bsstr =  haskey(io, :compact) ? "" : string(typeof(state.state))[11:length(string(typeof(state.state)))]
    
    print(io, bsstr*summary(state.state, ["|#1,","⟩ "])*( state.bonding_type == :bonding ? '+' : '-' )*" "*bsstr*summary(state.state, ["|#2,","⟩ "])
)
end

# custom summary function
function summary(bs::DelocalizedBasisState, brackets="()")
    return brackets[1]*"$(bs.bonding_type),$(bs.state),$(bs.state.ms>0 ? '↑' : '↓')"*brackets[2]
end


# export the basis type
export DelocalizedBasisState



#########################################
#   Pre-implemented delocalized basis   #
#########################################


"""
    getDelocalizedBasis(basis :: SPBasis{BS},site1::Int64, site2::Int64) :: SPBasis{DelocalizedBasisState{BS}} where {BS<:AbstractSPSSBasisState}

This function provides the pre-implemented single particle - two site bonding and antibonding XYZ basis for the t2g.
"""
function getDelocalizedBasis(basis :: SPBasis{BS},site1::Int64, site2::Int64) :: SPBasis{DelocalizedBasisState{BS}} where {BS<:AbstractSPSSBasisState}
    # make a list of multisite states
    delocalized_states = DelocalizedBasisState{BS}[]
    # push all states for all sites
    for b in states(basis)
        push!(delocalized_states, DelocalizedBasisState{BS}(b, :bonding, site1, site2))
        push!(delocalized_states, DelocalizedBasisState{BS}(b, :antibonding, site1, site2))
    end
    # return the multisite basis
    return SPBasis{DelocalizedBasisState{BS}}(delocalized_states)
end

export getDelocalizedBasis