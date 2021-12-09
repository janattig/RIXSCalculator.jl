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

function find_bsstr(io::IO, bs::BasisStateLS)
    return haskey(io, :compact) ? "" : "LS"
end
function find_bsstr(io::IO, bs::BasisStateXYZ)
    return haskey(io, :compact) ? "" : "XYZ"
end
function find_bsstr(io::IO, bs::BasisStateJ)
    return haskey(io, :compact) ? "" : "J"
end
function find_bsstr(io::IO, bs::BasisStateA1G)
    return haskey(io, :compact) ? "" : "A1G"
end

import Base.show
function Base.show(io::IO, state::DelocalizedBasisState{BS} where BS) 
    
    bsstr = find_bsstr(io, state.state)
    print(io, bsstr*summary(state.state, ["|#1,","⟩ "])*( state.bonding_type == :bonding ? '+' : '-' )*summary(state.state, ["|#2,","⟩ "])
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
    getT2GDelocalizedBasisXYZ() :: SPBasis{DelocalizedBasisState{SPSSBS} where SPSSBS <: AbstractSPSSBasisState} 

This function provides the pre-implemented single particle - two site bonding and antibonding XYZ basis for the t2g.
"""
function getT2GDelocalizedBasisXYZ() :: SPBasis{DelocalizedBasisState{SPSSBS} where SPSSBS <: AbstractSPSSBasisState} 
    # reset the basis object in the RIXSSite object
    XYZBasis=getT2GBasisXYZ()
    states = DelocalizedBasisState[
        # add all terms into the basis
        DelocalizedBasisState(XYZBasis[1],:bonding,1,2),
        DelocalizedBasisState(XYZBasis[1],:antibonding,1,2),
        DelocalizedBasisState(XYZBasis[2],:bonding,1,2),
        DelocalizedBasisState(XYZBasis[2],:antibonding,1,2),
        DelocalizedBasisState(XYZBasis[3],:bonding,1,2),
        DelocalizedBasisState(XYZBasis[3],:antibonding,1,2),
        DelocalizedBasisState(XYZBasis[4],:bonding,1,2),
        DelocalizedBasisState(XYZBasis[4],:antibonding,1,2),
        DelocalizedBasisState(XYZBasis[5],:bonding,1,2),
        DelocalizedBasisState(XYZBasis[5],:antibonding,1,2),
        DelocalizedBasisState(XYZBasis[6],:bonding,1,2),
        DelocalizedBasisState(XYZBasis[6],:antibonding,1,2),
    ]
    # return the basis
    return SPBasis{DelocalizedBasisState}(states)
end

export getT2GDelocalizedBasisXYZ