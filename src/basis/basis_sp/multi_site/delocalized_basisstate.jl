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
function Base.show(io::IO, state::DelocalizedBasisState{BS} where {BS <: BasisStateLS}) 
    
    bsstr = haskey(io, :compact) ? "" : "LS"
    
    s=bsstr*"|#1,"*string(state.state.ml)*","*( state.state.ms == 1//2 ? '↑' : '↓' )*"⟩ "
    s=s*( state.bonding_type == :bonding ? '+' : '-' )
    s=s*" "*bsstr*"|#2,"*string(state.state.ml)*","*( state.state.ms == 1//2 ? '↑' : '↓' )*"⟩"
    
    print(io, s)
end

# custom show function for XYZ delocalized basis states
function Base.show(io::IO, state::DelocalizedBasisState{BS} where {BS <: BasisStateXYZ}) 
    
    bsstr = haskey(io, :compact) ? "" : "XYZ"
    
    s=bsstr*"|#1,"*string(state.state.orbital)*","*( state.state.ms == 1//2 ? '↑' : '↓' )*"⟩ "
    s=s*( state.bonding_type == :bonding ? '+' : '-' )
    s=s*" "*bsstr*"|#2,"*string(state.state.orbital)*","*( state.state.ms == 1//2 ? '↑' : '↓' )*"⟩"
    
    print(io, s)
end

# custom show function for J delocalized basis states
function Base.show(io::IO, state::DelocalizedBasisState{BS} where {BS <: BasisStateJ}) 
    
    bsstr = haskey(io, :compact) ? "" : "J"
    
    s=bsstr*"|#1,"*string(state.state.j.num)*"/"*string(state.state.j.den)*","*(state.state.mj==0 ? " " : (state.state.mj>0 ? "+" : "-"))*string(abs(state.state.mj.num))*"/"*string(abs(state.state.mj.den))*"⟩ "
    s=s*( state.bonding_type == :bonding ? '+' : '-' )
    s=s*" "*bsstr*"|#2,"*string(state.state.j.num)*"/"*string(state.state.j.den)*","*(state.state.mj==0 ? " " : (state.state.mj>0 ? "+" : "-"))*string(abs(state.state.mj.num))*"/"*string(abs(state.state.mj.den))*"⟩ "
    
    print(io, s)
end

# custom show function for A1G delocalized basis states
function Base.show(io::IO, state::DelocalizedBasisState{BS} where {BS <: BasisStateA1G}) 
    
    bsstr = haskey(io, :compact) ? "" : "A1G"
    
    s=bsstr*"|#1,"*string(state.state.orbital)*","*( state.state.ms == 1//2 ? '↑' : '↓' )*"⟩ "
    s=s*( state.bonding_type == :bonding ? '+' : '-' )
    s=s*" "*bsstr*"|#2,"*string(state.state.orbital)*","*( state.state.ms == 1//2 ? '↑' : '↓' )*"⟩"
    
    print(io, s)
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