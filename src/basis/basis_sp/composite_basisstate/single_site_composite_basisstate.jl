####################################################
# Single Particle Single Site Composite Basis State #
####################################################

# basis state build on superposition
# |state> = prefactor_i |basisstate_i>
"""
    SPMSCompositeBasisState{SPB <: SPBasis{SPMSBasisState{SPSS}} where SPSS<:AbstractSPSSBasisState} <: AbstractSPBasisState

This object defines a single particle multi site composite basis state:

``\\left| state \\right> = \\sum_i prefactor_i \\left| state_i \\right>``

through the prefactors of the basis states `prefactors :: Vector{Complex{Float64}}` and its original basis `basis :: SPB`. 
"""
struct SPSSCompositeBasisState{SPB <: SPBasis{SPSS} where SPSS<:AbstractSPSSBasisState} <: AbstractSPSSBasisState
    # the prefactors of basis states
    prefactors :: Vector{Complex{Float64}}
    # the original basis
    basis      :: SPB
end

function SPSSCompositeBasisState(prefactors :: Vector{<:Number}, basis :: B) where {B <: SPBasis{SPSS} where SPSS<:AbstractSPSSBasisState}
    return SPSSCompositeBasisState{B}(prefactors, basis)
end

export SPSSCompositeBasisState

function CompositeBasisState(prefactors :: Vector{<:Number}, basis :: B) where {B <: SPBasis{SPSS} where SPSS<:AbstractSPSSBasisState}
    return SPSSCompositeBasisState(prefactors, basis)
end
export CompositeBasisState
