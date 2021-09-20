####################################################
# Single Particle Single Site Composite Basis State #
####################################################

# basis state build on superposition
# |state> = prefactor_i |basisstate_i>
struct SPSSCompositeBasisState{SPB <: SPBasis{SPSS} where SPSS<:AbstractSPSSBasisState} <: AbstractSPSSBasisState
    # the prefactors of basis states
    prefactors :: Vector{Complex{Float64}}
    # the original basis
    basis      :: SPB
end
export SPSSCompositeBasisState

function SPSSCompositeBasisState(prefactors :: Vector{<:Number}, basis :: B) where {B <: SPBasis{SPSS} where SPSS<:AbstractSPSSBasisState}
    return SPSSCompositeBasisState{B}(prefactors, basis)
end
function CompositeBasisState(prefactors :: Vector{<:Number}, basis :: B) where {B <: SPBasis{SPSS} where SPSS<:AbstractSPSSBasisState}
    return SPSSCompositeBasisState(prefactors, basis)
end
export SPSSCompositeBasisState, CompositeBasisState