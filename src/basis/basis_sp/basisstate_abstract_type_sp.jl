# Abstract super type for all SINGLE PARTICLE basis states
abstract type AbstractSPBasisState <: AbstractBasisState end
export  AbstractSPBasisState



# general overlap function for any single particle states
function overlap(state_1 :: BS1, state_2 :: BS2) :: Complex{Float64} where {BS1<:AbstractSPBasisState, BS2<:AbstractSPBasisState}
    @error "not implemented function 'overlap' for single particle basis states of types "*string(BS1)*" and "*string(BS2) stacktrace()
    return NaN+ NaN*im
end
export overlap
