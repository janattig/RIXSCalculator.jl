# Abstract super type for all SINGLE PARTICLE basis states
abstract type AbstractSPBasisState <: AbstractBasisState end
export  AbstractSPBasisState



# general overlap function for any single particle states
"""
    overlap(state_1 :: SPSSBS, state_2 :: SPSSBS) :: Complex{Float64}   where {SPSSBS <: AbstractSPSSBasisState}

This function computes the term ``\\left< state_1 | state_2\\right>``.

A complex number is returned.

"""
function overlap(state_1 :: BS1, state_2 :: BS2) :: Complex{Float64} where {BS1<:AbstractSPBasisState, BS2<:AbstractSPBasisState}
    @error "not implemented function 'overlap' for single particle basis states of types "*string(BS1)*" and "*string(BS2) stacktrace()
    return NaN+ NaN*im
end
export overlap
