# abstract type definition (up to NI particles interact)
abstract type AbstractMPInteractionHamiltonian{
    NI,
    MPB <: MPBasis{N,SPBS} where {N,SPBS <: AbstractSPBasisState}
} <: AbstractMPOperator{NI,MPB} end
