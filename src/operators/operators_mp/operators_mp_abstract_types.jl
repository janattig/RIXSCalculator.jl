# abstract type definition of single particle operators
abstract type AbstractMP1POperator{
        MPB <: MPBasis{N,SPBS} where {N,SPBS <: AbstractSPBasisState}
} <: AbstractMPOperator{1,MPB} end


# abstract type definition of single particle operators
abstract type AbstractMP2POperator{
        MPB <: MPBasis{N,SPBS} where {N,SPBS <: AbstractSPBasisState}
} <: AbstractMPOperator{2,MPB} end


# abstract type definition (up to NI particles interact)
abstract type AbstractMPInteractionHamiltonian{
    NI,
    MPB <: MPBasis{N,SPBS} where {N,SPBS <: AbstractSPBasisState}
} <: AbstractMPOperator{NI,MPB} end
