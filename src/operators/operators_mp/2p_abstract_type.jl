# abstract type definition of single particle operators
abstract type AbstractMP2POperator{
        MPB <: MPBasis{N,SPBS} where {N,SPBS <: AbstractSPBasisState}
} <: AbstractMPOperator{2,MPB} end

export AbstractMP2POperator