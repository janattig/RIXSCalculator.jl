# abstract type definition of single particle operators
abstract type AbstractMP1POperator{
        MPB <: MPBasis{N,SPBS} where {N,SPBS <: AbstractSPBasisState}
} <: AbstractMPOperator{1,MPB} end
