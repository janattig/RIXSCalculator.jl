# abstract operator super type to all MP operators (NI intracting particles)
abstract type AbstractMPOperator{
        NI,
        MPB <: MPBasis{N,SPBS} where {N,SPBS <: AbstractSPBasisState}
} <: AbstractOperator{MPB} end
export AbstractMPOperator




import Base.show
function Base.show(io::IO, op::OP) where {OP <: AbstractMPOperator}
    print(io, "Unknown multi-particle operator of type $(OP), acting on the following basis:\n")
    display(TextDisplay(io), basis(op))
end
