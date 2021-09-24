# abstract operator super type to all SP MS operators
abstract type AbstractSPMSOperator{
        SPB <: SPBasis{SPMSBS} where {SPMSBS<:Union{SPMSBasisState{SPSSBS} where SPSSBS, SPMSCompositeBasisState{B} where B}}
} <: AbstractSPOperator{SPB} end
export AbstractSPMSOperator




import Base.show
function Base.show(io::IO, op::OP) where {OP <: AbstractSPMSOperator}
    print(io, "Unknown SP MS operator of type $(OP), acting on the following basis:\n")
    display(TextDisplay(io), basis(op))
end
