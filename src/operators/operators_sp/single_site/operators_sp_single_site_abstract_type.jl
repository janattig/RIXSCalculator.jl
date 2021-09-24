# abstract operator super type to all SP SS operators
abstract type AbstractSPSSOperator{
        SPB <: SPBasis{SPSSBS} where {SPSSBS <: AbstractSPSSBasisState}
} <: AbstractSPOperator{SPB} end
export AbstractSPSSOperator





import Base.show
function Base.show(io::IO, op::OP) where {OP <: AbstractSPSSOperator}
    print(io, "Unknown SP SS operator of type $(OP), acting on the following basis:\n")
    display(TextDisplay(io), basis(op))
end
