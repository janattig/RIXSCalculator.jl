# abstract operator super type to all SP operatorss
abstract type AbstractSPOperator{
        SPB <: SPBasis{SPBS} where {SPBS <: AbstractSPBasisState}
} <: AbstractOperator{SPB} end
export AbstractSPOperator




# overwrite print function for custom types
import Base.show
function Base.show(io::IO, op::OP) where {OP <: AbstractSPOperator}
    print(io, "Unknown SP operator of type $(OP), acting on the following basis:\n")
    display(TextDisplay(io), basis(op))
end
