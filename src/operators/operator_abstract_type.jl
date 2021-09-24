# abstract operator super type to all operators
abstract type AbstractOperator{
        B <: AbstractBasis{BS} where {BS <: AbstractBasisState}
} end
export AbstractOperator




# overwrite print function for custom types
import Base.show
function Base.show(io::IO, op::OP) where {OP <: AbstractOperator}
    print(io, "Unknown operator of type $(OP), acting on the following basis:\n")
    display(TextDisplay(io), basis(op))
end




################################################################################
#
#   DEFINITION OF INTERFACE FUNCTIONS
#
################################################################################


# obtain the current basis
function basis(operator :: OP) where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, OP<:AbstractOperator{B}}
    @error "Interface function 'basis' not implemented for operator of type $(OP)" stacktrace()
end
export basis

# obtain the matrix representation
function matrix_representation(operator :: OP) where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, OP<:AbstractOperator{B}}
    @error "Interface function 'matrix_representation' not implemented for operator of type $(OP)" stacktrace()
end
export matrix_representation

# obtain the matrix element <state_1 | op | state_2>
function matrix_element(operator :: OP, state_1 :: Vector{<:Number}, state_2 :: Vector{<:Number}) :: Complex where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, OP<:AbstractOperator{B}}
    return dot(state_1, matrix_representation(operator) * state_2)
end
export matrix_element

# possibly recalculate the matrix representation
function recalculate!(operator :: OP, recursive::Bool=true, basis_change::Bool=true) where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, OP<:AbstractOperator{B}}
    @error "Interface function 'recalculate!' not implemented for operator of type $(OP)" stacktrace()
end
export recalculate!




# set a parameter (returns (found parameter?, changed matrix?))
function set_parameter!(operator :: OP, parameter :: Symbol, value; print_result::Bool=false, recalculate::Bool=true, kwargs...) where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, OP<:AbstractOperator{B}}
    @error "Interface function 'set_parameter' not implemented for operator of type $(OP)" stacktrace()
end
export set_parameter!

# get a parameter (returns parameter value or nothing)
function get_parameter(operator :: OP, parameter :: Symbol; print_result::Bool=false, kwargs...) where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, OP<:AbstractOperator{B}}
    @error "Interface function 'get_parameter' not implemented for operator of type $(OP)" stacktrace()
end
export get_parameter

# get a list of parameters
function get_parameters(operator :: OP; kwargs...) where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, OP<:AbstractOperator{B}}
    @error "Interface function 'get_parameters' not implemented for operator of type $(OP)" stacktrace()
end
export get_parameters
