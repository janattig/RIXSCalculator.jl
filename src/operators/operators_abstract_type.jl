################################################################################
#
#   DEFINITION OF DIFFERENT ABSTRACT OPERATOR TYPES
#   TO LAY OUT THE TYPE HIERACHY
#
################################################################################

# abstract operator super type to all operators
abstract type AbstractOperator{
        B <: AbstractBasis{BS} where {BS <: AbstractBasisState}
} end

# abstract operator super type to all SP operatorss
abstract type AbstractSPOperator{
        SPB <: SPBasis{SPBS} where {SPBS <: AbstractSPBasisState}
} <: AbstractOperator{SPB} end
# abstract operator super type to all MP operators (NI intracting particles)
abstract type AbstractMPOperator{
        NI,
        MPB <: MPBasis{N,SPBS} where {N,SPBS <: AbstractSPBasisState}
} <: AbstractOperator{MPB} end

# abstract operator super type to all SP SS operators
abstract type AbstractSPSSOperator{
        SPB <: SPBasis{SPSSBS} where {SPSSBS <: AbstractSPSSBasisState}
} <: AbstractSPOperator{SPB} end
# abstract operator super type to all SP MS operators
abstract type AbstractSPMSOperator{
        SPB <: SPBasis{SPMSBS} where {SPMSBS<:Union{SPMSBasisState{SPSSBS} where SPSSBS, SPMSCompositeBasisState{B} where B}}
} <: AbstractSPOperator{SPB} end


# export all operator types
export  AbstractOperator,
        AbstractSPOperator,
        AbstractMPOperator,
        AbstractSPSSOperator,
        AbstractSPMSOperator


# overwrite print function for custom types
import Base.show


# Base.show functions (FALLBACKS)
function Base.show(io::IO, op::OP) where {OP <: AbstractOperator}
    print(io, "Unknown operator of type $(OP), acting on the following basis:\n")
    display(TextDisplay(io), basis(op))
end
function Base.show(io::IO, op::OP) where {OP <: AbstractSPOperator}
    print(io, "Unknown SP operator of type $(OP), acting on the following basis:\n")
    display(TextDisplay(io), basis(op))
end
function Base.show(io::IO, op::OP) where {OP <: AbstractSPSSOperator}
    print(io, "Unknown SP SS operator of type $(OP), acting on the following basis:\n")
    display(TextDisplay(io), basis(op))
end
function Base.show(io::IO, op::OP) where {OP <: AbstractSPMSOperator}
    print(io, "Unknown SP MS operator of type $(OP), acting on the following basis:\n")
    display(TextDisplay(io), basis(op))
end
function Base.show(io::IO, op::OP) where {OP <: AbstractMPOperator}
    print(io, "Unknown multi-particle operator of type $(OP), acting on the following basis:\n")
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

# obtain the matrix representation
function matrix_representation(operator :: OP) where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, OP<:AbstractOperator{B}}
    @error "Interface function 'matrix_representation' not implemented for operator of type $(OP)" stacktrace()
end

# obtain the matrix element <state_1 | op | state_2>
function matrix_element(operator :: OP, state_1 :: Vector{<:Number}, state_2 :: Vector{<:Number}) :: Complex where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, OP<:AbstractOperator{B}}
    return dot(state_1, matrix_representation(operator) * state_2)
end

# possibly recalculate the matrix representation
function recalculate!(operator :: OP, recursive::Bool=true, basis_change::Bool=true) where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, OP<:AbstractOperator{B}}
    @error "Interface function 'recalculate!' not implemented for operator of type $(OP)" stacktrace()
end




# set a parameter (returns (found parameter?, changed matrix?))
function set_parameter!(operator :: OP, parameter :: Symbol, value; print_result::Bool=false, recalculate::Bool=true, kwargs...) where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, OP<:AbstractOperator{B}}
    @error "Interface function 'set_parameter' not implemented for operator of type $(OP)" stacktrace()
end

# get a parameter (returns parameter value or nothing)
function get_parameter(operator :: OP, parameter :: Symbol; print_result::Bool=false, kwargs...) where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, OP<:AbstractOperator{B}}
    @error "Interface function 'get_parameter' not implemented for operator of type $(OP)" stacktrace()
end

# get a list of parameters
function get_parameters(operator :: OP; kwargs...) where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, OP<:AbstractOperator{B}}
    @error "Interface function 'get_parameters' not implemented for operator of type $(OP)" stacktrace()
end



# export the inteface functions
export  basis,
        matrix_representation,
        matrix_element,
        recalculate!,
        set_parameter!,
        get_parameter