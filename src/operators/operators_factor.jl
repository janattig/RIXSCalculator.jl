################################################################################
#
#   Skalar OPERATOR
#
#   - type definition
#   - interface functions
#   - definition of plus
#
################################################################################




################################################################################
#   Type definition
################################################################################


# Type Definition of LS / SpinOrbit Operator
mutable struct ScalarProductOperator{B, O<:AbstractOperator{B}} <: AbstractOperator{B}
    # the basis
    basis :: B
    # the current matrix representation
    matrix_rep :: Matrix{Complex{Float64}}
    # the scalar factor
    factor :: Complex{Float64}
    # the two contained operators
    op :: O

    # Custom constructor (without explicit matrix rep)
    function ScalarProductOperator(factor::Number, op::O) where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, O<:AbstractOperator{B}}
        # construct new operator
        operator = new{B, O}(basis(op), zeros(Complex{Float64}, length(basis(op)), length(basis(op))), factor, op)
        # recalculate the matrix representation
        recalculate!(operator, false)
        # return the operator
        return operator
    end
end

# export operator type
export ScalarProductOperator

import Base.show
function Base.show(io::IO, op::ScalarProductOperator{B, O}) where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, O<:AbstractOperator{B}}
    if haskey(io, :compact)
        print(io, "($(op.factor)) * {")
        show(io, op.op)
        print(io, "}")
    else
        print(io, "$(op.factor) times ")
        show(io, op.op)
    end
end




################################################################################
#   Interface functions
################################################################################

# obtain the current basis
function basis(operator :: ScalarProductOperator{B, O}) :: B where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, O<:AbstractOperator{B}}
    return operator.basis
end

# obtain the matrix representation
function matrix_representation(operator :: ScalarProductOperator{B, O}) :: Matrix{Complex{Float64}} where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, O<:AbstractOperator{B}}
    return operator.matrix_rep
end

# possibly recalculate the matrix representation
function recalculate!(operator :: ScalarProductOperator{B, O}, recursive::Bool=true, basis_change::Bool=true) where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, O<:AbstractOperator{B}}
    # maybe recalculate recursively
    if recursive
        recalculate!(operator.op, true, basis_change)
    end
    # create new matrix
    operator.matrix_rep = matrix_representation(operator.op) .* operator.factor
end

# set a parameter (returns (found parameter?, changed matrix?))
function set_parameter!(operator :: ScalarProductOperator{B, O}, parameter :: Symbol, value; print_result::Bool=false, recalculate::Bool=true, kwargs...) where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, O<:AbstractOperator{B}}
    # check if it can be set in 1
    found_param, changed_matrix = set_parameter!(operator.op, parameter, value, print_result=print_result, recalculate=recalculate; kwargs...)
    if recalculate && changed_matrix
        recalculate!(operator, false, false)
    end
    return (found_param, changed_matrix)
end

# get a parameter (returns (found parameter?, parameter value or nothing))
function get_parameter(operator :: ScalarProductOperator{B, O}, parameter :: Symbol; print_result::Bool=false, kwargs...) where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, O<:AbstractOperator{B}}
    # check if parameter can be obtained
    return get_parameter(operator.op, parameter, print_result=print_result; kwargs...)
end

# get a parameter (returns (found parameter?, parameter value or nothing))
function get_parameters(operator :: ScalarProductOperator{B, O}; kwargs...) where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, O<:AbstractOperator{B}}
    # check if parameter can be obtained
    return get_parameters(operator.op; kwargs...)
end




################################################################################
#   (*) functions
################################################################################

# import of base function
import Base.(*)

# Define the (*) function for any operators
function *(f :: Number, op :: O)  where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, O<:AbstractOperator{B}}
    return ScalarProductOperator(f, op)
end
function *(op :: O, f :: Number)  where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, O<:AbstractOperator{B}}
    return ScalarProductOperator(f, op)
end
























################################################################################
#
#   Scalar OPERATOR with settable prefactor
#
#   - type definition
#   - interface functions
#   - definition of plus
#
################################################################################




################################################################################
#   Type definition
################################################################################


# Type Definition of LS / SpinOrbit Operator
mutable struct SettableScalarProductOperator{B, O<:AbstractOperator{B}} <: AbstractOperator{B}
    # the basis
    basis :: B
    # the current matrix representation
    matrix_rep :: Matrix{Complex{Float64}}
    # the scalar factor
    factor :: Complex{Float64}
    label  :: Symbol
    # the two contained operators
    op :: O

    # Custom constructor (without explicit matrix rep)
    function SettableScalarProductOperator(label::Symbol, factor::Number, op::O) where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, O<:AbstractOperator{B}}
        # construct new operator
        operator = new{B, O}(basis(op), zeros(Complex{Float64}, length(basis(op)), length(basis(op))), factor, label, op)
        # recalculate the matrix representation
        recalculate!(operator, false)
        # return the operator
        return operator
    end
    # Custom constructor (without explicit matrix rep)
    function SettableScalarProductOperator(label::Symbol, op::O) where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, O<:AbstractOperator{B}}
        # construct new operator
        operator = new{B, O}(basis(op), zeros(Complex{Float64}, length(basis(op)), length(basis(op))), 1.0, label, op)
        # recalculate the matrix representation
        recalculate!(operator, false)
        # return the operator
        return operator
    end
end

# export operator type
export SettableScalarProductOperator

import Base.show
function Base.show(io::IO, op::SettableScalarProductOperator{B, O}) where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, O<:AbstractOperator{B}}
    if haskey(io, :compact)
        print(io, "($(op.label)=$(op.factor)) * {")
        show(io, op.op)
        print(io, "}")
    else
        print(io, "($(op.label)=$(op.factor)) times ")
        show(io, op.op)
    end
end




################################################################################
#   Interface functions
################################################################################

# obtain the current basis
function basis(operator :: SettableScalarProductOperator{B, O}) :: B where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, O<:AbstractOperator{B}}
    return operator.basis
end

# obtain the matrix representation
function matrix_representation(operator :: SettableScalarProductOperator{B, O}) :: Matrix{Complex{Float64}} where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, O<:AbstractOperator{B}}
    return operator.matrix_rep
end

# possibly recalculate the matrix representation
function recalculate!(operator :: SettableScalarProductOperator{B, O}, recursive::Bool=true, basis_change::Bool=true) where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, O<:AbstractOperator{B}}
    # maybe recalculate recursively
    if recursive
        recalculate!(operator.op, true, basis_change)
    end
    # create new matrix
    operator.matrix_rep = matrix_representation(operator.op) .* operator.factor
end

# set a parameter (returns (found parameter?, changed matrix?))
function set_parameter!(operator :: SettableScalarProductOperator{B, O}, parameter :: Symbol, value; print_result::Bool=false, recalculate::Bool=true, kwargs...) where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, O<:AbstractOperator{B}}
    # check if it can be set in 1
    found_param, changed_matrix = set_parameter!(operator.op, parameter, value, print_result=print_result, recalculate=recalculate; kwargs...)
    if parameter == operator.label
        found_param_i = true
        operator.factor = value
    else
        found_param_i = false
    end
    if recalculate && (found_param_i || changed_matrix)
        recalculate!(operator, false, false)
    end
    return (found_param || found_param_i, found_param_i || changed_matrix)
end

# get a parameter (returns (found parameter?, parameter value or nothing))
function get_parameter(operator :: SettableScalarProductOperator{B, O}, parameter :: Symbol; print_result::Bool=false, kwargs...) where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, O<:AbstractOperator{B}}
    # check if parameter can be obtained
    param_1 = (parameter == operator.label) ? operator.factor : nothing
    param_2 = get_parameter(operator.op, parameter, print_result=print_result; kwargs...)
    # give warning if paramter found in both
    if param_1 != nothing && param_2 != nothing
        if print_result
            println("WARNING, parameter :$(parameter) found both in sub operator and factor operator. Returning value of factor operator")
        end
        return param_1
    elseif param_1 != nothing
        if print_result
            println("Parameter :$(parameter) found in factor operator, value $(param_1)")
        end
        return param_1
    elseif param_2 != nothing
        if print_result
            println("Parameter :$(parameter) found in sub operator, value $(param_2)")
        end
        return param_2
    else
        if print_result
            println("Parameter :$(parameter) not found in factor operator or sub operator")
        end
    end
end

# get a parameter (returns (found parameter?, parameter value or nothing))
function get_parameters(operator :: SettableScalarProductOperator{B, O}; kwargs...) where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, O<:AbstractOperator{B}}
    # check if parameter can be obtained
    l = get_parameters(operator.op; kwargs...)
    push!(l, operator.label)
    return unique(l)
end




################################################################################
#   (*) functions
################################################################################

# import of base function
import Base.(*)

# Define the (*) function for any operators
function *(f :: Symbol, op :: O)  where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, O<:AbstractOperator{B}}
    return SettableScalarProductOperator(f, op)
end
function *(op :: O, f :: Symbol)  where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, O<:AbstractOperator{B}}
    return SettableScalarProductOperator(f, op)
end


# Define the (*) function for any operators
function *(f :: Pair{Symbol,<:Number}, op :: O)  where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, O<:AbstractOperator{B}}
    return SettableScalarProductOperator(f[1], f[2], op)
end
function *(op :: O, f :: Pair{Symbol,<:Number})  where {BS<:AbstractBasisState, B<:AbstractBasis{BS}, O<:AbstractOperator{B}}
    return SettableScalarProductOperator(f[1], f[2], op)
end
