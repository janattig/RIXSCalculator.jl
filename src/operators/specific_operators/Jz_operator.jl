# Type Definition of Jz Operator
"""
    mutable struct JzOperator{SPB} <: AbstractSPSSOperator{SPB}

This object refers to the ``J_z`` Operator.

# Fields

- `basis :: SPB`, the single particle basis;
- `matrix_rep :: Matrix{Complex{Float64}}`, the matrix representation of the operator;
- `spin_quantization :: CoordinateFrame`, the spin quantization axis.

"""
mutable struct JzOperator{SPB} <: AbstractSPSSOperator{SPB}
    # the basis
    basis :: SPB
    # the current matrix representation (without prefactor)
    matrix_rep :: Matrix{Complex{Float64}}
    # spin quantization axis
    spin_quantization :: CoordinateFrame
end

# Custom constructor (without explicit matrix rep)
"""
    JzOperator(basis::SPB) where {SPSSBS<:AbstractSPSSBasisState, SPB<:SPBasis{SPSSBS}}

This function computes the matrix representation of the ``J_z`` Operator projected over the given `basis`.
"""
function JzOperator(basis::SPB) where {SPSSBS<:AbstractSPSSBasisState, SPB<:SPBasis{SPSSBS}}
    # construct new operator
    basis_internal = getT2GBasisLS()
    op = JzOperator{SPBasis{BasisStateLS}}(basis_internal, zeros(Complex{Float64}, length(basis_internal), length(basis_internal)), CoordinateFrame())
    # recalculate the matrix representation
    recalculate!(op)
    # build a projection operator around it
    op_proj = SPSSProjectorOperator(op, basis)
    # return the operator
    return op_proj
end

# Custom constructor (without explicit matrix rep)
function JzOperator(basis::SPB) where {SPB<:SPBasis{BasisStateLS}}
    # construct new operator
    op = JzOperator{SPB}(basis, zeros(Complex{Float64}, length(basis), length(basis)), CoordinateFrame())
    # recalculate the matrix representation
    recalculate!(op)
    # return the operator
    return op
end

# creating a Jz operator on a multi site basis
function JzOperator(basis::SPMSB, site::Int64) where {SPSSBS<:AbstractSPSSBasisState, SPMSB<:SPBasis{SPMSBasisState{SPSSBS}}}
    # construct new single site operator
    op = JzOperator(getSingleSiteBasis(basis, site))
    # construct new multi site operator and return it
    return SPLocalMSOperator(basis, op, site)
end

# export operator type
export  JzOperator




import Base.show
function Base.show(io::IO, op::OP) where {SPBS<:AbstractSPSSBasisState, SPB<:SPBasis{SPBS}, OP <: JzOperator{SPB}}
    if haskey(io, :compact)
        print(io, "("*replace(string(SPBS), "BasisState"=>"")*") Jz operator ")
    else
        print(io, "Jz operator (Lz+Sz)\n")
        print(io, "Basis consists of "*string(length(basis(op)))*" states of type $(SPBS)\n")
    end
end


################################################################################
#   Interface functions
################################################################################
import RIXSCalculator.basis, RIXSCalculator.matrix_representation, RIXSCalculator.recalculate!, RIXSCalculator.set_parameter!, RIXSCalculator.get_parameter, RIXSCalculator.get_parameters
# obtain the current basis
function basis(operator :: JzOperator{SPB}) :: SPB where {SPSSBS<:AbstractSPSSBasisState, SPB<:SPBasis{SPSSBS}}
    return operator.basis
end

# obtain the matrix representation
function matrix_representation(operator :: JzOperator{SPB}) :: Matrix{Complex{Float64}} where {SPSSBS<:AbstractSPSSBasisState, SPB<:SPBasis{SPSSBS}}
    return operator.matrix_rep
end

# possibly recalculate the matrix representation
function recalculate!(operator :: JzOperator{SPB}, recursive::Bool=true, basis_change::Bool=true) where {SPSSBS<:AbstractSPSSBasisState, SPB<:SPBasis{SPSSBS}}
    # check if the size of matrix is still okay
    if size(operator.matrix_rep,1) == length(basis(operator)) && size(operator.matrix_rep,2) == length(basis(operator))
        # size is okay, multiply matrix by 0 to erase all elements
        operator.matrix_rep .*= 0.0
    else
        # create new matrix
        operator.matrix_rep = zeros(Complex{Float64}, length(basis(operator)), length(basis(operator)))
    end
    # recalculate the matrix elements
    for alpha in 1:length(basis(operator))
    for beta in 1:length(basis(operator))
        # set the respective entry
        operator.matrix_rep[alpha, beta] = getMatrixElementJz(basis(operator)[alpha], basis(operator)[beta], operator.spin_quantization)
    end
    end
end

# set a parameter (returns (found parameter?, changed matrix?))
function set_parameter!(operator :: JzOperator{SPB}, parameter :: Symbol, value; print_result::Bool=false, recalculate::Bool=true, kwargs...) where {SPSSBS<:AbstractSPSSBasisState, SPB<:SPBasis{SPSSBS}}
    # check if parameter can be set
    if parameter == :spin_axis || parameter == :spin_quantization
        # set value
        if recalculate
            try
                operator.spin_quantization = value
            catch Exception
                @error "recognized parameter :$(parameter) but the type is not matching CoordinateFrame: $(typeof(value))" stacktrace()
            end
            recalculate!(operator)
        else
            try
                operator.spin_quantization = value
            catch Exception
                @error "recognized parameter :$(parameter) but the type is not matching CoordinateFrame: $(typeof(value))" stacktrace()
            end
        end
        # print
        if print_result
            println("Parameter :$(parameter) found and set to value")
        end
        return (true, true)
    else
        if print_result
            println("Parameter :$(parameter) not found")
        end
        return (false, false)
    end
end



# get a parameter (returns (found parameter?, parameter value or nothing))
function get_parameter(operator :: JzOperator{SPB}, parameter :: Symbol; print_result::Bool=false, kwargs...) where {SPSSBS<:AbstractSPSSBasisState, SPB<:SPBasis{SPSSBS}}
    # check if parameter can be set

    if print_result
        println("Parameter :$(parameter) not found")
    end
    return nothing
end

# get a list of parameters
function get_parameters(operator :: JzOperator{SPB}; kwargs...) where {SPSSBS<:AbstractSPSSBasisState, SPB<:SPBasis{SPSSBS}}
    return Symbol[]
end





################################################################################
#   Definition of matrix element functions
################################################################################

# Term in any basis
"""
    getMatrixElementJz(state_1::SPSSBS, state_2::SPSSBS) :: Complex{Float64} where {SPSSBS<:AbstractSPSSBasisState}
This function computes the matrix element `` \\left< state_1 \\left| J_z \\right| state_2 \\right>``.
"""
function getMatrixElementJz(state_1::SPSSBS, state_2::SPSSBS) :: Complex{Float64} where {SPSSBS<:AbstractSPSSBasisState}
    @error "Matrix element of Jz not yet implement for basis states of type $(SPSSBS)" stacktrace()
    return NaN + NaN*im
end
# Term in any basis
function getMatrixElementJz(state_1::SPSSBS, state_2::SPSSBS, cf::CoordinateFrame) :: Complex{Float64} where {SPSSBS<:AbstractSPSSBasisState}
    @error "Matrix element of Jz not yet implement for basis states of type $(SPSSBS)" stacktrace()
    return NaN + NaN*im
end

# Term in L,S Basis
# <s1|BS|s2>
function getMatrixElementJz(state_1::BasisStateLS, state_2::BasisStateLS) :: Complex{Float64}
    # build the return value
    return (
        operatorLz(state_1.l, state_1.ml, state_2.ml) * delta(state_1.ms,state_2.ms) +
        operatorSz(state_1.s, state_1.ms, state_2.ms) * delta(state_1.ml,state_2.ml)
    ) * delta(state_1.l,state_2.l)  * delta(state_1.s,state_2.s) 
end
# Term in L,S Basis
# <s1|BS|s2>
function getMatrixElementJz(state_1::BasisStateLS, state_2::BasisStateLS, cf::CoordinateFrame) :: Complex{Float64}
    # build the return value
    return (
        operatorLz(state_1.l, state_1.ml, state_2.ml)     * delta(state_1.ms,state_2.ms) +
        operatorSz(state_1.s, state_1.ms, state_2.ms, cf) * delta(state_1.ml,state_2.ml)
    ) * delta(state_1.l,state_2.l)  * delta(state_1.s,state_2.s) 
end

export getMatrixElementJz



##############################################################
#   Convenience functions for SS SP operators
##############################################################

# creating a magnetic field operator on a multi site basis
function JzOperator(basis::MPB, site::Int64; particle_type::Symbol=:hole) where {SPSSBS<:AbstractSPSSBasisState, N,MPB<:Union{MPBasis{N,SPMSBasisState{SPSSBS}}, MPBasis{N, TdSymBasisState}}}
    # construct new single site operator
    op = JzOperator(basis.single_particle_basis, site)
    # construct new multi particle operator out of that
    if particle_type == :electron
        return MPElectronGeneralizedSPOperator(basis, op)
    elseif particle_type == :hole
        return MPHoleGeneralizedSPOperator(basis, op)
    else 
        @error "Invalid particle type '$(particle_type)'; returned ':hole' version instead" stacktrace()
        return MPHoleGeneralizedSPOperator(basis, op)
    end
end