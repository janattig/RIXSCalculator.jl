################################################################################
#
#   SPIN ORBIT OPERATOR
#   (SINGLE PARTICLE / SINGLE SITE)
#
#   - type definition
#   - interface functions
#   - definition of matrix element functions
#
################################################################################




################################################################################
#   Type definition
################################################################################


# Type Definition of LS / SpinOrbit Operator
"""
    mutable struct SpinOrbitOperator{SPB} <: AbstractSPSSOperator{SPB}

This object refers to the Spin Orbit Operator.

# Fields

- `basis :: SPB`, the single particle basis;
- `matrix_rep :: Matrix{Complex{Float64}}`, the matrix representation of the operator;
- `lambda :: Float64`, the coupling strength;
- `spin_quantization :: CoordinateFrame`, the spin quantization axis.

"""
mutable struct SpinOrbitOperator{SPB} <: AbstractSPSSOperator{SPB}
    # the basis
    basis :: SPB
    # the current matrix representation (without prefactor)
    matrix_rep :: Matrix{Complex{Float64}}
    # the prefactor
    lambda :: Float64
    # spin quantisation axis
    spin_quantization :: CoordinateFrame
end



# Custom constructor (without explicit matrix rep)
"""
    SpinOrbitOperator(basis::SPB, lambda::Real) where {SPB<:SPBasis{BasisStateLS}}

This function computes the matrix representation of the Spin Orbit Operator with coupling strength `lambda` projected over the given `basis`.
"""
function SpinOrbitOperator(basis::SPB, lambda::Real) where {SPB<:SPBasis{BasisStateLS}}
    # construct new operator
    op = SpinOrbitOperator{SPB}(basis, zeros(Complex{Float64}, length(basis), length(basis)), lambda, CoordinateFrame())
    # recalculate the matrix representation
    recalculate!(op)
    # return the operator
    return op
end
function SpinOrbitOperator(basis::SPB, lambda::Real) where {SPSSBS<:AbstractSPSSBasisState, SPB<:SPBasis{SPSSBS}}
    # construct new operator
    basis_internal = getT2GBasisLS()
    op = SpinOrbitOperator{SPBasis{BasisStateLS}}(basis_internal, zeros(Complex{Float64}, length(basis_internal), length(basis_internal)), lambda, CoordinateFrame())
    # recalculate the matrix representation
    recalculate!(op)
    # build a projection operator around it
    op_proj = SPSSProjectorOperator(op, basis)
    # return the operator
    return op_proj
end

# creating a spin orbit operator on a multi site basis
function SpinOrbitOperator(basis::SPMSB, site::Int64, lambda::Real) where {SPSSBS<:AbstractSPSSBasisState, SPMSB<:SPBasis{SPMSBasisState{SPSSBS}}}
    # construct new single site operator
    op = SpinOrbitOperator(getSingleSiteBasis(basis, site), lambda)
    # construct new multi site operator and return it
    return SPLocalMSOperator(basis, op, site)
end

# export operator type
export  SpinOrbitOperator

import Base.show
function Base.show(io::IO, op::OP) where {SPBS<:AbstractSPSSBasisState, SPB<:SPBasis{SPBS}, OP <: SpinOrbitOperator{SPB}}
    if haskey(io, :compact)
        print(io, "("*replace(string(SPBS), "BasisState"=>"")*") Spin-Orbit operator  "*string(round(op.lambda, digits=2))*" L*S")
    else
        print(io, "Spin-orbit operator (L*S)\n")
        print(io, "Spin-orbit coupling lambda="*string(round.(op.lambda, digits=4))*"\n")
        print(io, "Basis consists of "*string(length(basis(op)))*" states of type $(SPBS)\n")
    end
end





################################################################################
#   Interface functions
################################################################################

# obtain the current basis
function basis(operator :: SpinOrbitOperator{SPB}) :: SPB where {SPSSBS<:AbstractSPSSBasisState, SPB<:SPBasis{SPSSBS}}
    return operator.basis
end

# obtain the matrix representation
function matrix_representation(operator :: SpinOrbitOperator{SPB}) :: Matrix{Complex{Float64}} where {SPSSBS<:AbstractSPSSBasisState, SPB<:SPBasis{SPSSBS}}
    return operator.matrix_rep .* operator.lambda
end

# possibly recalculate the matrix representation
function recalculate!(operator :: SpinOrbitOperator{SPB}, recursive::Bool=true, basis_change::Bool=true) where {SPSSBS<:AbstractSPSSBasisState, SPB<:SPBasis{SPSSBS}}
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
        operator.matrix_rep[alpha, beta] = getMatrixElementLDotS(basis(operator)[alpha], basis(operator)[beta], operator.spin_quantization)
    end
    end
end

# set a parameter (returns (found parameter?, changed matrix?))
function set_parameter!(operator :: SpinOrbitOperator{SPB}, parameter :: Symbol, value; print_result::Bool=false, recalculate::Bool=true, kwargs...) where {SPSSBS<:AbstractSPSSBasisState, SPB<:SPBasis{SPSSBS}}
    # check if parameter can be set
    if parameter == :lambda
        if recalculate && operator.lambda != value
            operator.lambda = value
            recalculate!(operator)
        else
            operator.lambda = value
        end
        if print_result
            println("Parameter :$(parameter) found and set to value $(operator.lambda)")
        end
        return (true, true)
    elseif parameter == :spin_axis || parameter == :spin_quantization
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
function get_parameter(operator :: SpinOrbitOperator{SPB}, parameter :: Symbol; print_result::Bool=false, kwargs...) where {SPSSBS<:AbstractSPSSBasisState, SPB<:SPBasis{SPSSBS}}
    # check if parameter can be set
    if parameter == :lambda
        if print_result
            println("Parameter :$(parameter) found and returned its value $(operator.lambda)")
        end
        return operator.lambda
    elseif parameter == :spin_quantization
        if print_result
            println("Parameter :$(spin_quantization) found and returned its value $(operator.spin_quantization)")
        end
        return operator.spin_quantization
    else
        if print_result
            println("Parameter :$(parameter) not found")
        end
        return nothing
    end
end

# get a list of parameters
function get_parameters(operator :: SpinOrbitOperator{SPB}; kwargs...) where {SPSSBS<:AbstractSPSSBasisState, SPB<:SPBasis{SPSSBS}}
    return Symbol[:lambda,]
end








################################################################################
#   Definition of matrix element functions
################################################################################

# Term in any basis
"""
    getMatrixElementLDotS(state_1::SPSSBS, state_2::SPSSBS) :: Complex{Float64} where {SPSSBS<:AbstractSPSSBasisState}


This function computes the matrix element `` \\left< state_1 \\left| \\boldsymbol{L} \\cdot \\boldsymbol{S} \\right| state_2 \\right> = \\sum_i \\left< state_1 \\left| L_i S_i \\right| state_2 \\right>``.
"""
function getMatrixElementLDotS(state_1::SPSSBS, state_2::SPSSBS) :: Complex{Float64} where {SPSSBS<:AbstractSPSSBasisState}
    @error "Matrix element of L*S not yet implement for basis states of type $(SPSSBS)" stacktrace()
    return NaN + NaN*im
end
# Term in any basis
function getMatrixElementLDotS(state_1::SPSSBS, state_2::SPSSBS, cf::CoordinateFrame) :: Complex{Float64} where {SPSSBS<:AbstractSPSSBasisState}
    @error "Matrix element of L*S not yet implement for basis states of type $(SPSSBS)" stacktrace()
    return NaN + NaN*im
end

# Term in L,S Basis
# <s1|LS|s2>
function getMatrixElementLDotS(state_1::BasisStateLS, state_2::BasisStateLS) :: Complex{Float64}
    # build the return value
    return (
        operatorLx(state_1.l, state_1.ml, state_2.ml) * operatorSx(state_1.s, state_1.ms, state_2.ms) +
        operatorLy(state_1.l, state_1.ml, state_2.ml) * operatorSy(state_1.s, state_1.ms, state_2.ms) +
        operatorLz(state_1.l, state_1.ml, state_2.ml) * operatorSz(state_1.s, state_1.ms, state_2.ms)
    ) * delta(state_1.l,state_2.l) * delta(state_1.s,state_2.s)
end

# Term in L,S Basis
# <s1|LS|s2>
# AND COORDINATE FRAME FOR SPIN QUANTIZATION
function getMatrixElementLDotS(state_1::BasisStateLS, state_2::BasisStateLS, cf::CoordinateFrame) :: Complex{Float64}
    # build the return value
    return (
        operatorLx(state_1.l, state_1.ml, state_2.ml) * operatorSx(state_1.s, state_1.ms, state_2.ms, cf) +
        operatorLy(state_1.l, state_1.ml, state_2.ml) * operatorSy(state_1.s, state_1.ms, state_2.ms, cf) +
        operatorLz(state_1.l, state_1.ml, state_2.ml) * operatorSz(state_1.s, state_1.ms, state_2.ms, cf)
    ) * delta(state_1.l,state_2.l) * delta(state_1.s,state_2.s)
end


##############################################################
#   Convenience functions for SS SP operators
##############################################################


# creating a spin orbit operator on a multi site basis
function SpinOrbitOperator(basis::MPB, site::Int64, lambda::Real) where {SPSSBS<:AbstractSPSSBasisState, N,MPB<:MPBasis{N,SPMSBasisState{SPSSBS}}}
    # construct new single site local operator
    op = SpinOrbitOperator(basis.single_particle_basis, site, lambda)
    # construct new multi particle operator out of that
    return MPGeneralizedSPOperator(basis, op)
end
