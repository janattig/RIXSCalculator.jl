################################################################################
#
#   MAGNETIC FIELD OPERATOR
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


# Type Definition of BS / Magnetic field Operator
"""
    MagneticFieldOperator{SPB} <: AbstractSPSSOperator{SPB}

This object refers to the Magnetic Field Operator.

It is characterized by the single particle `basis::SPB` it refers to, its `matrix_rep :: Matrix{Complex{Float64}}`, the field strength `B::Float64`, its direction `B_dir::Vector{Float64}` and the spin quantization axis `spin_quantization :: CoordinateFrame`.
"""
mutable struct MagneticFieldOperator{SPB} <: AbstractSPSSOperator{SPB}
    # the basis
    basis :: SPB
    # the current matrix representation (without prefactor)
    matrix_rep :: Matrix{Complex{Float64}}
    # the field strength
    B :: Float64
    # the field direction (normalized)
    B_dir :: Vector{Float64}
    # spin quantisation axis
    spin_quantization :: CoordinateFrame
end

# Custom constructor (without explicit matrix rep)
"""
    MagneticFieldOperator(basis::SPB, B::Vector{<:Real}) where {SPSSBS<:AbstractSPSSBasisState, SPB<:SPBasis{SPSSBS}}

This function computes the matrix representation of the Magnetic Field Operator projected over the given `basis`.
"""
function MagneticFieldOperator(basis::SPB, B::Vector{<:Real}) where {SPSSBS<:AbstractSPSSBasisState, SPB<:SPBasis{SPSSBS}}
    # construct new operator
    basis_internal = getT2GBasisLS()
    op = MagneticFieldOperator{SPBasis{BasisStateLS}}(basis_internal, zeros(Complex{Float64}, length(basis_internal), length(basis_internal)), norm(B), B./norm(B), CoordinateFrame())
    # recalculate the matrix representation
    recalculate!(op)
    # build a projection operator around it
    op_proj = SPSSProjectorOperator(op, basis)
    # return the operator
    return op_proj
end
# Custom constructor (without explicit matrix rep) separate strength and direction
function MagneticFieldOperator(basis::SPB, B::Real, B_dir::Vector{<:Real}=[0,0,1]) where {SPSSBS<:AbstractSPSSBasisState, SPB<:SPBasis{SPSSBS}}
    # construct new operator
    basis_internal = getT2GBasisLS()
    op = MagneticFieldOperator{SPBasis{BasisStateLS}}(basis_internal, zeros(Complex{Float64}, length(basis_internal), length(basis_internal)), B, B_dir./norm(B_dir), CoordinateFrame())
    # recalculate the matrix representation
    recalculate!(op)
    # build a projection operator around it
    op_proj = SPSSProjectorOperator(op, basis)
    # return the operator
    return op_proj
end



# Custom constructor (without explicit matrix rep)
function MagneticFieldOperator(basis::SPB, B::Vector{<:Real}) where {SPB<:SPBasis{BasisStateLS}}
    # construct new operator
    op = MagneticFieldOperator{SPB}(basis, zeros(Complex{Float64}, length(basis), length(basis)), norm(B), B./norm(B), CoordinateFrame())
    # recalculate the matrix representation
    recalculate!(op)
    # return the operator
    return op
end
# Custom constructor (without explicit matrix rep) separate strength and direction
function MagneticFieldOperator(basis::SPB, B::Real, B_dir::Vector{<:Real}=[0,0,1]) where {SPB<:SPBasis{BasisStateLS}}
    # construct new operator
    op = MagneticFieldOperator{SPB}(basis, zeros(Complex{Float64}, length(basis), length(basis)), B, B_dir./norm(B_dir), CoordinateFrame())
    # recalculate the matrix representation
    recalculate!(op)
    # return the operator
    return op
end


# creating a magnetic field operator on a multi site basis
function MagneticFieldOperator(basis::SPMSB, site::Int64, B::Real, B_dir::Vector{<:Real}=[0,0,1]) where {SPSSBS<:AbstractSPSSBasisState, SPMSB<:SPBasis{SPMSBasisState{SPSSBS}}}
    # construct new single site operator
    op = MagneticFieldOperator(getSingleSiteBasis(basis, site), B, B_dir)
    # construct new multi site operator and return it
    return SPLocalMSOperator(basis, op, site)
end



# export operator type
export  MagneticFieldOperator


import Base.show
function Base.show(io::IO, op::OP) where {SPBS<:AbstractSPSSBasisState, SPB<:SPBasis{SPBS}, OP <: MagneticFieldOperator{SPB}}
    if haskey(io, :compact)
        print(io, "("*replace(string(SPBS), "BasisState"=>"")*") Magnetic field operator  "*string(round(op.B, digits=2))*" B*S with B || "*string(round.(op.B_dir, digits=2)))
    else
        print(io, "Magnetic field operator (B*S)\n")
        print(io, "Magnetic field B || "*string(round.(op.B_dir, digits=4))*"\n")
        print(io, "Field Strength |B|="*string(round.(op.B, digits=4))*"\n")
        print(io, "Basis consists of "*string(length(basis(op)))*" states of type $(SPBS)\n")
    end
end




################################################################################
#   Interface functions
################################################################################

# obtain the current basis
function basis(operator :: MagneticFieldOperator{SPB}) :: SPB where {SPSSBS<:AbstractSPSSBasisState, SPB<:SPBasis{SPSSBS}}
    return operator.basis
end

# obtain the matrix representation
function matrix_representation(operator :: MagneticFieldOperator{SPB}) :: Matrix{Complex{Float64}} where {SPSSBS<:AbstractSPSSBasisState, SPB<:SPBasis{SPSSBS}}
    return operator.matrix_rep .* operator.B
end

# possibly recalculate the matrix representation
function recalculate!(operator :: MagneticFieldOperator{SPB}, recursive::Bool=true, basis_change::Bool=true) where {SPSSBS<:AbstractSPSSBasisState, SPB<:SPBasis{SPSSBS}}
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
        operator.matrix_rep[alpha, beta] = getMatrixElementBDotS(basis(operator)[alpha], basis(operator)[beta], operator.B_dir, operator.spin_quantization)
    end
    end
end

# set a parameter (returns (found parameter?, changed matrix?))
function set_parameter!(operator :: MagneticFieldOperator{SPB}, parameter :: Symbol, value; print_result::Bool=false, recalculate::Bool=true, kwargs...) where {SPSSBS<:AbstractSPSSBasisState, SPB<:SPBasis{SPSSBS}}
    # check if parameter can be set
    if parameter == :B
        # check if it is only strenth or if it includes direction
        if typeof(value) <: Vector{<:Real}
            # includes direction
            if recalculate && (operator.B != norm(value) || operator.B_dir != value ./ norm(value))
                operator.B = norm(value)
                operator.B_dir = value ./ norm(value)
                recalculate!(operator)
            else
                operator.B = norm(value)
                operator.B_dir = value ./ norm(value)
            end
            # print
            if print_result
                println("Parameter :$(parameter) found and set to value $(operator.B) with direction $(operator.B_dir)")
            end
        elseif typeof(value) <: Real
            # does not include direction
            if recalculate && operator.B != norm(value)
                operator.B = norm(value)
                recalculate!(operator)
            else
                operator.B = norm(value)
            end
            # print
            if print_result
                println("Parameter :$(parameter) found and set to value $(operator.B)")
            end
        else
            # parameter found but weird type, error
            @error "recognized parameter :$(parameter) but the type is not matching Vector{Real} or Real: $(typeof(value))" stacktrace()
            return (false, false)
        end
        return (true, true)# check if parameter can be set
    elseif parameter == :B_dir
        # set value
        if recalculate && operator.B_dir != value ./ norm(value)
            operator.B_dir = value ./ norm(value)
            recalculate!(operator)
        else
            operator.B_dir = value ./ norm(value)
        end
        # print
        if print_result
            println("Parameter :$(parameter) found and set to value $(operator.B_dir)")
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
function get_parameter(operator :: MagneticFieldOperator{SPB}, parameter :: Symbol; print_result::Bool=false, kwargs...) where {SPSSBS<:AbstractSPSSBasisState, SPB<:SPBasis{SPSSBS}}
    # check if parameter can be set
    if parameter == :B
        if print_result
            println("Parameter :$(parameter) found and returned its value $(operator.B)")
        end
        return operator.B
    elseif parameter == :B_dir
        if print_result
            println("Parameter :$(parameter) found and returned its value $(operator.B_dir)")
        end
        return operator.B_dir
    else
        if print_result
            println("Parameter :$(parameter) not found")
        end
        return nothing
    end
end

# get a list of parameters
function get_parameters(operator :: MagneticFieldOperator{SPB}; kwargs...) where {SPSSBS<:AbstractSPSSBasisState, SPB<:SPBasis{SPSSBS}}
    return Symbol[:B, :B_dir]
end









################################################################################
#   Definition of matrix element functions
################################################################################

# Term in any basis
"""
    getMatrixElementBDotS(state_1::SPSSBS, state_2::SPSSBS, B_dir::Vector{<:Real}) :: Complex{Float64} where {SPSSBS<:AbstractSPSSBasisState}

This function computes the matrix element `` \\left< state_1 \\left| \\boldsymbol{B} \\cdot \\boldsymbol{S} \\right| state_2 \\right> = \\sum_i \\left< state_1 \\left| B[i] S_i \\right| state_2 \\right>``.
"""
function getMatrixElementBDotS(state_1::SPSSBS, state_2::SPSSBS, B_dir::Vector{<:Real}) :: Complex{Float64} where {SPSSBS<:AbstractSPSSBasisState}
    @error "Matrix element of B*S not yet implement for basis states of type $(SPSSBS)" stacktrace()
    return NaN + NaN*im
end
# Term in any basis
function getMatrixElementBDotS(state_1::SPSSBS, state_2::SPSSBS, B_dir::Vector{<:Real}, cf::CoordinateFrame) :: Complex{Float64} where {SPSSBS<:AbstractSPSSBasisState}
    @error "Matrix element of B*S not yet implement for basis states of type $(SPSSBS)" stacktrace()
    return NaN + NaN*im
end

# Term in L,S Basis
# <s1|BS|s2>
function getMatrixElementBDotS(state_1::BasisStateLS, state_2::BasisStateLS, B_dir::Vector{<:Real}) :: Complex{Float64}
    # build the return value
    return (
        operatorSx(state_1.s, state_1.ms, state_2.ms) *  B_dir[1] +
        operatorSy(state_1.s, state_1.ms, state_2.ms) *  B_dir[2] +
        operatorSz(state_1.s, state_1.ms, state_2.ms) *  B_dir[3]
    ) * delta(state_1.l,state_2.l) * delta(state_1.ml,state_2.ml) * delta(state_1.s,state_2.s)
end
# Term in L,S Basis
# <s1|BS|s2>
function getMatrixElementBDotS(state_1::BasisStateLS, state_2::BasisStateLS, B_dir::Vector{<:Real}, cf::CoordinateFrame) :: Complex{Float64}
    # build the return value
    return (
        operatorSx(state_1.s, state_1.ms, state_2.ms, cf) *  B_dir[1] +
        operatorSy(state_1.s, state_1.ms, state_2.ms, cf) *  B_dir[2] +
        operatorSz(state_1.s, state_1.ms, state_2.ms, cf) *  B_dir[3]
    ) * delta(state_1.l,state_2.l) * delta(state_1.ml,state_2.ml) * delta(state_1.s,state_2.s)
end




##############################################################
#   Convenience functions for SS SP operators
##############################################################

# creating a magnetic field operator on a multi site basis
function MagneticFieldOperator(basis::MPB, site::Int64, B::Real, B_dir::Vector{<:Real}=[0,0,1]) where {SPSSBS<:AbstractSPSSBasisState, N,MPB<:MPBasis{N,SPMSBasisState{SPSSBS}}}
    # construct new single site operator
    op = MagneticFieldOperator(basis.single_particle_basis, site, B, B_dir)
    # construct new multi particle operator out of that
    return MPGeneralizedSPOperator(basis, op)
end
