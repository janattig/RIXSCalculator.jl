################################################################################
#
#   DISTORTION OPERATOR
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


# Type Definition of Ln^2 / Distortion Operator
"""
   mutable struct DistortionOperator{SPB} <: AbstractSPSSOperator{SPB}

This object refers to the Crystal Distortion Operator.

# Fields
- `basis::SPB`, the single particle basis it refers to;
- `matrix_rep :: Matrix{Complex{Float64}}`, its matrix representation;
- `Delta :: Float64`, the distortion strength;
- `n :: Vector{Float64}`, the (normalized) distortion direction.
"""
mutable struct DistortionOperator{SPB} <: AbstractSPSSOperator{SPB}
    # the basis
    basis :: SPB
    # the current matrix representation (without prefactor)
    matrix_rep :: Matrix{Complex{Float64}}
    # the distortion strength
    Delta :: Float64
    # the distortion direction (normalized)
    n :: Vector{Float64}
end


# Custom constructor (without explicit matrix rep) separate strength and direction
"""
    DistortionOperator(basis::SPB,
                       Delta::Real,
                       n::Vector{<:Real}=[0,0,1])
                       where {SPSSBS<:AbstractSPSSBasisState, SPB<:SPBasis{SPSSBS}}

This function computes the matrix representation of the Crystal Distortion Operator projected over the given `basis`. The crystal distortion direction `n` is given as [0,0,1] by default.
"""
function DistortionOperator(basis::SPB, Delta::Real, n::Vector{<:Real}=[0,0,1]) where {SPSSBS<:AbstractSPSSBasisState, SPB<:SPBasis{SPSSBS}}
    # construct new default basis
    basis_internal = getT2GBasisLS()
    op = DistortionOperator{SPBasis{BasisStateLS}}(basis_internal, zeros(Complex{Float64}, length(basis_internal), length(basis_internal)), Delta, n./norm(n))
    # recalculate the matrix representation
    recalculate!(op)
    # build a projection operator around it
    op_proj = SPSSProjectorOperator(op, basis)
    # return the operator
    return op_proj
end
function DistortionOperator(basis::SPB, Delta::Real, n::Vector{<:Real}=[0,0,1]) where {SPB<:SPBasis{BasisStateLS}}
    # construct new operator
    op = DistortionOperator{SPB}(basis, zeros(Complex{Float64}, length(basis), length(basis)), Delta, n./norm(n))
    # recalculate the matrix representation
    recalculate!(op)
    # return the operator
    return op
end

# creating a distortion operator on a multi site basis
function DistortionOperator(basis::SPMSB, site::Int64, Delta::Real, n::Vector{<:Real}=[0,0,1]) where {SPSSBS<:AbstractSPSSBasisState, SPMSB<:SPBasis{SPMSBasisState{SPSSBS}}}
    # construct new single site operator
    op = DistortionOperator(getSingleSiteBasis(basis, site), Delta, n)
    # construct new multi site operator and return it
    return SPLocalMSOperator(basis, op, site)
end


# export operator type
export  DistortionOperator




import Base.show
function Base.show(io::IO, op::OP) where {SPBS<:AbstractSPSSBasisState, SPB<:SPBasis{SPBS}, OP <: DistortionOperator{SPB}}
    if haskey(io, :compact)
        print(io, "("*replace(string(SPBS), "BasisState"=>"")*") Distortion operator  "*string(round(op.Delta, digits=2))*" (L*n)^2 with n || "*string(round.(op.n, digits=2)))
    else
        print(io, "Distortion operator ((L*n)^2)\n")
        print(io, "Distortion direction n || "*string(round.(op.n, digits=4))*"\n")
        print(io, "Prefactor Delta="*string(round.(op.Delta, digits=4))*"\n")
        print(io, "Basis consists of "*string(length(basis(op)))*" states of type $(SPBS)\n")
    end
end




################################################################################
#   Interface functions
################################################################################

# obtain the current basis
function basis(operator :: DistortionOperator{SPB}) :: SPB where {SPSSBS<:AbstractSPSSBasisState, SPB<:SPBasis{SPSSBS}}
    return operator.basis
end

# obtain the matrix representation
function matrix_representation(operator :: DistortionOperator{SPB}) :: Matrix{Complex{Float64}} where {SPSSBS<:AbstractSPSSBasisState, SPB<:SPBasis{SPSSBS}}
    return operator.matrix_rep .* operator.Delta
end

# possibly recalculate the matrix representation
function recalculate!(operator :: DistortionOperator{SPB}, recursive::Bool=true, basis_change::Bool=true) where {SPSSBS<:AbstractSPSSBasisState, SPB<:SPBasis{SPSSBS}}
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
        operator.matrix_rep[alpha, beta] = getMatrixElementLDotnSquared(basis(operator), basis(operator)[alpha], basis(operator)[beta], operator.n)
    end
    end
end

# set a parameter (returns (found parameter?, changed matrix?))
function set_parameter!(operator :: DistortionOperator{SPB}, parameter :: Symbol, value; print_result::Bool=false, recalculate::Bool=true, kwargs...) where {SPSSBS<:AbstractSPSSBasisState, SPB<:SPBasis{SPSSBS}}
    # check if parameter can be set
    if parameter == :n
        # check if it is only strenth or if it includes direction
        if typeof(value) <: Vector{<:Real}
            # includes direction
            if recalculate && operator.n != value ./ norm(value)
                operator.n = value ./ norm(value)
                recalculate!(operator)
            else
                operator.n = value ./ norm(value)
            end
            # print
            if print_result
                println("Parameter :$(parameter) found and set to value $(operator.n)")
            end
        else
            # parameter found but weird type, error
            @error "recognized parameter :$(parameter) but the type is not matching Vector{Real}: $(typeof(value))" stacktrace()
            return (false, false)
        end
        return (true, true)
    elseif parameter == :Delta
        # includes Prefactor
        if recalculate && operator.Delta != value
            operator.Delta = value
            recalculate!(operator)
        else
            operator.Delta = value
        end
        # print
        if print_result
            println("Parameter :$(parameter) found and set to value $(value)")
        end
        return (true, true)
    else
        if print_result
            println("Parameter :$(parameter) not found")
        end
        return (false, false)
    end
end

# get a parameter (returns (parameter value or nothing))
function get_parameter(operator :: DistortionOperator{SPB}, parameter :: Symbol; print_result::Bool=false, kwargs...) where {SPSSBS<:AbstractSPSSBasisState, SPB<:SPBasis{SPSSBS}}
    # check if parameter can be set
    if parameter == :Delta
        if print_result
            println("Parameter :$(parameter) found and returned its value $(operator.Delta)")
        end
        return operator.Delta
    elseif parameter == :n
        if print_result
            println("Parameter :$(parameter) found and returned its value $(value)")
        end
        return operator.n
    else
        if print_result
            println("Parameter :$(parameter) not found")
        end
        return nothing
    end
end

# get a list of parameters
function get_parameters(operator :: DistortionOperator{SPB}; kwargs...) where {SPSSBS<:AbstractSPSSBasisState, SPB<:SPBasis{SPSSBS}}
    return Symbol[:Delta, :n]
end











################################################################################
#   Definition of matrix element functions
################################################################################

# Simplification function to compute (L*n)^2 (for any basis)
# <s1|Ln^2|s2> = sum_3 <s1|Ln|s3> <s3|Ln|s2>
function getMatrixElementLDotnSquared(basis::SPB, state_1::SPSSBS, state_2::SPSSBS, n::Vector{<:Real}) :: Complex{Float64} where {SPSSBS<:AbstractSPSSBasisState, SPB<:SPBasis{SPSSBS}}
    # build as an explicit product
    term = 0.0 + 0.0im
    # sum over all possible intermediate states in the given basis
    for state_3 in basis
        # sum to the term
        term += getMatrixElementLDotn(state_1, state_3, n) * getMatrixElementLDotn(state_3, state_2, n)
    end
    # return the value of the term
    return term
end





# Term in any basis
"""
    getMatrixElementLDotn(state_1::SPSSBS, state_2::SPSSBS, n::Vector{<:Real}) :: Complex{Float64} where {SPSSBS<:AbstractSPSSBasisState}

This function computes the matrix element `` \\left< state_1 \\left| \\boldsymbol{L} \\cdot \\boldsymbol{n} \\right| state_2 \\right> = \\sum_i \\left< state_1 \\left| L_i n[i] \\right| state_2 \\right>``.
"""
function getMatrixElementLDotn(state_1::SPSSBS, state_2::SPSSBS, n::Vector{<:Real}) :: Complex{Float64} where {SPSSBS<:AbstractSPSSBasisState}
    @error "Matrix element of L*n not yet implement for basis states of type $(SPSSBS)" stacktrace()
    return NaN + NaN*im
end

# Term in L,S Basis
# <s1|Ln|s2> = sum_i <s1|L_i * n[i]|s2>
function getMatrixElementLDotn(state_1::BasisStateLS, state_2::BasisStateLS, n::Vector{<:Real}) :: Complex{Float64}
    # return the value
    return (
          operatorLx(state_1.l, state_1.ml, state_2.ml) * n[1]
        + operatorLy(state_1.l, state_1.ml, state_2.ml) * n[2]
        + operatorLz(state_1.l, state_1.ml, state_2.ml) * n[3]
    ) * delta(state_1.l,state_2.l) * delta(state_1.s,state_2.s) * delta(state_1.ms,state_2.ms)
end
export getMatrixElementLDotn



##############################################################
#   Convenience functions for SS SP operators
##############################################################

# creating a distortion operator on a multi site basis
function DistortionOperator(basis::MPB, site::Int64, Delta::Real, n::Vector{<:Real}=[0,0,1]; particle_type::Symbol=:hole) where {SPSSBS<:AbstractSPSSBasisState, N,MPB<:MPBasis{N,SPMSBasisState{SPSSBS}}}
    # construct new single site operator
    op = DistortionOperator(basis.single_particle_basis, site, Delta, n)
    # construct new multi particle operator out of that
    if particle_type == :electron
        return MPElectronGeneralizedSPOperator(basis, op)
    elseif particle_type == :hole
        return MPHoleGeneralizedSPOperator(basis, op)
    else 
        @error "Invalid particle type '$(particle_type)'; returned ':hole' instead" stacktrace()
        return MPHoleGeneralizedSPOperator(basis, op)
    end
        
end
