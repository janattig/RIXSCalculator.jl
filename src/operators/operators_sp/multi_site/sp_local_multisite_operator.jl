# Type definition of local MS operators (SPSS operator only acting on a single site of many)
"""
    SPLocalMSOperator{
        SPSSBS <: AbstractSPSSBasisState,
        SPMSB <: SPBasis{SPMSBS} where {SPMSBS<:SPMSBasisState{SPSSBS}},
        SPO <: AbstractSPSSOperator{SPBasis{SPSSBS}}
    } <: AbstractSPMSOperator{SPMSB}

This object refers to the multisite description of a single particle local operator.

It is characterized by the single particle multi site`basis::SPMSB` it refers to, the said single site operator `operator::SPO`, the site `site::Int64` it acts on, the basis `basis_ss   :: SPBasis{SPSSBS}` and its indices `indices_ss :: Vector{Int64}` in the multi site basis and
its `matrix_rep :: Matrix{Complex{Float64}}`.
"""
mutable struct SPLocalMSOperator{
        SPSSBS <: AbstractSPSSBasisState,
        SPMSB <: SPBasis{SPMSBS} where {SPMSBS<:SPMSBasisState{SPSSBS}},
        SPO <: AbstractSPSSOperator{SPBasis{SPSSBS}}
    } <: AbstractSPMSOperator{SPMSB}

    # the MS basis
    basis :: SPMSB
    # the single site operator
    operator :: SPO
    # the site that it acts on
    site :: Int64
    # the basis and its indices in the multi site basis
    basis_ss   :: SPBasis{SPSSBS}
    indices_ss :: Vector{Int64}
    # the matrix representation
    matrix_rep :: Matrix{Complex{Float64}}
end

# custom constructor
"""
    SPLocalMSOperator(basis :: SPMSB, operator :: SPO, site :: Int64) where {SPSSBS <: AbstractSPSSBasisState, SPMSBS <: SPMSBasisState{SPSSBS}, SPMSB <: SPBasis{SPMSBS}, SPO <: AbstractSPSSOperator{SPBasis{SPSSBS}}}

This function computes the matrix representation of the single particle - multi site local Operator.
"""
function SPLocalMSOperator(basis :: SPMSB, operator :: SPO, site :: Int64) where {SPSSBS <: AbstractSPSSBasisState, SPMSBS <: SPMSBasisState{SPSSBS}, SPMSB <: SPBasis{SPMSBS}, SPO <: AbstractSPSSOperator{SPBasis{SPSSBS}}}
    # create a new operator
    op = SPLocalMSOperator{SPSSBS, SPMSB, SPO}(basis, operator, site, getSingleSiteBasis(basis, site), Int64[i     for i in 1:length(basis) if basis[i].site == site], zeros(length(basis), length(basis)))
    # recalculate the matrix representation
    recalculate!(op, true)
    # return the operator
    return op
end

# export operator type
export  SPLocalMSOperator

import Base.show
function Base.show(io::IO, op::SPLocalMSOperator{SPSSBS, SPMSB, SPO}) where {SPSSBS <: AbstractSPSSBasisState, SPMSBS <: SPMSBasisState{SPSSBS}, SPMSB <: SPBasis{SPMSBS}, SPO <: AbstractSPSSOperator{SPBasis{SPSSBS}}}
    if haskey(io, :compact)
        print(io, "@site("*string(op.site)*") ")
        show(io, op.operator)
    else
        print(io, "Multi-site operator on site "*string(op.site)*" of type ")
        show(io, op.operator)
        print(io, "Multi-site basis contains "*string(length(basis(op)))*" states in total, defined on sites "*string(unique([b.site for b in basis(op)]))*"\n")
    end
end




################################################################################
#   Interface functions
################################################################################

# obtain the current basis
function basis(operator :: SPLocalMSOperator{SPSSBS, SPMSB, SPO}) :: SPMSB where {SPSSBS <: AbstractSPSSBasisState, SPMSBS <: SPMSBasisState{SPSSBS}, SPMSB <: SPBasis{SPMSBS}, SPO <: AbstractSPSSOperator{SPBasis{SPSSBS}}}
    return operator.basis
end

# obtain the matrix representation
function matrix_representation(operator :: SPLocalMSOperator{SPSSBS, SPMSB, SPO}) :: Matrix{Complex{Float64}} where {SPSSBS <: AbstractSPSSBasisState, SPMSBS <: SPMSBasisState{SPSSBS}, SPMSB <: SPBasis{SPMSBS}, SPO <: AbstractSPSSOperator{SPBasis{SPSSBS}}}
    return operator.matrix_rep
end

# possibly recalculate the matrix representation
function recalculate!(operator :: SPLocalMSOperator{SPSSBS, SPMSB, SPO}, recursive::Bool=true, basis_change::Bool=true) where {SPSSBS <: AbstractSPSSBasisState, SPMSBS <: SPMSBasisState{SPSSBS}, SPMSB <: SPBasis{SPMSBS}, SPO <: AbstractSPSSOperator{SPBasis{SPSSBS}}}
    # maybe recursive
    if recursive
        # maybe change the basis of the SPSS operator
        if basis_change
            # get the single site basis
            operator.basis_ss = getSingleSiteBasis(basis(operator), operator.site)
            # set the single site basis in the operator of the operator
            operator.operator.basis = operator.basis_ss
            # get the indices of single multisite states that belong to single site basis
            #operator.indices_ss = Int64[i for i in 1:length(basis(operator)) if basis(operator)[i].site == operator.site]
            operator.indices_ss = zeros(Int64, length(operator.basis_ss)) .- 1
            for i in 1:length(operator.indices_ss)
                # look for sp ss state i among the states
                for j in 1:length(basis(operator))
                    if basis(operator)[j].site == operator.site && basis(operator)[j].state == operator.basis_ss[i]
                        operator.indices_ss[i] = j
                        break
                    end
                end
            end
        end
        # let operator recalculate
        recalculate!(operator.operator, recursive, basis_change)
    end
    # create new matrix
    operator.matrix_rep = zeros(Complex{Float64}, length(basis(operator)), length(basis(operator)))
    # get matrix representation of spss operator
    matrix_rep_ss = matrix_representation(operator.operator)
    # recalculate the matrix elements of the ms operator
    for a in 1:length(operator.basis_ss)
    for b in 1:length(operator.basis_ss)
        # get alpha and beta
        alpha = operator.indices_ss[a]
        beta  = operator.indices_ss[b]
        # set the respective entry
        operator.matrix_rep[alpha, beta] = matrix_rep_ss[a,b]
    end
    end
end

# set a parameter (returns (found parameter?, changed matrix?))
function set_parameter!(operator :: SPLocalMSOperator{SPSSBS, SPMSB, SPO}, parameter :: Symbol, value; print_result::Bool=false, recalculate::Bool=true, site::Union{Int64, Symbol}=-1, kwargs...) where {SPSSBS <: AbstractSPSSBasisState, SPMSBS <: SPMSBasisState{SPSSBS}, SPMSB <: SPBasis{SPMSBS}, SPO <: AbstractSPSSOperator{SPBasis{SPSSBS}}}
    # pass on to containing operator
    if site == operator.site || site == :all
        found_param, changed_matrix = set_parameter!(operator.operator, parameter, value, print_result=print_result, recalculate=recalculate; kwargs...)
        if recalculate && changed_matrix
            recalculate!(operator, false, false)
        end
        return (found_param, changed_matrix)
    else
        if print_result
            println("site $(site) is not matching to site of operator: site=$(operator.site)")
        end
        return (false, false)
    end
end

# get a parameter (returns (found parameter?, parameter value or nothing))
function get_parameter(operator :: SPLocalMSOperator{SPSSBS, SPMSB, SPO}, parameter :: Symbol; print_result::Bool=false, site::Union{Int64, Symbol}=-1, kwargs...) where {SPSSBS <: AbstractSPSSBasisState, SPMSBS <: SPMSBasisState{SPSSBS}, SPMSB <: SPBasis{SPMSBS}, SPO <: AbstractSPSSOperator{SPBasis{SPSSBS}}}
    # pass on to containing operator
    if site == operator.site || site == :all
        return get_parameter(operator.operator, parameter, print_result=print_result; kwargs...)
    else
        if print_result
            println("site $(site) is not matching to site of operator: site=$(operator.site)")
        end
        return nothing
    end
end

# get a list of parameters
function get_parameters(operator :: SPLocalMSOperator{SPSSBS, SPMSB, SPO}; site::Union{Int64, Symbol}=:all, kwargs...) where {SPSSBS <: AbstractSPSSBasisState, SPMSBS <: SPMSBasisState{SPSSBS}, SPMSB <: SPBasis{SPMSBS}, SPO <: AbstractSPSSOperator{SPBasis{SPSSBS}}}
    if site == operator.site || site == :all
        return get_parameters(operator.operator; kwargs...)
    else
        return Symbol[]
    end
end
