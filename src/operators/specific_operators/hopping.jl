# include the geometry (only temporary here)
#include("geometry_unitcell.jl")


# abstract super type for hopping operators
abstract type AbstractSPHoppingOperator{SPB} <: AbstractSPMSOperator{SPB} end

export AbstractSPHoppingOperator





# concrete hopping operator
"""
    mutable struct SPOrbitalHoppingOperator{
            SPMSB <: SPBasis{SPMSBS} where {SPMSBS<:SPMSBasisState{SPSSBS} where SPSSBS <: AbstractSPSSBasisState}
        } <: AbstractSPHoppingOperator{SPMSB}
The object defines the single particle orbital hopping operator.
# Fields
- `basis :: SPMSB`, the single-particle multi-site basis;
- `hopping_processes :: Vector{Tuple{Int64, Int64, Symbol}}`, the list of hopping orbitals;
- `hopping_strengths :: Dict{Symbol, Complex{Float64}}`, the hopping strengths;
- `matrix_rep :: Matrix{Complex{Float64}}`, the matrix representation.
"""
mutable struct SPOrbitalHoppingOperator{
        SPMSB <: SPBasis{SPMSBS} where {SPMSBS<:SPMSBasisState{SPSSBS} where SPSSBS <: AbstractSPSSBasisState}
    } <: AbstractSPHoppingOperator{SPMSB}

    # the MS basis
    basis :: SPMSB

    # the list of hopping orbitals
    hopping_processes :: Vector{Tuple{Int64, Int64, Symbol}}
    hopping_strengths :: Dict{Symbol, Complex{Float64}}

    # the matrix representation
    matrix_rep :: Matrix{Complex{Float64}}

    # custom constructor
    function SPOrbitalHoppingOperator(basis :: SPMSB) where {SPSSBS <: AbstractSPSSBasisState, SPMSBS <: SPMSBasisState{SPSSBS}, SPMSB <: SPBasis{SPMSBS}}
        # create a new operator
        op = new{SPMSB}(basis, Tuple{Int64, Int64, Symbol}[], Dict{Symbol, Complex{Float64}}(), zeros(Complex{Float64}, length(basis), length(basis)))
        # recalculate the matrix representation
        recalculate!(op)
        # return the operator
        return op
    end
end


function SPOrbitalHoppingOperator(basis :: SPMSB) where { 
    SPMSBS <: DelocalizedBasisStateXYZ, 
    SPMSB <: SPBasis{SPMSBS}}
    # create a new operator
    basis_internal = getMultiSiteBasis(getT2GBasisXYZ(),2)
    op = SPOrbitalHoppingOperator(basis_internal)
    # recalculate the matrix representation
    recalculate!(op)
    # build a projection operator around it
    op_proj = SPMSProjectorOperator(op, basis)
    # return the operator
    return op_proj
end



# export operator type
export SPOrbitalHoppingOperator

# creating a magnetic field operator on a multi site basis
function SPOrbitalHoppingOperator(basis::MPB) where {SPSSBS<:AbstractSPSSBasisState, N,MPB<:MPBasis{N,SPMSBasisState{SPSSBS}}}
    # construct new single site operator
    op = SPOrbitalHoppingOperator(basis.single_particle_basis)
    # construct new multi particle operator out of that
    return MPGeneralizedSPOperator(basis, op)
end





# show function
import Base.show
function Base.show(io::IO, op::SPOrbitalHoppingOperator{SPMSB}) where {SPSSBS <: AbstractSPSSBasisState, SPMSBS <: Union{SPMSBasisState{SPSSBS}, DelocalizedBasisStateXYZ}, SPMSB <: SPBasis{SPMSBS}}
    if haskey(io, :compact)
        print(io, "orbital hopping operator (", length(op.hopping_processes), " hopping elements, couplings ", keys(op.hopping_strengths), ") ")
    else
        print(io, "Orbital hopping operator\n")
        print(io, length(op.hopping_processes), " hopping elements, couplings ", string.(keys(op.hopping_strengths)), "\n")
        for h in op.hopping_processes
            print(io, "  - ", summary(basis(op)[h[1]]), " --> ", summary(basis(op)[h[2]]), ", :", string(h[3]), " = ", op.hopping_strengths[h[3]] , "\n")
        end
        print(io, "Multi-site basis contains "*string(length(basis(op)))*" states in total, defined on sites "*string(unique([b.site for b in basis(op)]))*"\n")
    end
end





################################################################################
#   Interface functions
################################################################################

# obtain the current basis
function basis(operator :: SPOrbitalHoppingOperator{SPMSB}) :: SPMSB where {SPSSBS <: AbstractSPSSBasisState, SPMSBS <: Union{SPMSBasisState{SPSSBS}, DelocalizedBasisStateXYZ}, SPMSB <: SPBasis{SPMSBS}}
    return operator.basis
end

# obtain the matrix representation
function matrix_representation(operator :: SPOrbitalHoppingOperator{SPMSB}) :: Matrix{Complex{Float64}} where {SPSSBS <: AbstractSPSSBasisState, SPMSBS <: Union{SPMSBasisState{SPSSBS}, DelocalizedBasisStateXYZ}, SPMSB <: SPBasis{SPMSBS}}
    return operator.matrix_rep
end

# possibly recalculate the matrix representation
function recalculate!(operator :: SPOrbitalHoppingOperator{SPMSB}, recursive::Bool=true, basis_change::Bool=true) where {SPSSBS <: AbstractSPSSBasisState, SPMSBS <: Union{SPMSBasisState{SPSSBS}, DelocalizedBasisStateXYZ}, SPMSB <: SPBasis{SPMSBS}}
    # create new matrix
    operator.matrix_rep = zeros(Complex{Float64}, length(basis(operator)), length(basis(operator)))
    # recalculate the matrix elements
    for h in operator.hopping_processes
        # obtain the strength
        t = operator.hopping_strengths[h[3]] :: Complex{Float64}
        # set the respective entries of the matrix
        operator.matrix_rep[h[2], h[1]] += t
        operator.matrix_rep[h[1], h[2]] += conj(t)
    end
end

# set a parameter (returns (found parameter?, changed matrix?))
function set_parameter!(operator :: SPOrbitalHoppingOperator{SPMSB}, parameter :: Symbol, value; print_result::Bool=false, recalculate::Bool=true, site::Union{Int64, Symbol}=-1, kwargs...) where {SPSSBS <: AbstractSPSSBasisState, SPMSBS <: Union{SPMSBasisState{SPSSBS}, DelocalizedBasisStateXYZ}, SPMSB <: SPBasis{SPMSBS}}
    if parameter in keys(operator.hopping_strengths)
        if operator.hopping_strengths[parameter] == value
            if print_result
                println("Parameter :$(parameter) found, but no changes as value $(value) still the same")
            end
            return (true, false)
        else
            operator.hopping_strengths[parameter] = value
            if recalculate
                recalculate!(operator, false, false)
                if print_result
                    println("Parameter :$(parameter) found, value changed to $(value), also operator recalculated")
                end
            else
                if print_result
                    println("Parameter :$(parameter) found, value changed to $(value)")
                end
            end
            return (true, true)
        end
    else
        if print_result
            println("Parameter :$(parameter) not found")
        end
        return (false, false)
    end
end

# get a parameter (returns (found parameter?, parameter value or nothing))
function get_parameter(operator :: SPOrbitalHoppingOperator{SPMSB}, parameter :: Symbol; print_result::Bool=false, site::Union{Int64, Symbol}=-1, kwargs...) where {SPSSBS <: AbstractSPSSBasisState, SPMSBS <: Union{SPMSBasisState{SPSSBS}, DelocalizedBasisStateXYZ}, SPMSB <: SPBasis{SPMSBS}}
    if parameter in keys(operator.hopping_strengths)
        if print_result
            println("Parameter :$(parameter) found and returned its value $(operator.hopping_strengths[parameter])")
        end
        return operator.hopping_strengths[parameter]
    else
        if print_result
            println("Parameter :$(parameter) not found")
        end
        return nothing
    end
end

# get a parameter (returns (found parameter?, parameter value or nothing))
function get_parameters(operator :: SPOrbitalHoppingOperator{SPMSB}; kwargs...) where {SPSSBS <: AbstractSPSSBasisState, SPMSBS <: Union{SPMSBasisState{SPSSBS}, DelocalizedBasisStateXYZ}, SPMSB <: SPBasis{SPMSBS}}
    return Symbol[k for k in keys(operator.hopping_strengths)]
end







################################################################################
#   hopping specifying functions
################################################################################

# add hopping from one orbital to another
function add_hopping!(operator :: SPOrbitalHoppingOperator{SPMSB}, state_1 :: Int64, state_2 :: Int64, coupling :: Symbol) where {SPSSBS <: AbstractSPSSBasisState, SPMSBS <: Union{SPMSBasisState{SPSSBS}}, SPMSB <: SPBasis{SPMSBS}}
    # push into the respective lists
    push!(operator.hopping_processes, (state_1, state_2, coupling))
    if !(coupling in keys(operator.hopping_strengths))
        operator.hopping_strengths[coupling] = 0.0
    end
    # return nothing
    return nothing
end
function add_hopping!(operator :: SPOrbitalHoppingOperator{SPMSB}, state_1 :: Int64, state_2 :: Int64, coupling :: Symbol, value :: Number) where {SPSSBS <: AbstractSPSSBasisState, SPMSBS <: Union{SPMSBasisState{SPSSBS}}, SPMSB <: SPBasis{SPMSBS}}
    # push into the respective lists
    add_hopping!(operator, state_1, state_2, coupling)
    operator.hopping_strengths[coupling] = value
    # return nothing
    return nothing
end

# add hopping from one orbital to another by specifying sites explicitly
function add_hopping!(operator :: SPOrbitalHoppingOperator{SPMSB}, state_1 :: Int64, site_1 :: Int64, state_2 :: Int64, site_2 :: Int64, coupling :: Symbol) where {SPSSBS <: AbstractSPSSBasisState, SPMSBS <: Union{SPMSBasisState{SPSSBS}}, SPMSB <: SPBasis{SPMSBS}}
    # find the respective index
    state_1_ms = 0
    for j in 1:length(basis(operator))
        if basis(operator)[j].site == site_1
            state_1 -= 1
        end
        if state_1 == 0
            state_1_ms = j
            break
        end
    end
    state_2_ms = 0
    for j in 1:length(basis(operator))
        if basis(operator)[j].site == site_2
            state_2 -= 1
        end
        if state_2 == 0
            state_2_ms = j
            break
        end
    end
    # add the hopping
    add_hopping!(operator, state_1_ms, state_2_ms, coupling)
end
function add_hopping!(operator :: SPOrbitalHoppingOperator{SPMSB}, state_1 :: Int64, site_1 :: Int64, state_2 :: Int64, site_2 :: Int64, coupling :: Symbol, value :: Number) where {SPSSBS <: AbstractSPSSBasisState, SPMSBS <: Union{SPMSBasisState{SPSSBS}}, SPMSB <: SPBasis{SPMSBS}}
    # push into the respective lists
    add_hopping!(operator, state_1, site_1, state_2, site_2, coupling)
    operator.hopping_strengths[coupling] = value
    # return nothing
    return nothing
end

# add hopping from one orbital to another by specifying strings of orbitals
function add_hopping!(operator :: SPOrbitalHoppingOperator{SPMSB}, state_1 :: String, state_2 :: String, coupling :: Symbol) where {SPSSBS <: AbstractSPSSBasisState, SPMSBS <: Union{SPMSBasisState{SPSSBS}}, SPMSB <: SPBasis{SPMSBS}}
    # find the respective indices
    state_1_ms = 0
    brackets = string(state_1[1]) * string(state_1[end])
    for j in 1:length(basis(operator))
        if summary(basis(operator)[j], brackets) == state_1
            state_1_ms = j
            break
        end
    end
    state_2_ms = 0
    brackets = string(state_1[1]) * string(state_1[end])
    for j in 1:length(basis(operator))
        if summary(basis(operator)[j], brackets) == state_2
            state_2_ms = j
            break
        end
    end
    # add the hopping
    add_hopping!(operator, state_1_ms, state_2_ms, coupling)
end
function add_hopping!(operator :: SPOrbitalHoppingOperator{SPMSB}, state_1 :: String, state_2 :: String, coupling :: Symbol, value :: Number) where {SPSSBS <: AbstractSPSSBasisState, SPMSBS <: Union{SPMSBasisState{SPSSBS}}, SPMSB <: SPBasis{SPMSBS}}
    # push into the respective lists
    add_hopping!(operator, state_1, state_2, coupling)
    operator.hopping_strengths[coupling] = value
    # return nothing
    return nothing
end
function add_diagonal_hopping!(operator :: SPOrbitalHoppingOperator{SPMSB}, site_1 :: Int64, site_2 :: Int64, coupling :: Symbol, value :: Number) where {SPSSBS <: AbstractSPSSBasisState, SPMSBS <: Union{SPMSBasisState{SPSSBS}}, SPMSB <: SPBasis{SPMSBS}}
    # find the number of site 1 and 2 states
    b1 = getSingleSiteBasis(basis(operator), site_1)
    b2 = getSingleSiteBasis(basis(operator), site_2)
    @assert length(b1) == length(b2)
    # push into the respective lists
    for i in 1:length(b1)
        add_hopping!(operator, i, site_1, i, site_2, coupling, value)
    end
    # return nothing
    return nothing
end

##############################################################
#   Convenience functions for SP MS hopping operators
##############################################################

# overwrite all add_hopping functions

# add_hopping! for DelocalizedBasisState

function add_hopping!(operator :: SPMSProjectorOperator{SPMSB_IN, SPMSB_OUT, SPO}, args...) where {SPSSBS <: AbstractSPSSBasisState, SPMSBS <: Union{SPMSBasisState{SPSSBS}, DelocalizedBasisStateXYZ}, SPMSB_IN<: SPBasis{SPMSBS}, SPMSB_OUT<: SPBasis{SPMSBS}, SPO <: AbstractSPOperator}
    # pass further into the hopping operator
    add_hopping!(operator.operator, args...)
    # return nothing
    return nothing
end
function add_diagonal_hopping!(operator :: SPMSProjectorOperator{SPMSB_IN, SPMSB_OUT, SPO}, args...) where {SPSSBS <: AbstractSPSSBasisState, SPMSBS <: Union{SPMSBasisState{SPSSBS}, DelocalizedBasisStateXYZ}, SPMSB_IN<: SPBasis{SPMSBS}, SPMSB_OUT<: SPBasis{SPMSBS}, SPO <: AbstractSPOperator}
    # pass further into the hopping operator
    add_diagonal_hopping!(operator.operator, args...)
    # return nothing
    return nothing
end

# add_hopping! for MPGeneralizedSPOperator

function add_hopping!(operator :: MPGeneralizedSPOperator{SPMSBS, MPB, SPOHO}, args...) where {SPSSBS <: AbstractSPSSBasisState, SPMSBS <: Union{SPMSBasisState{SPSSBS}}, SPMSB <: SPBasis{SPMSBS}, SPOHO<:SPOrbitalHoppingOperator{SPMSB}, N,MPB<:MPBasis{N,SPMSBasisState{SPSSBS}}}
    # pass further into the hopping operator
    add_hopping!(operator.operator, args...)
    # return nothing
    return nothing
end
function add_diagonal_hopping!(operator :: MPGeneralizedSPOperator{SPMSBS, MPB, SPOHO}, args...) where {SPSSBS <: AbstractSPSSBasisState, SPMSBS <: Union{SPMSBasisState{SPSSBS}}, SPMSB <: SPBasis{SPMSBS}, SPOHO<:SPOrbitalHoppingOperator{SPMSB}, N,MPB<:MPBasis{N,SPMSBasisState{SPSSBS}}}
    # pass further into the hopping operator
    add_diagonal_hopping!(operator.operator, args...)
    # return nothing
    return nothing
end

export add_hopping!, add_diagonal_hopping!
