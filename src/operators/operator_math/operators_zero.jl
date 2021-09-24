################################################################################
#   Type definition
################################################################################


# Type Definition of zero Operator
mutable struct ZeroOperator{B} <: AbstractOperator{B}
    # the basis
    basis :: B
    # the current matrix representation
    matrix_rep :: Matrix{Complex{Float64}}

    # Custom constructor (without explicit matrix rep)
    function ZeroOperator(b :: B) where {BS<:AbstractBasisState, B<:AbstractBasis{BS}}
        # construct new operator
        op = new{B}(b, zeros(Complex{Float64}, length(b), length(b)))
        # return the operator
        return op
    end
end

# export operator type
export  ZeroOperator


function Base.show(io::IO, op::ZeroOperator{B}) where {BS<:AbstractBasisState, B<:AbstractBasis{BS}}
    if haskey(io, :compact)
        print(io, "zero operator (0)")
    else
        print(io, "zero operator (0)\n")
        printBasisInformation(io, op)
    end
end

function printBasisInformation(io::IO, op::ZeroOperator{B}) where {BS<:AbstractBasisState, B<:AbstractBasis{BS}}
    print(io, "Basis contains "*string(length(basis(op)))*" states in total, states are of type "*string(BS)*"\n")
end
function printBasisInformation(io::IO, op::ZeroOperator{B}) where {BS<:AbstractSPSSBasisState, B<:AbstractBasis{BS}}
    print(io, "Basis is single-particle / single-site basis and contains "*string(length(basis(op)))*" states of type "*string(BS)*"\n")
end
function printBasisInformation(io::IO, op::ZeroOperator{B}) where {SPBS<:AbstractSPSSBasisState, BS<:SPMSBasisState{SPBS}, B<:AbstractBasis{BS}}
    print(io, "Basis is single-particle / multi-site basis and contains "*string(length(basis(op)))*" states with SPSS type "*string(SPBS)*"\n")
    print(io, "Basis contains sites "*string(unique([b.site for b in basis(op)]))*"\n")
end
function printBasisInformation(io::IO, op::ZeroOperator{B}) where {N, MPBS<:MPBasisState{N}, SPBS<:AbstractSPSSBasisState, B<:MPBasis{MPBS,SPBS}}
    print(io, "Basis is multi-particle / single-site basis and contains "*string(length(basis(op)))*" states with "*string(N)*" particles each\n")
    print(io, "SPSS basis contains "*string(length(basis(op.single_particle_basis)))*" states of type "*string(SPBS)*"\n")
end
function printBasisInformation(io::IO, op::ZeroOperator{B}) where {N, MPBS<:MPBasisState{N}, SPSSBS<:AbstractSPSSBasisState, SPBS<:SPMSBasisState{SPSSBS}, B<:MPBasis{N,SPBS}}
    print(io, "Basis is multi-particle / multi-site basis and contains "*string(length(basis(op)))*" states with "*string(N)*" particles each\n")
    print(io, "SPSS basis contains "*string(length(basis(op).single_particle_basis))*" states of type "*string(SPSSBS)*" on sites "*string(unique([b.site for b in basis(op).single_particle_basis]))*"\n")
end



################################################################################
#   Interface functions
################################################################################

# obtain the current basis
function basis(operator :: ZeroOperator{B}) :: B where {BS<:AbstractBasisState, B<:AbstractBasis{BS}}
    return operator.basis
end

# obtain the matrix representation
function matrix_representation(operator :: ZeroOperator{B}) :: Matrix{Complex{Float64}} where {BS<:AbstractBasisState, B<:AbstractBasis{BS}}
    return operator.matrix_rep
end

# possibly recalculate the matrix representation
function recalculate!(operator :: ZeroOperator{B}, recursive::Bool=true, basis_change::Bool=true)  where {BS<:AbstractBasisState, B<:AbstractBasis{BS}}
    # create new matrix
    operator.matrix_rep = zeros(Complex{Float64}, length(basis(operator)), length(basis(operator)))
end

# set a parameter (returns (found parameter?, changed matrix?))
function set_parameter!(operator :: ZeroOperator{B}, parameter :: Symbol, value; print_result::Bool=false, recalculate::Bool=true, kwargs...) where {BS<:AbstractBasisState, B<:AbstractBasis{BS}}
    if print_result
        println("Parameter :$(parameter) not found")
    end
    return (false, false)
end

# get a parameter (returns (found parameter?, parameter value or nothing))
function get_parameter(operator :: ZeroOperator{B}, parameter :: Symbol; print_result::Bool=false, kwargs...) where {BS<:AbstractBasisState, B<:AbstractBasis{BS}}
    if print_result
        println("Parameter :$(parameter) not found")
    end
    return nothing
end

# get a parameter (returns (found parameter?, parameter value or nothing))
function get_parameters(operator :: ZeroOperator{B}; kwargs...) where {BS<:AbstractBasisState, B<:AbstractBasis{BS}}
    return Symbol[]
end
