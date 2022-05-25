#########################################################
####################  HOLE PICTURE  #####################
#########################################################

# define a generalized SP operator type for MP basis states
mutable struct MPHoleGeneralizedSPOperator{
        SPBS <: AbstractSPBasisState,
        MPB <: MPBasis{N,SPBS} where {N},
        SPO <: AbstractSPOperator{SPBasis{SPBS}}
    } <: AbstractMP1POperator{MPB}

    # the MP basis
    basis :: MPB
    # the single particle operator
    operator :: SPO
    # the matrix representation
    matrix_rep :: Matrix{Complex{Float64}}

    # custom constructor
    function MPHoleGeneralizedSPOperator(basis :: MPB, operator :: SPO) where {
                N,
                SPBS <: AbstractSPBasisState,
                MPB <: MPBasis{N,SPBS},
                SPO <: AbstractSPOperator{SPBasis{SPBS}}
            }
        # create a new operator
        op = new{SPBS, MPB, SPO}(basis, operator, zeros(length(basis), length(basis)))
        # recalculate the matrix representation
        recalculate!(op)
        # return the operator
        return op
    end
end
export MPHoleGeneralizedSPOperator




import Base.show
function Base.show(io::IO, op::MPHoleGeneralizedSPOperator{SPBS, MPB, SPO}) where {
            N,
            SPBS <: AbstractSPBasisState,
            MPB <: MPBasis{N,SPBS},
            SPO <: AbstractSPOperator{SPBasis{SPBS}}
        }
    if haskey(io, :compact)
        print(io, "1-particle ")
        show(io, op.operator)
    else
        print(io, "Multi-particle ")
        show(io, op.operator)
        print(io, "Multi-particle basis contains "*string(length(basis(op)))*" states in total, with "*string(N))
        printstyled(io, " holes", color=:light_red)
        print(io, " per state\n")    
    end
end




##############################################################
#   Interface functions
##############################################################

# obtain the current basis
function basis(operator :: MPHoleGeneralizedSPOperator{SPBS, MPB, SPO}) :: MPB where {
            N,
            SPBS <: AbstractSPBasisState,
            MPB <: MPBasis{N,SPBS},
            SPO <: AbstractSPOperator{SPBasis{SPBS}}
        }
    return operator.basis
end

# obtain the matrix representation
function matrix_representation(operator :: MPHoleGeneralizedSPOperator{SPBS, MPB, SPO}) :: Matrix{Complex{Float64}} where {
            N,
            SPBS <: AbstractSPBasisState,
            MPB <: MPBasis{N,SPBS},
            SPO <: AbstractSPOperator{SPBasis{SPBS}}
        }
    return operator.matrix_rep
end

# possibly recalculate the matrix representation
function recalculate!(operator :: MPHoleGeneralizedSPOperator{SPBS, MPB, SPO}, recursive::Bool=true, basis_change::Bool=true) where {
            N,
            SPBS <: AbstractSPBasisState,
            MPB <: MPBasis{N,SPBS},
            SPO <: AbstractSPOperator{SPBasis{SPBS}}
        }
    # maybe calculate recursively
    if recursive
        if basis_change
            # reset the basis in the single particle operator
            operator.operator.basis = operator.basis.single_particle_basis
            # let operator recalculate
            recalculate!(operator.operator, true, true)
        else
            # let operator recalculate
            recalculate!(operator.operator, true, false)
        end
    end
    # get matrix representation
    matrix_rep_sp = matrix_representation(operator.operator)
    # get the important matrix elements
    relevant_sp = map(x->abs(x)>1e-8, matrix_rep_sp)
    # create new matrix
    operator.matrix_rep = zeros(Complex{Float64}, length(basis(operator)), length(basis(operator)))
    # allocate a buffer state
    state_buffer = deepcopy(basis(operator)[1])
    state_buffer.basis_index = -1
    state_buffer.basis_sign  = 0
    # recalculate the own matrix elements
    for a in 1:length(operator.basis.single_particle_basis)
    for b in 1:length(operator.basis.single_particle_basis)
        # check if relevant
        if !relevant_sp[a,b]
            continue
        end
        # get the element of the single particle hamiltonian
        op_sp_ab = matrix_rep_sp[a,b]
        # generate all element contributions to the many body hamiltonian
        for alpha in basis(operator).lookup_sp_states[a]
        for beta  in basis(operator).lookup_sp_states[b]

            operator.matrix_rep[alpha, beta] += expectation_value_ca!(basis(operator), basis(operator)[alpha], a,b, basis(operator)[beta], state_buffer) * op_sp_ab

        end
        end
    end
    end
end

# set a parameter (returns (found parameter?, changed matrix?))
function set_parameter!(operator :: MPHoleGeneralizedSPOperator{SPBS, MPB, SPO}, parameter :: Symbol, value; print_result::Bool=true, recalculate::Bool=true, kwargs...) where {
            N,
            SPBS <: AbstractSPBasisState,
            MPB <: MPBasis{N,SPBS},
            SPO <: AbstractSPOperator{SPBasis{SPBS}}
        }
    # pass on to contained operator
    found_param, changed_matrix = set_parameter!(operator.operator, parameter, value, print_result=print_result, recalculate=recalculate; kwargs...)
    if recalculate && changed_matrix
        recalculate!(operator, false, false)
    end
    return (found_param, changed_matrix)
end

# get a parameter (returns (found parameter?, parameter value or nothing))
function get_parameter(operator :: MPHoleGeneralizedSPOperator{SPBS, MPB, SPO}, parameter :: Symbol; print_result::Bool=true, site::Union{Int64, Symbol}=-1, kwargs...) where {
            N,
            SPBS <: AbstractSPBasisState,
            MPB <: MPBasis{N,SPBS},
            SPO <: AbstractSPOperator{SPBasis{SPBS}}
        }
    # pass on to containing operator
    return get_parameter(operator.operator, parameter, print_result=print_result; site=site, kwargs...)
end

# get a parameter (returns (found parameter?, parameter value or nothing))
function get_parameters(operator :: MPHoleGeneralizedSPOperator{SPBS, MPB, SPO}; site::Union{Int64, Symbol}=-1, kwargs...) where {
            N,
            SPBS <: AbstractSPBasisState,
            MPB <: MPBasis{N,SPBS},
            SPO <: AbstractSPOperator{SPBasis{SPBS}}
        }
    # pass on to containing operator
    return get_parameters(operator.operator; site=site, kwargs...)
end


#########################################################
##################  ELECTRON PICTURE  ###################
#########################################################



# define a generalized SP operator type for MP basis states
mutable struct MPElectronGeneralizedSPOperator{
        SPBS <: AbstractSPBasisState,
        MPB <: MPBasis{N,SPBS} where {N},
        SPO <: AbstractSPOperator{SPBasis{SPBS}}
    } <: AbstractMP1POperator{MPB}

    # the MP basis
    basis :: MPB
    # the single particle operator
    operator :: SPO
    # the matrix representation
    matrix_rep :: Matrix{Complex{Float64}}

    # custom constructor
    function MPElectronGeneralizedSPOperator(basis :: MPB, operator :: SPO) where {
                N,
                SPBS <: AbstractSPBasisState,
                MPB <: MPBasis{N,SPBS},
                SPO <: AbstractSPOperator{SPBasis{SPBS}}
            }
        # create a new operator
        op = new{SPBS, MPB, SPO}(basis, operator, zeros(length(basis), length(basis)))
        # recalculate the matrix representation
        recalculate!(op)
        # return the operator
        return op
    end
end
export MPElectronGeneralizedSPOperator




import Base.show
function Base.show(io::IO, op::MPElectronGeneralizedSPOperator{SPBS, MPB, SPO}) where {
            N,
            SPBS <: AbstractSPBasisState,
            MPB <: MPBasis{N,SPBS},
            SPO <: AbstractSPOperator{SPBasis{SPBS}}
        }
    if haskey(io, :compact)
        print(io, "1-particle ")
        show(io, op.operator)
    else
        print(io, "Multi-particle ")
        show(io, op.operator)
        print(io, "Multi-particle basis contains "*string(length(basis(op)))*" states in total, with "*string(N))
        printstyled(io, " electrons", color=:light_blue )
        print(io, " per state\n")
    end
end




##############################################################
#   Interface functions
##############################################################

# obtain the current basis
function basis(operator :: MPElectronGeneralizedSPOperator{SPBS, MPB, SPO}) :: MPB where {
            N,
            SPBS <: AbstractSPBasisState,
            MPB <: MPBasis{N,SPBS},
            SPO <: AbstractSPOperator{SPBasis{SPBS}}
        }
    return operator.basis
end

# obtain the matrix representation
function matrix_representation(operator :: MPElectronGeneralizedSPOperator{SPBS, MPB, SPO}) :: Matrix{Complex{Float64}} where {
            N,
            SPBS <: AbstractSPBasisState,
            MPB <: MPBasis{N,SPBS},
            SPO <: AbstractSPOperator{SPBasis{SPBS}}
        }
    return operator.matrix_rep
end

# possibly recalculate the matrix representation
function recalculate!(operator :: MPElectronGeneralizedSPOperator{SPBS, MPB, SPO}, recursive::Bool=true, basis_change::Bool=true) where {
            N,
            SPBS <: AbstractSPBasisState,
            MPB <: MPBasis{N,SPBS},
            SPO <: AbstractSPOperator{SPBasis{SPBS}}
        }
    # maybe calculate recursively
    if recursive
        if basis_change
            # reset the basis in the single particle operator
            operator.operator.basis = operator.basis.single_particle_basis
            # let operator recalculate
            recalculate!(operator.operator, true, true)
        else
            # let operator recalculate
            recalculate!(operator.operator, true, false)
        end
    end
    # get matrix representation
    matrix_rep_sp = matrix_representation(operator.operator)
    # get the important matrix elements
    relevant_sp = map(x->abs(x)>1e-8, matrix_rep_sp)
    # create new matrix
    operator.matrix_rep = zeros(Complex{Float64}, length(basis(operator)), length(basis(operator)))
    # allocate a buffer state
    state_buffer = deepcopy(basis(operator)[1])
    state_buffer.basis_index = -1
    state_buffer.basis_sign  = 0
    # recalculate the own matrix elements
    for a in 1:length(operator.basis.single_particle_basis)
    for b in 1:length(operator.basis.single_particle_basis)
        # check if relevant; if Hab_sp=0, skip cycle
        if !relevant_sp[a,b]
            continue
        end
        # get the element of the single particle hamiltonian
        op_sp_ab = matrix_rep_sp[a,b]
        # generate all element contributions to the many body hamiltonian
        
        for alpha in 1:length(basis(operator))#.lookup_sp_states[a]
        for beta  in 1:length(basis(operator)) #basis(operator).lookup_sp_states[b]
            
            operator.matrix_rep[alpha, beta] += ( delta(a,b)*overlap(basis(operator),basis(operator)[alpha],basis(operator),basis(operator)[beta])- expectation_value_ca!(basis(operator), basis(operator)[alpha], b,a, basis(operator)[beta], state_buffer)) * op_sp_ab
        end
        end
    end
    end
    
end

# set a parameter (returns (found parameter?, changed matrix?))
function set_parameter!(operator :: MPElectronGeneralizedSPOperator{SPBS, MPB, SPO}, parameter :: Symbol, value; print_result::Bool=true, recalculate::Bool=true, kwargs...) where {
            N,
            SPBS <: AbstractSPBasisState,
            MPB <: MPBasis{N,SPBS},
            SPO <: AbstractSPOperator{SPBasis{SPBS}}
        }
    # pass on to contained operator
    found_param, changed_matrix = set_parameter!(operator.operator, parameter, value, print_result=print_result, recalculate=recalculate; kwargs...)
    if recalculate && changed_matrix
        recalculate!(operator, false, false)
    end
    return (found_param, changed_matrix)
end

# get a parameter (returns (found parameter?, parameter value or nothing))
function get_parameter(operator :: MPElectronGeneralizedSPOperator{SPBS, MPB, SPO}, parameter :: Symbol; print_result::Bool=true, site::Union{Int64, Symbol}=-1, kwargs...) where {
            N,
            SPBS <: AbstractSPBasisState,
            MPB <: MPBasis{N,SPBS},
            SPO <: AbstractSPOperator{SPBasis{SPBS}}
        }
    # pass on to containing operator
    return get_parameter(operator.operator, parameter, print_result=print_result; site=site, kwargs...)
end

# get a parameter (returns (found parameter?, parameter value or nothing))
function get_parameters(operator :: MPElectronGeneralizedSPOperator{SPBS, MPB, SPO}; site::Union{Int64, Symbol}=-1, kwargs...) where {
            N,
            SPBS <: AbstractSPBasisState,
            MPB <: MPBasis{N,SPBS},
            SPO <: AbstractSPOperator{SPBasis{SPBS}}
        }
    # pass on to containing operator
    return get_parameters(operator.operator; site=site, kwargs...)
end
