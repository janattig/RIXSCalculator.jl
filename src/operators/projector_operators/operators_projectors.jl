# Single particle projectors

################################################################################
#   Type definition
################################################################################

# Type definition of local MS operators (SPSS operator only acting on a single site of many)
"""
    mutable struct SPSSProjectorOperator{
            SPSSBS_IN  <: AbstractSPSSBasisState,
            SPSSBS_OUT <: AbstractSPSSBasisState,
            SPSSB_IN   <: SPBasis{SPSSBS_IN},
            SPSSB_OUT  <: SPBasis{SPSSBS_OUT},
            SPO <: AbstractSPSSOperator{SPSSB_IN}
        } <: AbstractSPSSOperator{SPBasis{SPSSBS_OUT}}

This object defines the local multi-site projector operator. This operator only acts on a single site of many.

# Fields

- `basis_in  :: SPSSB_IN`, `basis :: SPSSB_OUT`, the multi-site bases
- `operator  :: SPO`, the single site basis;
- `projector_in_out :: Matrix{Complex{Float64}}`, the matrix representation of the projector;
- `matrix_rep :: Matrix{Complex{Float64}}`, the matrix representation.

"""
mutable struct SPSSProjectorOperator{
        SPSSBS_IN  <: AbstractSPSSBasisState,
        SPSSBS_OUT <: AbstractSPSSBasisState,
        SPSSB_IN   <: SPBasis{SPSSBS_IN},
        SPSSB_OUT  <: SPBasis{SPSSBS_OUT},
        SPO <: AbstractSPSSOperator{SPSSB_IN}
    } <: AbstractSPSSOperator{SPBasis{SPSSBS_OUT}}

    # the MS basis
    basis_in  :: SPSSB_IN
    basis :: SPSSB_OUT
    # the single site operator
    operator  :: SPO
    # the matrix representation of the projector
    projector_in_out :: Matrix{Complex{Float64}}
    # the matrix representation
    matrix_rep :: Matrix{Complex{Float64}}

    # custom constructor
    function SPSSProjectorOperator(operator :: SPO, basis_to :: SPSSB_OUT) where {
                SPSSBS_IN  <: AbstractSPSSBasisState,
                SPSSBS_OUT <: AbstractSPSSBasisState,
                SPSSB_IN   <: SPBasis{SPSSBS_IN},
                SPSSB_OUT  <: SPBasis{SPSSBS_OUT},
                SPO <: AbstractSPSSOperator{SPSSB_IN}
            }
        # create a new operator
        op = new{SPSSBS_IN,SPSSBS_OUT, SPSSB_IN,SPSSB_OUT, SPO}(basis(operator), basis_to, operator, projector_matrix(basis(operator), basis_to), zeros(length(basis_to), length(basis_to)))
        # recalculate the matrix representation
        recalculate!(op, true)
        # return the operator
        return op
    end
end

# export operator type
export  SPSSProjectorOperator



function ProjectorOperator(operator :: SPO, basis_to :: SPSSB_OUT) where {
            SPSSBS_IN  <: AbstractSPSSBasisState,
            SPSSBS_OUT <: AbstractSPSSBasisState,
            SPSSB_IN   <: SPBasis{SPSSBS_IN},
            SPSSB_OUT  <: SPBasis{SPSSBS_OUT},
            SPO <: AbstractSPSSOperator{SPSSB_IN}
        }

    return SPSSProjectorOperator(operator, basis_to)
end

export ProjectorOperator



# overwritten show function
import Base.show
function Base.show(io::IO, op::OP) where {
            SPSSBS_IN  <: AbstractSPSSBasisState,
            SPSSBS_OUT <: AbstractSPSSBasisState,
            SPSSB_IN   <: SPBasis{SPSSBS_IN},
            SPSSB_OUT  <: SPBasis{SPSSBS_OUT},
            SPO <: AbstractSPSSOperator{SPSSB_IN},
            OP <: SPSSProjectorOperator{SPSSBS_IN,SPSSBS_OUT, SPSSB_IN,SPSSB_OUT, SPO}
        }
    if haskey(io, :compact)
        print(io, "(projected) ")
        show(io, op.operator)
    else
        print(io, "Projection to new basis of operator:\n")
        show(io, op.operator)
        print(io, "New (outer basis) consists of "*string(length(basis(op)))*" states of type $(SPSSBS_OUT)\n")
    end
end


################################################################################
#   Interface functions
################################################################################

# obtain the current basis
function basis(operator :: SPSSProjectorOperator{SPSSBS_IN,SPSSBS_OUT, SPSSB_IN,SPSSB_OUT, SPO}) :: SPSSB_OUT  where {
            SPSSBS_IN  <: AbstractSPSSBasisState,
            SPSSBS_OUT <: AbstractSPSSBasisState,
            SPSSB_IN   <: SPBasis{SPSSBS_IN},
            SPSSB_OUT  <: SPBasis{SPSSBS_OUT},
            SPO <: AbstractSPSSOperator{SPSSB_IN}
        }
    return operator.basis
end

# obtain the matrix representation
function matrix_representation(operator :: SPSSProjectorOperator{SPSSBS_IN,SPSSBS_OUT, SPSSB_IN,SPSSB_OUT, SPO}) :: Matrix{Complex{Float64}}  where {
            SPSSBS_IN  <: AbstractSPSSBasisState,
            SPSSBS_OUT <: AbstractSPSSBasisState,
            SPSSB_IN   <: SPBasis{SPSSBS_IN},
            SPSSB_OUT  <: SPBasis{SPSSBS_OUT},
            SPO <: AbstractSPSSOperator{SPSSB_IN}
        }
    return operator.matrix_rep
end

# possibly recalculate the matrix representation
function recalculate!(operator :: SPSSProjectorOperator{SPSSBS_IN,SPSSBS_OUT, SPSSB_IN,SPSSB_OUT, SPO}, recursive::Bool=true, basis_change::Bool=true) where {
            SPSSBS_IN  <: AbstractSPSSBasisState,
            SPSSBS_OUT <: AbstractSPSSBasisState,
            SPSSB_IN   <: SPBasis{SPSSBS_IN},
            SPSSB_OUT  <: SPBasis{SPSSBS_OUT},
            SPO <: AbstractSPSSOperator{SPSSB_IN}
        }
    # maybe recalculate the inner operator
    if recursive
        recalculate!(operator.operator, recursive, basis_change)
    end
    # calculate new projector P_12
    if basis_change
        operator.projector_in_out = projector_matrix(operator.basis_in, operator.basis)
    end
    # calculate new matrix h_2 = P_12' * h_1 * P_12
    operator.matrix_rep = operator.projector_in_out' * matrix_representation(operator.operator) * operator.projector_in_out
    # return nothing
    return nothing
end

# set a parameter (returns (found parameter?, changed matrix?))
function set_parameter!(operator :: SPSSProjectorOperator{SPSSBS_IN,SPSSBS_OUT, SPSSB_IN,SPSSB_OUT, SPO}, parameter :: Symbol, value; print_result::Bool=false, recalculate::Bool=true, kwargs...) where {
            SPSSBS_IN  <: AbstractSPSSBasisState,
            SPSSBS_OUT <: AbstractSPSSBasisState,
            SPSSB_IN   <: SPBasis{SPSSBS_IN},
            SPSSB_OUT  <: SPBasis{SPSSBS_OUT},
            SPO <: AbstractSPSSOperator{SPSSB_IN}
        }
    # pass to inner operator
    (fp, cm) = set_parameter!(operator.operator, parameter, value; print_result=print_result, recalculate=recalculate, kwargs...)
    # maybe recalculate itself
    if recalculate && fp && cm
        recalculate!(operator, false, false)
    end
    # return the same variables
    return (fp, cm)
end

# get a parameter (returns (parameter value or nothing))
function get_parameter(operator :: SPSSProjectorOperator{SPSSBS_IN,SPSSBS_OUT, SPSSB_IN,SPSSB_OUT, SPO}, parameter :: Symbol; print_result::Bool=false, kwargs...) where {
            SPSSBS_IN  <: AbstractSPSSBasisState,
            SPSSBS_OUT <: AbstractSPSSBasisState,
            SPSSB_IN   <: SPBasis{SPSSBS_IN},
            SPSSB_OUT  <: SPBasis{SPSSBS_OUT},
            SPO <: AbstractSPSSOperator{SPSSB_IN}
        }
    # check if parameter can be set
    return get_parameter(operator.operator, parameter; print_result=print_result, kwargs...)
end

# get a list of parameters
function get_parameters(operator :: SPSSProjectorOperator{SPSSBS_IN,SPSSBS_OUT, SPSSB_IN,SPSSB_OUT, SPO}; kwargs...) where {
            SPSSBS_IN  <: AbstractSPSSBasisState,
            SPSSBS_OUT <: AbstractSPSSBasisState,
            SPSSB_IN   <: SPBasis{SPSSBS_IN},
            SPSSB_OUT  <: SPBasis{SPSSBS_OUT},
            SPO <: AbstractSPSSOperator{SPSSB_IN}
        }
    return get_parameters(operator.operator; kwargs...)
end
















# Single particle multi site projectors

################################################################################
#   Type definition
################################################################################

# Type definition of local MS operators (SPSS operator only acting on a single site of many)
"""
    mutable struct SPMSProjectorOperator{
            SPMSB_IN   <: SPBasis{SPMSBS_IN} where {SPMSBS_IN<:Union{SPMSBasisState{BS} where BS, SPMSCompositeBasisState{B} where B, DelocalizedBasisState{BS} where BS}},
            SPMSB_OUT  <: SPBasis{SPMSBS_OUT} where {SPMSBS_OUT<:Union{SPMSBasisState{BS} where BS, SPMSCompositeBasisState{B} where B, DelocalizedBasisState{BS} where BS}},
            SPO <: AbstractSPMSOperator{SPMSB_IN}
        } <: AbstractSPMSOperator{SPMSB_OUT}

This object defines the local multi-site operator.

# Fields

- `basis_in  :: SPMSB_IN`, `basis :: SPMSB_OUT`, the multi-site bases;
- `operator  :: SPO`, the single particle multi site operator;
- `projector_in_out :: Matrix{Complex{Float64}}`, the matrix representation of the projector;
- `matrix_rep :: Matrix{Complex{Float64}}`, the matrix representation.

"""
mutable struct SPMSProjectorOperator{
        SPMSB_IN   <: SPBasis{SPMSBS_IN} where {SPMSBS_IN<:Union{SPMSBasisState{BS} where BS, SPMSCompositeBasisState{B} where B, DelocalizedBasisState{BS} where BS}},
        SPMSB_OUT  <: SPBasis{SPMSBS_OUT} where {SPMSBS_OUT<:Union{SPMSBasisState{BS} where BS, SPMSCompositeBasisState{B} where B, DelocalizedBasisState{BS} where BS}},
        SPO <: AbstractSPMSOperator{SPMSB_IN}
    } <: AbstractSPMSOperator{SPMSB_OUT}

    # the MS basis
    basis_in  :: SPMSB_IN
    basis :: SPMSB_OUT
    # the single particle multi site operator
    operator  :: SPO
    # the matrix representation of the projector
    projector_in_out :: Matrix{Complex{Float64}}
    # the matrix representation
    matrix_rep :: Matrix{Complex{Float64}}

    # custom constructor
    function SPMSProjectorOperator(operator :: SPO, basis_to :: SPMSB_OUT) where {
                SPMSB_IN   <: SPBasis{SPMSBS_IN} where {SPMSBS_IN<:Union{SPMSBasisState{BS} where BS, SPMSCompositeBasisState{B} where B, DelocalizedBasisState{BS} where BS}},
                SPMSB_OUT  <: SPBasis{SPMSBS_OUT} where {SPMSBS_OUT<:Union{SPMSBasisState{BS} where BS, SPMSCompositeBasisState{B} where B, DelocalizedBasisState{BS} where BS}},
                SPO <: AbstractSPMSOperator{SPMSB_IN}
            }
        # create a new operator
        op = new{SPMSB_IN,SPMSB_OUT, SPO}(basis(operator), basis_to, operator, projector_matrix(basis(operator), basis_to), zeros(length(basis_to), length(basis_to)))
        # recalculate the matrix representation
        recalculate!(op, true)
        # return the operator
        return op
    end
end

# export operator type
export  SPMSProjectorOperator



function ProjectorOperator(operator :: SPO, basis_to :: SPMSB_OUT) where {
            SPMSB_IN   <: SPBasis{SPMSBS_IN} where {SPMSBS_IN<:Union{SPMSBasisState{BS} where BS, SPMSCompositeBasisState{B} where B, DelocalizedBasisState{BS} where BS}},
            SPMSB_OUT  <: SPBasis{SPMSBS_OUT} where {SPMSBS_OUT<:Union{SPMSBasisState{BS} where BS, SPMSCompositeBasisState{B} where B, DelocalizedBasisState{BS} where BS}},
            SPO <: AbstractSPMSOperator{SPMSB_IN}
        }

    return SPMSProjectorOperator(operator, basis_to)
end

export ProjectorOperator



# overwritten show function
import Base.show
function Base.show(io::IO, op::OP) where {
            SPMSB_IN   <: SPBasis{SPMSBS_IN} where {SPMSBS_IN<:Union{SPMSBasisState{BS} where BS, SPMSCompositeBasisState{B} where B, DelocalizedBasisState{BS} where BS}},
            SPMSB_OUT  <: SPBasis{SPMSBS_OUT} where {SPMSBS_OUT<:Union{SPMSBasisState{BS} where BS, SPMSCompositeBasisState{B} where B, DelocalizedBasisState{BS} where BS}},
            SPO <: AbstractSPMSOperator{SPMSB_IN},
            OP <: SPMSProjectorOperator{SPMSB_IN,SPMSB_OUT, SPO}
        }
    if haskey(io, :compact)
        print(io, "(projected) ")
        show(io, op.operator)
    else
        print(io, "Projection to new basis of operator:\n")
        show(io, op.operator)
        print(io, "New (outer basis) consists of "*string(length(basis(op)))*" states of type $(SPMSBS_OUT)\n")
    end
end


################################################################################
#   Interface functions
################################################################################

# obtain the current basis
function basis(operator :: SPMSProjectorOperator{SPMSB_IN,SPMSB_OUT, SPO}) :: SPMSB_OUT  where {
            SPMSB_IN   <: SPBasis{SPMSBS_IN} where {SPMSBS_IN<:Union{SPMSBasisState{BS} where BS, SPMSCompositeBasisState{B} where B, DelocalizedBasisState{BS} where BS}},
            SPMSB_OUT  <: SPBasis{SPMSBS_OUT} where {SPMSBS_OUT<:Union{SPMSBasisState{BS} where BS, SPMSCompositeBasisState{B} where B, DelocalizedBasisState{BS} where BS}},
            SPO <: AbstractSPMSOperator{SPMSB_IN}
        }
    return operator.basis
end

# obtain the matrix representation
function matrix_representation(operator :: SPMSProjectorOperator{SPMSB_IN,SPMSB_OUT, SPO}) :: Matrix{Complex{Float64}}  where {
            SPMSB_IN   <: SPBasis{SPMSBS_IN} where {SPMSBS_IN<:Union{SPMSBasisState{BS} where BS, SPMSCompositeBasisState{B} where B, DelocalizedBasisState{BS} where BS}},
            SPMSB_OUT  <: SPBasis{SPMSBS_OUT} where {SPMSBS_OUT<:Union{SPMSBasisState{BS} where BS, SPMSCompositeBasisState{B} where B, DelocalizedBasisState{BS} where BS}},
            SPO <: AbstractSPMSOperator{SPMSB_IN}
        }
    return operator.matrix_rep
end

# possibly recalculate the matrix representation
function recalculate!(operator :: SPMSProjectorOperator{SPMSB_IN,SPMSB_OUT, SPO}, recursive::Bool=true, basis_change::Bool=true) where {
            SPMSB_IN   <: SPBasis{SPMSBS_IN} where {SPMSBS_IN<:Union{SPMSBasisState{BS} where BS, SPMSCompositeBasisState{B} where B, DelocalizedBasisState{BS} where BS}},
            SPMSB_OUT  <: SPBasis{SPMSBS_OUT} where {SPMSBS_OUT<:Union{SPMSBasisState{BS} where BS, SPMSCompositeBasisState{B} where B, DelocalizedBasisState{BS} where BS}},
            SPO <: AbstractSPMSOperator{SPMSB_IN}
        }
    # maybe recalculate the inner operator
    if recursive
        recalculate!(operator.operator, recursive, basis_change)
    end
    # calculate new projector P_12
    if basis_change
        operator.projector_in_out = projector_matrix(operator.basis_in, operator.basis)
    end
    # calculate new matrix h_2 = P_12' * h_1 * P_12
    operator.matrix_rep = operator.projector_in_out' * matrix_representation(operator.operator) * operator.projector_in_out
    # return nothing
    return nothing
end

# set a parameter (returns (found parameter?, changed matrix?))
function set_parameter!(operator :: SPMSProjectorOperator{SPMSB_IN,SPMSB_OUT, SPO}, parameter :: Symbol, value; print_result::Bool=false, recalculate::Bool=true, kwargs...) where {
            SPMSB_IN   <: SPBasis{SPMSBS_IN} where {SPMSBS_IN<:Union{SPMSBasisState{BS} where BS, SPMSCompositeBasisState{B} where B, DelocalizedBasisState{BS} where BS}},
            SPMSB_OUT  <: SPBasis{SPMSBS_OUT} where {SPMSBS_OUT<:Union{SPMSBasisState{BS} where BS, SPMSCompositeBasisState{B} where B, DelocalizedBasisState{BS} where BS}},
            SPO <: AbstractSPMSOperator{SPMSB_IN}
        }
    # pass to inner operator
    (fp, cm) = set_parameter!(operator.operator, parameter, value; print_result=print_result, recalculate=recalculate, kwargs...)
    # maybe recalculate itself
    if recalculate && fp && cm
        recalculate!(operator, false, false)
    end
    # return the same variables
    return (fp, cm)
end

# get a parameter (returns (parameter value or nothing))
function get_parameter(operator :: SPMSProjectorOperator{SPMSB_IN,SPMSB_OUT, SPO}, parameter :: Symbol; print_result::Bool=false, kwargs...) where {
            SPMSB_IN   <: SPBasis{SPMSBS_IN} where {SPMSBS_IN<:Union{SPMSBasisState{BS} where BS, SPMSCompositeBasisState{B} where B, DelocalizedBasisState{BS} where BS}},
            SPMSB_OUT  <: SPBasis{SPMSBS_OUT} where {SPMSBS_OUT<:Union{SPMSBasisState{BS} where BS, SPMSCompositeBasisState{B} where B, DelocalizedBasisState{BS} where BS}},
            SPO <: AbstractSPMSOperator{SPMSB_IN}
        }
    # check if parameter can be set
    return get_parameter(operator.operator, parameter; print_result=print_result, kwargs...)
end

# get a list of parameters
function get_parameters(operator :: SPMSProjectorOperator{SPMSB_IN,SPMSB_OUT, SPO}; kwargs...) where {
            SPMSB_IN   <: SPBasis{SPMSBS_IN} where {SPMSBS_IN<:Union{SPMSBasisState{BS} where BS, SPMSCompositeBasisState{B} where B, DelocalizedBasisState{BS} where BS}},
            SPMSB_OUT  <: SPBasis{SPMSBS_OUT} where {SPMSBS_OUT<:Union{SPMSBasisState{BS} where BS, SPMSCompositeBasisState{B} where B, DelocalizedBasisState{BS} where BS}},
            SPO <: AbstractSPMSOperator{SPMSB_IN}
        }
    return get_parameters(operator.operator; kwargs...)
end






















# Single particle multi site projectors

################################################################################
#   Type definition
################################################################################

# Type definition of local MS operators (SPSS operator only acting on a single site of many)
mutable struct MPProjectorOperator{
        N, NO,
        SPBS_IN  <: AbstractSPBasisState,
        SPBS_OUT <: AbstractSPBasisState,
        MPB_IN  <: MPBasis{N,SPBS_IN},
        MPB_OUT <: MPBasis{N,SPBS_OUT},
        MPO <: AbstractMPOperator{NO,MPB_IN}
    } <: AbstractMPOperator{NO, MPB_OUT}

    # the MS basis
    basis_in  :: MPB_IN
    basis     :: MPB_OUT
    # the multi particle operator
    operator  :: MPO
    # the matrix representation of the projector
    projector_in_out :: Matrix{Complex{Float64}}
    # the matrix representation
    matrix_rep :: Matrix{Complex{Float64}}

    # custom constructor
    function MPProjectorOperator(operator :: MPO, basis_to :: MPB_OUT) where {
                N, NO,
                SPBS_IN  <: AbstractSPBasisState,
                SPBS_OUT <: AbstractSPBasisState,
                MPB_IN  <: MPBasis{N,SPBS_IN},
                MPB_OUT <: MPBasis{N,SPBS_OUT},
                MPO <: AbstractMPOperator{NO,MPB_IN}
            }
        # create a new operator
        op = new{N,NO, SPBS_IN,SPBS_OUT, MPB_IN,MPB_OUT, MPO}(basis(operator), basis_to, operator, projector_matrix(basis(operator), basis_to), zeros(length(basis_to), length(basis_to)))
        # recalculate the matrix representation
        recalculate!(op, true)
        # return the operator
        return op
    end
end

# export operator type
export  MPProjectorOperator



function ProjectorOperator(operator :: MPO, basis_to :: MPB_OUT) where {
            N, NO,
            SPBS_IN  <: AbstractSPBasisState,
            SPBS_OUT <: AbstractSPBasisState,
            MPB_IN  <: MPBasis{N,SPBS_IN},
            MPB_OUT <: MPBasis{N,SPBS_OUT},
            MPO <: AbstractMPOperator{NO,MPB_IN}
        }

    return MPProjectorOperator(operator, basis_to)
end

export ProjectorOperator



# overwritten show function
function Base.show(io::IO, op::OP) where {
            N, NO,
            SPBS_IN  <: AbstractSPBasisState,
            SPBS_OUT <: AbstractSPBasisState,
            MPB_IN  <: MPBasis{N,SPBS_IN},
            MPB_OUT <: MPBasis{N,SPBS_OUT},
            MPO <: AbstractMPOperator{NO,MPB_IN},
            OP <: MPProjectorOperator{N,NO, SPBS_IN,SPBS_OUT, MPB_IN,MPB_OUT, MPO}
        }
    if haskey(io, :compact)
        print(io, "(MP projected) ")
        show(io, op.operator)
    else
        print(io, "Projection to new MP basis of operator:\n")
        show(io, op.operator)
        print(io, "New (outer basis) consists of "*string(length(basis(op)))*" states of type $(SPSSBS_OUT)\n")
    end
end


################################################################################
#   Interface functions
################################################################################

# obtain the current basis
function basis(operator :: MPProjectorOperator{N,NO, SPBS_IN,SPBS_OUT, MPB_IN,MPB_OUT, MPO}) :: MPB_OUT  where {
            N, NO,
            SPBS_IN  <: AbstractSPBasisState,
            SPBS_OUT <: AbstractSPBasisState,
            MPB_IN  <: MPBasis{N,SPBS_IN},
            MPB_OUT <: MPBasis{N,SPBS_OUT},
            MPO <: AbstractMPOperator{NO,MPB_IN}
        }
    return operator.basis
end

# obtain the matrix representation
function matrix_representation(operator :: MPProjectorOperator{N,NO, SPBS_IN,SPBS_OUT, MPB_IN,MPB_OUT, MPO}) :: Matrix{Complex{Float64}}  where {
            N, NO,
            SPBS_IN  <: AbstractSPBasisState,
            SPBS_OUT <: AbstractSPBasisState,
            MPB_IN  <: MPBasis{N,SPBS_IN},
            MPB_OUT <: MPBasis{N,SPBS_OUT},
            MPO <: AbstractMPOperator{NO,MPB_IN}
        }
    return operator.matrix_rep
end

# possibly recalculate the matrix representation
function recalculate!(operator :: MPProjectorOperator{N,NO, SPBS_IN,SPBS_OUT, MPB_IN,MPB_OUT, MPO}, recursive::Bool=true, basis_change::Bool=true) where {
            N, NO,
            SPBS_IN  <: AbstractSPBasisState,
            SPBS_OUT <: AbstractSPBasisState,
            MPB_IN  <: MPBasis{N,SPBS_IN},
            MPB_OUT <: MPBasis{N,SPBS_OUT},
            MPO <: AbstractMPOperator{NO,MPB_IN}
        }
    # maybe recalculate the inner operator
    if recursive
        recalculate!(operator.operator, recursive, basis_change)
    end
    # calculate new projector P_12
    if basis_change
        operator.projector_in_out = projector_matrix(operator.basis_in, operator.basis)
    end
    # calculate new matrix h_2 = P_12' * h_1 * P_12
    operator.matrix_rep = operator.projector_in_out' * matrix_representation(operator.operator) * operator.projector_in_out
    # return nothing
    return nothing
end

# set a parameter (returns (found parameter?, changed matrix?))
function set_parameter!(operator :: MPProjectorOperator{N,NO, SPBS_IN,SPBS_OUT, MPB_IN,MPB_OUT, MPO}, parameter :: Symbol, value; print_result::Bool=false, recalculate::Bool=true, kwargs...) where {
            N, NO,
            SPBS_IN  <: AbstractSPBasisState,
            SPBS_OUT <: AbstractSPBasisState,
            MPB_IN  <: MPBasis{N,SPBS_IN},
            MPB_OUT <: MPBasis{N,SPBS_OUT},
            MPO <: AbstractMPOperator{NO,MPB_IN}
        }
    # pass to inner operator
    (fp, cm) = set_parameter!(operator.operator, parameter, value; print_result=print_result, recalculate=recalculate, kwargs...)
    # maybe recalculate itself
    if recalculate && fp && cm
        recalculate!(operator, false, false)
    end
    # return the same variables
    return (fp, cm)
end

# get a parameter (returns (parameter value or nothing))
function get_parameter(operator :: MPProjectorOperator{N,NO, SPBS_IN,SPBS_OUT, MPB_IN,MPB_OUT, MPO}, parameter :: Symbol; print_result::Bool=false, kwargs...) where {
            N, NO,
            SPBS_IN  <: AbstractSPBasisState,
            SPBS_OUT <: AbstractSPBasisState,
            MPB_IN  <: MPBasis{N,SPBS_IN},
            MPB_OUT <: MPBasis{N,SPBS_OUT},
            MPO <: AbstractMPOperator{NO,MPB_IN}
        }
    # check if parameter can be set
    return get_parameter(operator.operator, parameter; print_result=print_result, kwargs...)
end

# get a list of parameters
function get_parameters(operator :: MPProjectorOperator{N,NO, SPBS_IN,SPBS_OUT, MPB_IN,MPB_OUT, MPO}; kwargs...) where {
            N, NO,
            SPBS_IN  <: AbstractSPBasisState,
            SPBS_OUT <: AbstractSPBasisState,
            MPB_IN  <: MPBasis{N,SPBS_IN},
            MPB_OUT <: MPBasis{N,SPBS_OUT},
            MPO <: AbstractMPOperator{NO,MPB_IN}
        }
    return get_parameters(operator.operator; kwargs...)
end




















# Single particle multi site projectors

################################################################################
#   Type definition
################################################################################

# Type definition of local MS operators (SPSS operator only acting on a single site of many)
mutable struct GeneralProjectorOperator{
        AB_IN  <: AbstractBasis,
        AB_OUT <: AbstractBasis,
        AO <: AbstractOperator{AB_IN}
    } <: AbstractOperator{AB_OUT}

    # the MS basis
    basis_in  :: AB_IN
    basis     :: AB_OUT
    # the multi particle operator
    operator  :: AO
    # the matrix representation of the projector
    projector_in_out :: Matrix{Complex{Float64}}
    # the matrix representation
    matrix_rep :: Matrix{Complex{Float64}}

    # custom constructor
    function GeneralProjectorOperator(operator :: AO, basis_to :: AB_OUT) where {
                AB_IN  <: AbstractBasis,
                AB_OUT <: AbstractBasis,
                AO <: AbstractOperator{AB_IN}
            }
        # create a new operator
        op = new{AB_IN, AB_OUT, AO}(basis(operator), basis_to, operator, projector_matrix(basis(operator), basis_to), zeros(length(basis_to), length(basis_to)))
        # recalculate the matrix representation
        recalculate!(op, true)
        # return the operator
        return op
    end
end

# export operator type
export  GeneralProjectorOperator



function ProjectorOperator(operator :: AO, basis_to :: AB_OUT) where {
            AB_IN  <: AbstractBasis,
            AB_OUT <: AbstractBasis,
            AO <: AbstractOperator{AB_IN}
        }

    return GeneralProjectorOperator(operator, basis_to)
end

export ProjectorOperator



# overwritten show function
import Base.show
function Base.show(io::IO, op::OP) where {
            AB_IN  <: AbstractBasis,
            AB_OUT <: AbstractBasis,
            AO <: AbstractOperator{AB_IN},
            OP <: GeneralProjectorOperator{AB_IN, AB_OUT, AO}
        }
    if haskey(io, :compact)
        print(io, "(A projected) ")
        show(io, op.operator)
    else
        print(io, "Projection to new abstract basis of operator:\n")
        show(io, op.operator)
        print(io, "New (outer basis) consists of "*string(length(basis(op)))*"\n")
    end
end


################################################################################
#   Interface functions
################################################################################

# obtain the current basis
function basis(operator :: GeneralProjectorOperator{AB_IN, AB_OUT, AO}) :: AB_OUT  where {
            AB_IN  <: AbstractBasis,
            AB_OUT <: AbstractBasis,
            AO <: AbstractOperator{AB_IN}
        }
    return operator.basis
end

# obtain the matrix representation
function matrix_representation(operator :: GeneralProjectorOperator{AB_IN, AB_OUT, AO}) :: Matrix{Complex{Float64}}  where {
            AB_IN  <: AbstractBasis,
            AB_OUT <: AbstractBasis,
            AO <: AbstractOperator{AB_IN}
        }
    return operator.matrix_rep
end

# possibly recalculate the matrix representation
function recalculate!(operator :: GeneralProjectorOperator{AB_IN, AB_OUT, AO}, recursive::Bool=true, basis_change::Bool=true) where {
            AB_IN  <: AbstractBasis,
            AB_OUT <: AbstractBasis,
            AO <: AbstractOperator{AB_IN}
        }
    # maybe recalculate the inner operator
    if recursive
        recalculate!(operator.operator, recursive, basis_change)
    end
    # calculate new projector P_12
    if basis_change
        operator.projector_in_out = projector_matrix(operator.basis_in, operator.basis)
    end
    # calculate new matrix h_2 = P_12' * h_1 * P_12
    operator.matrix_rep = operator.projector_in_out' * matrix_representation(operator.operator) * operator.projector_in_out
    # return nothing
    return nothing
end

# set a parameter (returns (found parameter?, changed matrix?))
function set_parameter!(operator :: GeneralProjectorOperator{AB_IN, AB_OUT, AO}, parameter :: Symbol, value; print_result::Bool=false, recalculate::Bool=true, kwargs...) where {
            AB_IN  <: AbstractBasis,
            AB_OUT <: AbstractBasis,
            AO <: AbstractOperator{AB_IN}
        }
    # pass to inner operator
    (fp, cm) = set_parameter!(operator.operator, parameter, value; print_result=print_result, recalculate=recalculate, kwargs...)
    # maybe recalculate itself
    if recalculate && fp && cm
        recalculate!(operator, false, false)
    end
    # return the same variables
    return (fp, cm)
end

# get a parameter (returns (parameter value or nothing))
function get_parameter(operator :: GeneralProjectorOperator{AB_IN, AB_OUT, AO}, parameter :: Symbol; print_result::Bool=false, kwargs...) where {
            AB_IN  <: AbstractBasis,
            AB_OUT <: AbstractBasis,
            AO <: AbstractOperator{AB_IN}
        }
    # check if parameter can be set
    return get_parameter(operator.operator, parameter; print_result=print_result, kwargs...)
end

# get a list of parameters
function get_parameters(operator :: GeneralProjectorOperator{AB_IN, AB_OUT, AO}; kwargs...) where {
            AB_IN  <: AbstractBasis,
            AB_OUT <: AbstractBasis,
            AO <: AbstractOperator{AB_IN}
        }
    return get_parameters(operator.operator; kwargs...)
end
