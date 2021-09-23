##############################################################
#
#   MP operators for general interaction Hamiltonians
#   - type definition Perkins Woelfle
#   - interface functions
#   - convenience functions
#
##############################################################


##############################################################
#   Type definition
##############################################################

# define a Perkins Woelfe Hamiltonian for 2 particle interation of ELECTRONS
"""
    mutable struct MPElectronPerkinsWoelfleHamiltonian{MPB <: MPBasis{N,SPBS} where {N,SPBS}} <: AbstractMPInteractionHamiltonian{2,MPB}

This object defines the interaction Hamiltonian described in the Perskin-Sizyuk-Woelfle "paper" [ https://arxiv.org/abs/1311.0852v2 ].

The hamiltonian is defined as follows in the electron picture:

`` \\sum_{i,\\alpha} n_{i\\alpha\\uparrow}n_{i\\alpha\\downarrow} +  
   \\sum_{i, \\sigma,\\alpha \\neq \\alpha^\\prime} n_{i\\alpha\\sigma}n_{i\\alpha^\\prime\\sigma} +
   \\sum_{i,\\alpha\\neq\\alpha^\\prime} n_{i\\alpha\\uparrow}n_{i\\alpha^\\prime\\downarrow} +
   \\sum_{i,\\alpha\\neq\\alpha^\\prime} d^\\dagger_{i\\alpha\\uparrow} d^\\dagger_{i\\alpha\\downarrow}d_{i\\alpha^\\prime\\downarrow}d_{i\\alpha^\\prime\\uparrow} + 
   \\sum_{i,\\alpha\\neq\\alpha^\\prime} d^\\dagger_{i\\alpha\\uparrow} d_{i\\alpha\\downarrow}d^\\dagger_{i\\alpha^\\prime\\downarrow}d_{i\\alpha^\\prime\\uparrow}``

"""
mutable struct MPElectronPerkinsWoelfleHamiltonian{
        MPB <: MPBasis{N,SPBS} where {N,SPBS}
    } <: AbstractMPInteractionHamiltonian{2,MPB}

    # the MP basis
    basis :: MPB
    # the matrix representation
    matrix_rep :: Matrix{Complex{Float64}}

    # the site
    site :: Int64

    # the parameters
    u1 :: Float64
    u2 :: Float64
    jH :: Float64

    # the sub hamiltonians
    # density-density interaction
    op_den_den_same_orb  :: MPElectronDensityDensityOperator{MPB}
    op_den_den_same_spin :: MPElectronDensityDensityOperator{MPB}
    op_den_den_remainder :: MPElectronDensityDensityOperator{MPB}
    # two-particle scattering
    op_2p_sc_spin_cons   :: MPElectron2PScatteringOperator{MPB}
    op_2p_sc_spin_flip   :: MPElectron2PScatteringOperator{MPB}

    # Simplified constructor
    function MPElectronPerkinsWoelfleHamiltonian(basis::MPB, site::Int64, u1::Real, u2::Real, jH::Real) where {
                N, SPBS <: AbstractSPBasisState,
                MPB <: MPBasis{N,SPBS}
            }
        # define new operators
        op_den_den_same_orb  = generateDensityDensityElectronInteractionSameOrbital(basis, site, u1)
        op_den_den_same_spin = generateDensityDensityElectronInteractionSameSpin(basis, site, 0.5*(u2-jH))
        op_den_den_remainder = generateDensityDensityElectronInteractionRemainder(basis, site, u2)
        op_2p_sc_spin_cons   = generate2PScatteringElectronInteractionSpinConserve(basis, site, jH)
        op_2p_sc_spin_flip   = generate2PScatteringElectronInteractionSpinFlip(basis, site, -jH)
        # construct new operator with those suboperators
        op = new{MPB}(basis, zeros(Complex{Float64}, length(basis),length(basis)), site, u1,u2,jH, op_den_den_same_orb, op_den_den_same_spin, op_den_den_remainder, op_2p_sc_spin_cons, op_2p_sc_spin_flip)
        # let it recalculate
        recalculate!(op, false)
        # return it
        return op
    end
end
export MPElectronPerkinsWoelfleHamiltonian

# define a Perkins Woelfe Hamiltonian for 2 particle interation of ELECTRONS
mutable struct MPHolePerkinsWoelfleHamiltonian{
        MPB <: MPBasis{N,SPBS} where {N,SPBS}
    } <: AbstractMPInteractionHamiltonian{2,MPB}

    # the MP basis
    basis :: MPB
    # the matrix representation
    matrix_rep :: Matrix{Complex{Float64}}

    # the site
    site :: Int64

    # the parameters
    u1 :: Float64
    u2 :: Float64
    jH :: Float64

    # the sub hamiltonians
    # density-density interaction
    op_den_den_same_orb  :: MPHoleDensityDensityOperator{MPB}
    op_den_den_same_spin :: MPHoleDensityDensityOperator{MPB}
    op_den_den_remainder :: MPHoleDensityDensityOperator{MPB}
    # two-particle scattering
    op_2p_sc_spin_cons   :: MPHole2PScatteringOperator{MPB}
    op_2p_sc_spin_flip   :: MPHole2PScatteringOperator{MPB}

    # Simplified constructor
    function MPHolePerkinsWoelfleHamiltonian(basis::MPB, site::Int64, u1::Real, u2::Real, jH::Real) where {
                N, SPBS <: AbstractSPBasisState,
                MPB <: MPBasis{N,SPBS}
            }
        # define new operators
        op_den_den_same_orb  = generateDensityDensityHoleInteractionSameOrbital(basis, site, u1)
        op_den_den_same_spin = generateDensityDensityHoleInteractionSameSpin(basis, site, 0.5*(u2-jH))
        op_den_den_remainder = generateDensityDensityHoleInteractionRemainder(basis, site, u2)
        op_2p_sc_spin_cons   = generate2PScatteringHoleInteractionSpinConserve(basis, site, jH)
        op_2p_sc_spin_flip   = generate2PScatteringHoleInteractionSpinFlip(basis, site, -jH)
        # construct new operator with those suboperators
        op = new{MPB}(basis, zeros(Complex{Float64}, length(basis),length(basis)), site, u1,u2,jH, op_den_den_same_orb, op_den_den_same_spin, op_den_den_remainder, op_2p_sc_spin_cons, op_2p_sc_spin_flip)
        # let it recalculate
        recalculate!(op, false)
        # return it
        return op
    end
end
export MPHolePerkinsWoelfleHamiltonian

import Base.show
function Base.show(io::IO, op::MPHolePerkinsWoelfleHamiltonian{MPB}) where {
            N, SPBS <: AbstractSPBasisState,
            MPB <: MPBasis{N,SPBS}
        }
    if haskey(io, :compact)
        if (op.u2 == op.u1-2*op.jH)
            print(io, "2-particle @site("*string(op.site)*") (hole) Perkins-Woelfle hamiltonian with U="*string(op.u1)*", J_H="*string(op.jH)*"")
        else
            print(io, "2-particle @site("*string(op.site)*") (hole) Perkins-Woelfle hamiltonian with U1="*string(op.u1)*", U2="*string(op.u2)*", J_H="*string(op.jH)*"")
        end
    else
        print(io, "2-particle (hole) Perkins-Woelfle hamiltonian @site("*string(op.site)*")\n")
        if (op.u2 == op.u1-2*op.jH)
            print(io, "--> parameter   U = "*string(op.u1)*"\n")
            print(io, "--> parameter J_H = "*string(op.jH)*"\n")
        else
            print(io, "--> parameter  U1 = "*string(op.u1)*"\n")
            print(io, "--> parameter  U2 = "*string(op.u2)*"\n")
            print(io, "--> parameter J_H = "*string(op.jH)*"\n")
        end
        print(io, "Multi-particle basis contains "*string(length(basis(op)))*" states in total, with "*string(N)*" particles per state\n")
    end
end


import Base.show
function Base.show(io::IO, op::MPElectronPerkinsWoelfleHamiltonian{MPB}) where {
            N, SPBS <: AbstractSPBasisState,
            MPB <: MPBasis{N,SPBS}
        }
    if haskey(io, :compact)
        if (op.u2 == op.u1-2*op.jH)
            print(io, "2-particle @site("*string(op.site)*") (electron) Perkins-Woelfle hamiltonian with U="*string(op.u1)*", J_H="*string(op.jH)*"")
        else
            print(io, "2-particle @site("*string(op.site)*") (electron) Perkins-Woelfle hamiltonian with U1="*string(op.u1)*", U2="*string(op.u2)*", J_H="*string(op.jH)*"")
        end
    else
        print(io, "2-particle (electron) Perkins-Woelfle hamiltonian @site("*string(op.site)*")\n")
        if (op.u2 == op.u1-2*op.jH)
            print(io, "--> parameter   U ="*string(op.u1)*"\n")
            print(io, "--> parameter J_H ="*string(op.jH)*"\n")
        else
            print(io, "--> parameter  U1 ="*string(op.u1)*"\n")
            print(io, "--> parameter  U2 ="*string(op.u2)*"\n")
            print(io, "--> parameter J_H ="*string(op.jH)*"\n")
        end
        print(io, "Multi-particle basis contains "*string(length(basis(op)))*" states in total, with "*string(N)*" particles per state\n")
    end
end




##############################################################
#   Interface functions
##############################################################

# obtain the current basis (ELECTRON & HOLE)
function basis(operator :: MPIH) :: MPB where {
            N, SPBS <: AbstractSPBasisState,
            MPB <: MPBasis{N,SPBS},
            MPIH <: AbstractMPInteractionHamiltonian{2,MPB}
        }
    return operator.basis
end

# obtain the matrix representation (ELECTRON & HOLE)
function matrix_representation(operator :: MPIH) :: Matrix{Complex{Float64}} where {
            N, SPBS <: AbstractSPBasisState,
            MPB <: MPBasis{N,SPBS},
            MPIH <: AbstractMPInteractionHamiltonian{2,MPB}
        }
    return operator.matrix_rep
end


# possibly recalculate the matrix representation (ELECTRON & HOLE) (Fallback for non XYZ)
function recalculate!(operator :: MPIH, recursive::Bool=true) where {
            N, SPBS <: AbstractSPBasisState,
            MPB <: MPBasis{N,SPBS},
            MPIH <: AbstractMPInteractionHamiltonian{2,MPB}
        }
    @error "currently only recalculation of density-density operator implemented for XYZ basis states, not basis states of type $(SPBS)" stacktrace()
end

# possibly recalculate the matrix representation (ELECTRON & HOLE)
function recalculate!(operator :: MPIH, recursive::Bool=true, basis_change::Bool=true) where {
            N, SPBS <: Union{SPMSBasisState{BasisStateXYZ}, BasisStateXYZ},
            MPB <: MPBasis{N,SPBS},
            MPIH <: AbstractMPInteractionHamiltonian{2,MPB}
        }
    # set all prefactors
    operator.op_den_den_same_orb.prefactor  = operator.u1
    operator.op_den_den_same_spin.prefactor = 0.5*(operator.u2-operator.jH)
    operator.op_den_den_remainder.prefactor = operator.u2
    operator.op_2p_sc_spin_cons.prefactor   = operator.jH
    operator.op_2p_sc_spin_flip.prefactor   = -operator.jH
    # maybe recalculate recursively
    if recursive
        recalculate!(operator.op_den_den_same_orb)
        recalculate!(operator.op_den_den_same_spin)
        recalculate!(operator.op_den_den_remainder)
        recalculate!(operator.op_2p_sc_spin_cons)
        recalculate!(operator.op_2p_sc_spin_flip)
    end
    # create new matrix by summing all contributions of suboperators
    operator.matrix_rep =
        matrix_representation(operator.op_den_den_same_orb)  .+
        matrix_representation(operator.op_den_den_same_spin) .+
        matrix_representation(operator.op_den_den_remainder) .+
        matrix_representation(operator.op_2p_sc_spin_cons)   .+
        matrix_representation(operator.op_2p_sc_spin_flip)
    # return nothing
    return nothing
end




# setting interaction parameter (U, J_H)
function set_parameter_ujH(operator :: MPIH, u::Real, jH::Real) where {
            N, SPBS <: AbstractSPBasisState,
            MPB <: MPBasis{N,SPBS},
            MPIH <: AbstractMPInteractionHamiltonian{2,MPB}
        }
    # set the parameters
    operator.jH = jH
    operator.u1 = u
    operator.u2 = u-2jH
    # recalculate, but not recursive
    recalculate!(operator, true)
end
# setting interaction parameter (U1, U2, J_H)
function set_parameter_u1u2jH(operator :: MPIH, u1::Real, u2::Real, jH::Real) where {
            N, SPBS <: AbstractSPBasisState,
            MPB <: MPBasis{N,SPBS},
            MPIH <: AbstractMPInteractionHamiltonian{2,MPB}
        }
    # set the parameters
    operator.jH = jH
    operator.u1 = u1
    operator.u2 = u2
    # recalculate, but not recursive
    recalculate!(operator, true)
end



# set parameter interface (returns (found parameter?, changed matrix?))
function set_parameter!(operator :: MPIH, parameter :: Symbol, value; print_result::Bool=true, recalculate::Bool=true, site::Union{Int64, Symbol}=-1, kwargs...) where {
            N, SPBS <: Union{SPMSBasisState{BasisStateXYZ}, BasisStateXYZ},
            MPB <: MPBasis{N,SPBS},
            MPIH <: AbstractMPInteractionHamiltonian{2,MPB}
        }
    # pass on to containing operator
    if site == operator.site || site == :all
        # check which parameter is adressed
        if parameter == :U
            if recalculate && (operator.u1 != value || operator.u2 != operator.u1-2*operator.jH)
                operator.u1 = value
                operator.u2 = operator.u1-2*operator.jH
                recalculate!(operator, true, false)
            else
                operator.u1 = value
                operator.u2 = operator.u1-2*operator.jH
            end
            if print_result
                println("Parameter :$(parameter) found and (u1,u2) set to values ($(operator.u1),$(operator.u2))")
            end
            return (true, true)
        elseif parameter == :U1
            if recalculate && operator.u1 != value
                operator.u1 = value
                recalculate!(operator, true, false)
            else
                operator.u1 = value
            end
            if print_result
                println("Parameter :$(parameter) found and u1 set to value $(operator.u1)")
            end
            return (true, true)
        elseif parameter == :U2
            if recalculate && operator.u2 != value
                operator.u2 = value
                recalculate!(operator, true, false)
            else
                operator.u2 = value
            end
            if print_result
                println("Parameter :$(parameter) found and u2 set to value $(operator.u2)")
            end
            return (true, true)
        elseif parameter == :J_H
            if operator.u2 != operator.u1-2*operator.jH
                if recalculate && operator.jH != value
                    operator.jH = value
                    recalculate!(operator, true, false)
                else
                    operator.jH = value
                end
                if print_result
                    println("Parameter :$(parameter) found and jH set to value $(operator.jH)")
                end
            else
                if recalculate && (operator.jH != value || operator.u2 != operator.u1-2*operator.jH)
                    operator.jH = value
                    operator.u2 = operator.u1-2*operator.jH
                    recalculate!(operator, true, false)
                else
                    operator.jH = value
                    operator.u2 = operator.u1-2*operator.jH
                end
                if print_result
                    println("Parameter :$(parameter) found and (u2,jH) set to values ($(operator.u2),$(operator.jH))")
                end
            end
            return (true, true)
        else
            if print_result
                println("Parameter :$(parameter) not found")
            end
            return (false, false)
        end
    else
        if print_result
            println("site $(site) is not matching to site of operator: site=$(operator.site)")
        end
        return (false, false)
    end
end


# set parameter interface (returns (found parameter?, changed matrix?))
function get_parameter(operator :: MPIH, parameter :: Symbol; print_result::Bool=true, site::Union{Int64, Symbol}=-1, kwargs...) where {
            N, SPBS <: Union{SPMSBasisState{BasisStateXYZ}, BasisStateXYZ},
            MPB <: MPBasis{N,SPBS},
            MPIH <: AbstractMPInteractionHamiltonian{2,MPB}
        }
    # pass on to containing operator
    if site == operator.site || site == :all
        # check which parameter is adressed
        if parameter == :U
            if (operator.u2 == operator.u1-2*operator.jH)
                if print_result
                    println("Parameter :$(parameter) found and (u1,u2) were set in accordance, returning parameter")
                end
                return operator.u1
            else
                if print_result
                    println("Parameter :$(parameter) found BUT (u1,u2) were NOT set in accordance, returning nothing")
                end
                return nothing
            end
        elseif parameter == :U1
            if print_result
                println("Parameter :$(parameter) found and returned")
            end
            return operator.u1
        elseif parameter == :U2
            if print_result
                println("Parameter :$(parameter) found and returned")
            end
            return operator.u2
        elseif parameter == :J_H
            if print_result
                println("Parameter :$(parameter) found and returned")
            end
            return operator.jH
        else
            if print_result
                println("Parameter :$(parameter) not found")
            end
            return nothing
        end
    else
        if print_result
            println("site $(site) is not matching to site of operator: site=$(operator.site)")
        end
        return nothing
    end
end

# set parameter interface (returns (found parameter?, changed matrix?))
function get_parameters(operator :: MPIH; site::Union{Int64, Symbol}=:all, kwargs...) where {
            N, SPBS <: Union{SPMSBasisState{BasisStateXYZ}, BasisStateXYZ},
            MPB <: MPBasis{N,SPBS},
            MPIH <: AbstractMPInteractionHamiltonian{2,MPB}
        }
    # pass on to containing operator
    if site == operator.site || site == :all
        # check which parameter is adressed
        if (operator.u2 == operator.u1-2*operator.jH)
            return Symbol[:U, :J_H]
        else
            return Symbol[:U1, :U2, :J_H]
        end
    else
        Symbol[]
    end
end
