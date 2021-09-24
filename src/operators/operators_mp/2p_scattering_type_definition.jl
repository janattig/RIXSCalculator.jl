# abstract supertype for two-particle scattering operators
abstract type AbstractMP2PScatteringOperator{MPB} <: AbstractMP2POperator{MPB} end

# define a MP operator for 2 particle scattering interactions of ELECTRONS
mutable struct MPElectron2PScatteringOperator{
        MPB <: MPBasis{N,SPBS} where {N,SPBS}
    } <: AbstractMP2PScatteringOperator{MPB}

    # the MP basis
    basis :: MPB
    # the matrix representation
    matrix_rep :: Matrix{Complex{Float64}}
    # the list of interacting orbitals ((a,b),(c,d))
    # meaning: b->a, d->c scattering
    interacting_orbitals :: Vector{Tuple{Tuple{Int64,Int64},Tuple{Int64,Int64}}}
    # the prefactor
    prefactor :: Float64

    # custom constructor
    function MPElectron2PScatteringOperator(basis :: MPB, interacting_orbitals :: Vector{Tuple{Tuple{Int64,Int64},Tuple{Int64,Int64}}}, prefactor :: Real) where {
                N, SPBS <: AbstractSPBasisState,
                MPB <: MPBasis{N,SPBS}
            }
        # create a new operator
        op = new{MPB}(basis, zeros(length(basis), length(basis)), interacting_orbitals, prefactor)
        # recalculate the matrix representation
        recalculate!(op)
        # return the operator
        return op
    end
end
export MPElectron2PScatteringOperator

# define a MP operator for 2 particle scattering interactions of HOLES
mutable struct MPHole2PScatteringOperator{
        MPB <: MPBasis{N,SPBS} where {N,SPBS}
    } <: AbstractMP2PScatteringOperator{MPB}

    # the MP basis
    basis :: MPB
    # the matrix representation
    matrix_rep :: Matrix{Complex{Float64}}
    # the list of interacting orbitals ((a,b),(c,d))
    # meaning: b->a, d->c scattering
    interacting_orbitals :: Vector{Tuple{Tuple{Int64,Int64},Tuple{Int64,Int64}}}
    # the prefactor
    prefactor :: Float64

    # custom constructor
    function MPHole2PScatteringOperator(basis :: MPB, interacting_orbitals :: Vector{Tuple{Tuple{Int64,Int64},Tuple{Int64,Int64}}}, prefactor :: Real) where {
                N, SPBS <: AbstractSPBasisState,
                MPB <: MPBasis{N,SPBS}
            }
        # create a new operator
        op = new{MPB}(basis, zeros(length(basis), length(basis)), interacting_orbitals, prefactor)
        # recalculate the matrix representation
        recalculate!(op)
        # return the operator
        return op
    end
end
export MPHole2PScatteringOperator


import Base.show
function Base.show(io::IO, op::MPHole2PScatteringOperator{MPB}) where {
            N,
            SPBS <: AbstractSPBasisState,
            MPB <: MPBasis{N,SPBS}
        }
    if haskey(io, :compact)
        print(io, "2-particle (hole-hole) scattering with "*string(length(op.interacting_orbitals))*" (c*a)*(c*a) terms")
    else
        print(io, "2-particle (hole-hole) scattering\n")
        print(io, "Interactions include "*string(length(op.interacting_orbitals))*" (c*a)*(c*a) terms:\n")
        for orbital_pair in op.interacting_orbitals
            orb_1_1 = summary(basis(op).single_particle_basis[orbital_pair[1][1]], "{}")
            orb_1_2 = summary(basis(op).single_particle_basis[orbital_pair[1][2]], "{}")
            orb_2_1 = summary(basis(op).single_particle_basis[orbital_pair[2][1]], "{}")
            orb_2_2 = summary(basis(op).single_particle_basis[orbital_pair[2][2]], "{}")
            print(io, " + (c_"*orb_1_1*"*a_"*orb_1_2*")*(c_"*orb_2_1*"*a_"*orb_2_2*")\n")
        end
        print(io, "Overall prefactor is "*string(op.prefactor)*"\n")
        print(io, "Multi-particle basis contains "*string(length(basis(op)))*" states in total, with "*string(N)*" particles per state\n")
    end
end
function Base.show(io::IO, op::MPElectron2PScatteringOperator{MPB}) where {
            N,
            SPBS <: AbstractSPBasisState,
            MPB <: MPBasis{N,SPBS}
        }
    if haskey(io, :compact)
        print(io, "2-particle (electron-electron) scattering with "*string(length(op.interacting_orbitals))*" (c*a)*(c*a) terms")
    else
        print(io, "2-particle (electron-electron) scattering\n")
        print(io, "Interactions include "*string(length(op.interacting_orbitals))*" (c*a)*(c*a) terms:\n")
        for orbital_pair in op.interacting_orbitals
            orb_1_1 = summary(basis(op).single_particle_basis[orbital_pair[1][1]], "{}")
            orb_1_2 = summary(basis(op).single_particle_basis[orbital_pair[1][2]], "{}")
            orb_2_1 = summary(basis(op).single_particle_basis[orbital_pair[2][1]], "{}")
            orb_2_2 = summary(basis(op).single_particle_basis[orbital_pair[2][2]], "{}")
            print(io, " + (c_"*orb_1_1*"*a_"*orb_1_2*")*(c_"*orb_2_1*"*a_"*orb_2_2*")\n")
        end
        print(io, "Overall prefactor is "*string(op.prefactor)*"\n")
        print(io, "Multi-particle basis contains "*string(length(basis(op)))*" states in total, with "*string(N)*" particles per state\n")
    end
end



##############################################################
#   Interface functions
##############################################################

# obtain the current basis (ELECTRON & HOLE)
function basis(operator :: MPDDOP) :: MPB where {
            N, SPBS <: AbstractSPBasisState,
            MPB <: MPBasis{N,SPBS},
            MPDDOP <: AbstractMP2PScatteringOperator{MPB}
        }
    return operator.basis
end

# obtain the matrix representation (ELECTRON & HOLE)
function matrix_representation(operator :: MPDDOP) :: Matrix{Complex{Float64}} where {
            N, SPBS <: AbstractSPBasisState,
            MPB <: MPBasis{N,SPBS},
            MPDDOP <: AbstractMP2PScatteringOperator{MPB}
        }
    return operator.matrix_rep .* operator.prefactor
end


# possibly recalculate the matrix representation (ELECTRON & HOLE) (Fallback for non XYZ)
function recalculate!(operator :: MPDDOP, recursive::Bool=true, basis_change::Bool=true) where {
            N, SPBS <: AbstractSPBasisState,
            MPB <: MPBasis{N,SPBS},
            MPDDOP <: AbstractMP2PScatteringOperator{MPB}
        }
    @error "currently only recalculation of density-density operator implemented for XYZ basis states, not basis states of type $(SPBS)" stacktrace()
end

# possibly recalculate the matrix representation (ELECTRON & HOLE)
function recalculate!(operator :: MPDDOP, recursive::Bool=true, basis_change::Bool=true) where {
            N, SPBS <: Union{SPMSBasisState{BasisStateXYZ}, BasisStateXYZ},
            MPB <: MPBasis{N,SPBS},
            MPDDOP <: AbstractMP2PScatteringOperator{MPB}
        }
    # create new matrix
    operator.matrix_rep = zeros(Complex{Float64}, length(basis(operator)), length(basis(operator)))
    # iterate over all interacting orbitals
    for orb in operator.interacting_orbitals
        # let it fill its contribution
        fill_orbital_contribution!(operator, orb[1], orb[2])
    end
end


# contribution of one orbital (ELEKTRON)
function fill_orbital_contribution!(operator :: MPElectron2PScatteringOperator{MPB}, scattering_1::Tuple{Int64,Int64}, scattering_2::Tuple{Int64,Int64}) where {
            N, SPBS <: AbstractSPBasisState,
            MPB <: MPBasis{N,SPBS}
        }
    # iterate over all suitable orbitals
    for alpha in basis(operator).lookup_sp_states[scattering_1[1]]
    for beta  in basis(operator).lookup_sp_states[scattering_2[2]]
        # fill matrix entry naively
        operator.matrix_rep[alpha, beta] += expectation_value_caca(
            basis(operator),
            basis(operator)[alpha],
            scattering_1[1], scattering_1[2], scattering_2[1], scattering_2[2],
            basis(operator)[beta]
        )
    end
    end
end

# contribution of one orbital (HOLE)
function fill_orbital_contribution!(operator :: MPHole2PScatteringOperator{MPB}, scattering_1::Tuple{Int64,Int64}, scattering_2::Tuple{Int64,Int64}) where {
            N, SPBS <: AbstractSPBasisState,
            MPB <: MPBasis{N,SPBS}
        }
    # iterate over all suitable orbitals
    for alpha in basis(operator).lookup_sp_states[scattering_2[2]]
    for beta  in basis(operator).lookup_sp_states[scattering_1[1]]
        # fill matrix entry naively
        operator.matrix_rep[alpha, beta] += expectation_value_caca(
            basis(operator),
            basis(operator)[alpha],
            scattering_1[2], scattering_1[1], scattering_2[2], scattering_2[1],
            basis(operator)[beta]
        )
    end
    end
end


export fill_orbital_contribution!
