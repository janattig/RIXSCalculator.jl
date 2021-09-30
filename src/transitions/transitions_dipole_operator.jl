# Type Definition of Dipole Operator
"""
    DipoleOperator <: AbstractSPMSOperator{SPBasis{SPMSBasisState{BasisStateXYZ}}}

`DipoleOperator` is a `mutable struct` defined by:

a single particle basis `basis :: SPBasis{SPMSBasisState{BasisStateXYZ}}` 

the current edge `edge :: Int64`

the core hole basis `basis_core :: SPBasis{SPMSBasisState{SPSSCompositeBasisState{SPBasis{BasisStateLS}}}}` (p orbitals or core with j=1/2 or j=3/2)

the spin quantization frame `spin_quantization:: CoordinateFrame` 

the core- hole projection `core_hole_projection :: Matrix{Complex{Float64}}`

the site `site :: Int64` and its position `position :: Vector{Float64}`

the ingoing polarization `eps_in :: Vector{Float64}`, the ingoing momentum `q_in :: Vector{Float64}`

the out going polarization `eps_out :: Vector{Float64}`, the ingoing momentum `q_out :: Vector{Float64}`

dipole matrices `D_x`, `D_y`, `D_z`, all of type `Matrix{Complex{Float64}}`

and the current matrix representation `matrix_rep :: Matrix{Complex{Float64}}` without prefactor.
"""
mutable struct DipoleOperator <: AbstractSPMSOperator{SPBasis{SPMSBasisState{BasisStateXYZ}}}

    # the inner basis (XYZ type basis)
    basis :: SPBasis{SPMSBasisState{BasisStateXYZ}}

    # the current edge
    edge :: Int64
    # the core hole basis (p orbitals or core with j=1/2 or j=3/2)
    basis_core :: SPBasis{SPMSBasisState{SPSSCompositeBasisState{SPBasis{BasisStateLS}}}}

    # spin quantization frame
    spin_quantization :: CoordinateFrame

    # the core hole projection (L=2 / L=3)
    core_hole_projection :: Matrix{Complex{Float64}}

    # the site index
    site :: Int64

    # the site position
    position :: Vector{Float64}

    # the ingoing polarization
    eps_in  :: Vector{Float64}
    # the ingoing momentum
    q_in    :: Vector{Float64}

    # the outgoing polarization
    eps_out :: Vector{Float64}
    # the outgoing momentum
    q_out   :: Vector{Float64}

    # dipole matrices
    D_x :: Matrix{Complex{Float64}}
    D_y :: Matrix{Complex{Float64}}
    D_z :: Matrix{Complex{Float64}}

    # the current matrix representation (without prefactor)
    matrix_rep :: Matrix{Complex{Float64}}

end

"""
    DipoleOperator(
                ;
                edge::Int64=3,
                site::Int64=1,
                site_position::Vector{<:Real}=[0,0,0],
                eps_in  ::Vector{<:Real}=[1,0,0],
                eps_out ::Vector{<:Real}=[1,0,0],
                q_in    ::Vector{<:Real}=[0,0,0],
                q_out   ::Vector{<:Real}=[0,0,0]
            ) :: DipoleOperator

Creates an object of type `DipoleOperator` with the input values.
"""
function DipoleOperator(
            ;
            edge::Int64=3,
            site::Int64=1,
            site_position::Vector{<:Real}=[0,0,0],
            eps_in  ::Vector{<:Real}=[1,0,0],
            eps_out ::Vector{<:Real}=[1,0,0],
            q_in    ::Vector{<:Real}=[0,0,0],
            q_out   ::Vector{<:Real}=[0,0,0]
        ) :: DipoleOperator

    # create a new object
    op = DipoleOperator(
        SPBasis{SPMSBasisState{BasisStateXYZ}}([ SPMSBasisState{BasisStateXYZ}(state, site) for state in states(getT2GBasisXYZ()) ]),
        edge,
        SPBasis{SPMSBasisState{SPSSCompositeBasisState{SPBasis{BasisStateLS}}}}([]),
        CoordinateFrame(),
        zeros(Complex{Float64},2,2),
        site,
        site_position,
        eps_in,
        q_in,
        eps_out,
        q_out,
        zeros(Complex{Float64},2,2),
        zeros(Complex{Float64},2,2),
        zeros(Complex{Float64},2,2),
        zeros(Complex{Float64},2,2),
    )

    # recalculate everything
    recalculate!(op, true, true)

    # return the object
    return op
end

function DipoleOperator(
        basis :: SPBasis{SPMSBasisState{SPSS}}
        ;
        kwargs...
    ) where {SPSS <: AbstractSPSSBasisState}

    # return the projection to that basis
    return ProjectorOperator(DipoleOperator(;kwargs...), basis)
end

function DipoleOperator(
        basis :: MPB
        ;
        kwargs...
    ) where {N, SPSS <: AbstractSPSSBasisState, SPBS <: SPMSBasisState{SPSS}, MPB <: MPBasis{N,SPBS}}

    # return the 1p generalization of projection to that basis
    return MPGeneralizedSPOperator(basis, ProjectorOperator(DipoleOperator(;kwargs...), basis.single_particle_basis))
end