# General projector as MATRICES

# projection matrix P: |basis_to> = P |basis_from>
# |basis_to><basis_from|
"""
    projector_matrix(
        basis_to   :: B1,
        basis_from :: B2
    ) :: Matrix{Complex{Float64}} where {BS1<:AbstractSPBasisState, BS2<:AbstractSPBasisState, B1 <: AbstractBasis{BS1}, B2 <: AbstractBasis{BS2}}

This function computes the matrix projector P, defined as:

`` \\left| basis_{to} \\right> = P \\left| basis_{from} \\right> \\longrightarrow P=\\left|basis_{to} \\right> \\left< basis_{from} \\right| ``
"""
function projector_matrix(
            basis_to   :: B1,
            basis_from :: B2
        ) :: Matrix{Complex{Float64}} where {BS1<:AbstractSPBasisState, BS2<:AbstractSPBasisState, B1 <: AbstractBasis{BS1}, B2 <: AbstractBasis{BS2}}
    # create new matrix
    matrix = zeros(Complex{Float64}, length(basis_to),length(basis_from))
    # fill the matrix
    for i in 1:length(basis_to)
    for j in 1:length(basis_from)
        # fill in the element i,j (NOTE THE CONVENTION!!!)
        matrix[i,j] = overlap(basis_to[i], basis_from[j])
    end
    end
    # return the matrix
    return matrix
end
export projector_matrix


# projection matrix P: |basis_to> = P |basis_from>
# |basis_to><basis_from|
function projector_matrix(
            basis_to   :: B1,
            basis_from :: B2
        ) :: Matrix{Complex{Float64}} where {N,SPBS1<:AbstractSPBasisState,SPBS2<:AbstractSPBasisState, B1 <: MPBasis{N,SPBS1}, B2 <: MPBasis{N,SPBS2}}
    # create new matrix
    matrix = zeros(Complex{Float64}, length(basis_to),length(basis_from))
    # fill the matrix
    for i in 1:length(basis_to)
    for j in 1:length(basis_from)
        # fill in the element i,j (NOTE THE CONVENTION!!!)
        matrix[i,j] = overlap(basis_to, basis_to[i], basis_from, basis_from[j])
    end
    end
    # return the matrix
    return matrix
end
export projector_matrix
