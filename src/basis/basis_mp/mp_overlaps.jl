###########
#  <B|K>  #
###########

# overlap for multiparticle states including basis (fallback)
function overlap(basis :: MPBasis{N,SPBS}, state_1 :: MPBasisState{N1}, state_2 :: MPBasisState{N2}) :: Complex{Float64} where {N,N1,N2,SPBS<:AbstractSPBasisState}
    @error "Dimensions of multi particle overlap do not agree: $(N); $(N1), $(N2)" stacktrace()
    return NaN+ NaN*im
end


# overlap for multiparticle states including basis
function overlap(basis :: MPBasis{N,SPBS}, state_1 :: MPBasisState{N}, state_2 :: MPBasisState{N}) :: Complex{Float64} where {N,SPBS<:AbstractSPBasisState}
    # compute all permutations of state 2 and build delta functions according to
    index_and_sign!(basis, state_1)
    index_and_sign!(basis, state_2)
    # if both are included in the basis, simply return delta function
    if state_1.basis_index == state_2.basis_index && state_1.basis_index > 0
        return complex(float(state_1.basis_sign*state_2.basis_sign))
    else
        return 0.0 + 0.0im
    end
end



# overlap for multiparticle states including basis
function overlap(basis_1 :: MPBasis{N,SPBS1}, state_1 :: MPBasisState{N}, basis_2 :: MPBasis{N,SPBS2},  state_2 :: MPBasisState{N}) :: Complex{Float64} where {N,SPBS1<:AbstractSPBasisState,SPBS2<:AbstractSPBasisState}
    # compute all permutations of state 2 and build delta functions according to
    #sign1, index1 = index_and_sign(basis_1, state_1)
    #sign2, index2 = index_and_sign(basis_2, state_2)
    # get the single particle basis
    basis_sp_1 = basis_1.single_particle_basis
    basis_sp_2 = basis_2.single_particle_basis
    # init overlap
    overlap_mp = 0.0 + 0.0im
    # get all permutations of state 2
    for permutated_occ in permutations(state_2.occupation)
        # build the overlap according to permutation
        ovl = prod(Complex{Float64}[
            overlap(basis_sp_1[state_1.occupation[i]], basis_sp_2[permutated_occ[i]]) for i in 1:N
        ])
        psgn = permutation_sign!(permutated_occ)
        # add with sign to the complete overlap
        overlap_mp += ovl * psgn
    end
    return overlap_mp
end
