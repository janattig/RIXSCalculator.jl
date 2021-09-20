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





#################
#  <B|c^d c|K>  #
#################

# <state1| c_a^dagger c_b |state2>
# creation annihilation
function expectation_value_ca(basis :: MPBasis{N,SPBS}, state_1 :: MPBasisState{N}, a::Int64, b::Int64, state_2 :: MPBasisState{N}) :: Complex{Float64} where {N,SPBS<:AbstractSPBasisState}
    # check all elements of state 2 explicitly
    if b in state_2.occupation
        # build a new state and return the overlap
        return overlap(basis, state_1, MPBasisState{N}(replace(state_2.occupation, b=>a)))
    else
        # otherwise, return 0
        return 0.0 + 0.0im
    end
end
export expectation_value_ca

# <state1| c_a^dagger c_b |state2>
# creation annihilation
# including buffer state
function expectation_value_ca!(basis :: MPBasis{N,SPBS}, state_1 :: MPBasisState{N}, a::Int64, b::Int64, state_2 :: MPBasisState{N}, state_buffer :: MPBasisState{N}) :: Complex{Float64} where {N,SPBS<:AbstractSPBasisState}
    # if no replace takes place, basis index is > 0
    # therefore sign will be returned in the overlap function
    # and product of signs is 0 which is then returned
    state_buffer.basis_index = -2
    state_buffer.basis_sign  = 0
    # set the elements of the buffer state respectively to state 2 elements
    @simd for i in 1:length(state_2.occupation)
        if state_2.occupation[i] == b
            state_buffer.occupation[i] = a
            state_buffer.basis_index   = -1
        else
            state_buffer.occupation[i] = state_2.occupation[i]
        end
    end
    # maybe 0 overlap if state 2 did not contain b
    if state_buffer.basis_index == -2
        return 0.0 + 0.0im
    end
    # return the overlap
    return overlap(basis, state_1, state_buffer)
end
export expectation_value_ca!











#######################
#  <B|c^d c c^d c|K>  #
#######################

# <state1| c_a^dagger c_b  c_c^dagger c_d |state2>
# creation annihilation
function expectation_value_caca(basis :: MPBasis{N,SPBS}, state_1 :: MPBasisState{N}, a::Int64, b::Int64, c::Int64, d::Int64, state_2 :: MPBasisState{N}) :: Complex{Float64} where {N,SPBS<:AbstractSPBasisState}
    # check all elements of state 2 explicitly
    if b == c && d in state_2.occupation
        # build a new state
        s2 = MPBasisState{N}(replace(state_2.occupation, d=>a))
        return overlap(basis, state_1, s2)
    elseif d in state_2.occupation && b in state_2.occupation
        # build a new state
        s2 = MPBasisState{N}(replace(state_2.occupation, d=>c, b=>a))
        return overlap(basis, state_1, s2)
    else
        # otherwise, return 0
        return 0.0 + 0.0im
    end
end
export expectation_value_caca
