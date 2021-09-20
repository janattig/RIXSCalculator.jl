function getSPCompositeState(
        basis_of_state :: SPBasis{SPSSBS},
        descriptor :: Tuple{Any, Any}
    ) where {SPSSBS <: AbstractSPSSBasisState}

    # find out the state
    state = getSPBasisState(SPSSBS, descriptor)
    # find out the index within the basis
    idx = 0
    for i in 1:length(basis_of_state)
        if state == basis_of_state[i]
            idx = i
        end
    end
    # check if index is actually within basis
    if idx < 1 || idx > length(basis_of_state)
        error("Could not find state $(descriptor) within basis $(basis_of_state)")
    end
    # construct a state vector
    state_vector = zeros(Complex{Float64}, length(basis_of_state))
    # increment entry of basis state
    state_vector[idx] = 1
    # return the vector
    return CompositeBasisState(state_vector, basis_of_state)
end
export getSPCompositeState