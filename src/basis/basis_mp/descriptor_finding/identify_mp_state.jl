function getMPState(
            basis :: MPBasis{N, SPSS},
            descriptors :: Union{
                Tuple{Integer, String, Real},
                Tuple{Integer, Symbol, Real},
                Tuple{Integer, Real, Real},
                Tuple{Integer, String, Symbol},
                Tuple{Integer, Symbol, Symbol},
                Tuple{Integer, Real, Symbol},
                Tuple{Integer, String, String},
                Tuple{Integer, Symbol, String},
                Tuple{Integer, Real, String}
            } ...
        ) where {N, SPSS<:AbstractSPBasisState}

    # get the basis state
    mp_basis_state = getMPBasisState(basis, descriptors...)
    # return new vector where there is a 1 at the fitting position
    mp_state = zeros(Complex{Float64}, length(basis))
    mp_state[mp_basis_state.basis_index] = 1.0 * mp_basis_state.basis_sign
    return mp_state
end
function getMPBasisState(
            basis :: MPBasis{N, SPSS},
            descriptors :: Union{
                Tuple{Integer, String, Real},
                Tuple{Integer, Symbol, Real},
                Tuple{Integer, Real, Real},
                Tuple{Integer, String, Symbol},
                Tuple{Integer, Symbol, Symbol},
                Tuple{Integer, Real, Symbol},
                Tuple{Integer, String, String},
                Tuple{Integer, Symbol, String},
                Tuple{Integer, Real, String}
            } ...
        ) where {N, SPSS<:AbstractSPBasisState}

    # check that the correct number of spss states defined
    @assert N==length(descriptors) "Wrong number of SPSS states described, $(N) needed, $(length(descriptors)) given"
    # get all basis states
    spss_states = SPSS[
        getSPBasisState(SPSS, d)
        for d in descriptors
    ]
    # find the SPSS indices in the basis
    spss_indices = Int64[
        index(basis.single_particle_basis, spss) for spss in spss_states
    ]
    # create new MP State
    mp_state_new = MPBasisState{N}(spss_indices)
    # set the index and sign
    index_and_sign!(basis, mp_state_new)
    # return the state
    return mp_state_new
end

#export
export getMPState,getMPBasisState