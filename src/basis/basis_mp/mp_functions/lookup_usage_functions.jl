# resetting the index lookup
function resetLookupIndex!(basis :: MPBasis{N,SPBS}) where {N,SPBS<:AbstractSPBasisState}
    # set a new dictonary
    basis.lookup_index = Dict{Vector{Int64}, Tuple{Int64,Int64}}()
    # use the index function on all basis states
    @simd for bs in basis
        # remove index
        bs.basis_index = -1
        bs.basis_sign  =  0
        # recalculate index
        index_and_sign!(basis, bs)
    end
end
# resetting the SP state lookup
function resetLookupSPStates!(basis :: MPBasis{N,SPBS}) where {N,SPBS<:AbstractSPBasisState}
    # set new lists i in which all states are belonging to site i
    basis.lookup_sp_states = Vector{Int64}[
        Int64[j for j in 1:length(basis) if i in basis[j].occupation]
        for i in 1:length(basis.single_particle_basis)
    ]
end


# calculate index and sign
function calculate_index_and_sign(basis :: MPBasis{N,SPBS}, state :: MPBasisState{N}) ::Tuple{Int64, Int64} where {N,SPBS<:AbstractSPBasisState}
    # obtain the sorted permutation of the occupation as well as the permutation sign
    (s, perm_occupation) = permutated(state.occupation)
    # search the element in the basis
    i = -1
    for j in 1:length(basis)
        @inbounds if basis[j].occupation == perm_occupation
            i = j
            break
        end
    end
    # if not present in basis, set sign to 0
    if i == -1
        s = 0
    end
    # return (sign, index)
    (s, i)
end


# using the index lookup (also for sites that dont belong to the basis)
function index_and_sign!(basis :: MPBasis{N,SPBS}, state :: MPBasisState{N}) where {N,SPBS<:AbstractSPBasisState}
    # check if the state has a set index
    if state.basis_index > 0
        return nothing
    end
    # otherwise check if in lookup table, otherwise add it
    # https://docs.julialang.org/en/v1/base/collections/#Base.get!-Tuple{Function,Any,Any}
    try
        # obtain sign and index and set in the object
        state.basis_sign, state.basis_index = basis.lookup_index[state.occupation]
        # return nothing since the values are set
        return nothing
    catch KeyError
        # print the occupation that is missing
        #println(state.occupation)
        # create the lookup reference
        basis.lookup_index[copy(state.occupation)] = calculate_index_and_sign(basis, state)
        # obtain sign and index and set in the object
        state.basis_sign, state.basis_index = basis.lookup_index[state.occupation]
        # return nothing since the values are set
        return nothing
    end
end
