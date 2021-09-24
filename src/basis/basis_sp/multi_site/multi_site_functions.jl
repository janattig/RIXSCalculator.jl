#####################################
#  CONVERSION between               #
#  single particle -single site and #
#  single particle-multi site       #
#####################################


# convert any single particle single site basis to multiple sites
"""
    getMultiSiteBasis(basis :: SPBasis{BS}, n :: Integer) :: SPBasis{SPMSBasisState{BS}} where {BS<:AbstractSPSSBasisState}

This function converts any single particle - single site basis to a single particle - multi site one. The integer input number `n` represents the number of sites of the new basis.
"""
function getMultiSiteBasis(basis :: SPBasis{BS}, n :: Integer) :: SPBasis{SPMSBasisState{BS}} where {BS<:AbstractSPSSBasisState}
    # make a list of multisite states
    multisite_states = SPMSBasisState{BS}[]
    # push all states for all sites
    for i in 1:n
    for b in states(basis)
        push!(multisite_states, SPMSBasisState{BS}(b, i))
    end
    end
    # return the multisite basis
    return SPBasis{SPMSBasisState{BS}}(multisite_states)
end
export getMultiSiteBasis


# convert a single particle multi site basis to a single site basis of site i
# by selection of basis elements
"""
    getSingleSiteBasis(basis :: SPBasis{SPMSBasisState{BS}}, i :: Integer) :: SPBasis{BS} where {BS<:AbstractSPSSBasisState}

This function converts a single particle - multi site basis to a single particle - single site basis of site `i`.
"""
function getSingleSiteBasis(basis :: SPBasis{SPMSBasisState{BS}}, i :: Integer) :: SPBasis{BS} where {BS<:AbstractSPSSBasisState}
    # make a list of all relevant single site states
    singlesite_states = BS[]
    # push all states for all sites
    for b in states(basis)
        if b.site == i
            push!(singlesite_states, b.state)
        end
    end
    # return the multisite basis
    return SPBasis{BS}(singlesite_states)
end
export getSingleSiteBasis
