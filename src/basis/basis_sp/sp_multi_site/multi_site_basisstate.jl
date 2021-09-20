#########################################################
# Definition of SINGLE PARTICLE MULTI SITE Basis states #
#########################################################


struct SPMSBasisState{B<:AbstractSPSSBasisState} <: AbstractSPBasisState
    # the basis state itself
    state :: B
    # the site index
    site  :: Int64
end
export SPMSBasisState


# custom print function
import Base.show
function Base.show(io::IO, state::SPMSBasisState{BS}) where {BS}
    show(io, state.state)
    print(io, " @ site(" * string(state.site) * ")")
end
# custom summary function
function summary(bs::SPMSBasisState{BS}, brackets="()") where {BS}
    spssbss = summary(bs.state, brackets)
    spssbss = replace(spssbss, brackets[1]=>brackets[1]*"#$(bs.site),")
    return spssbss
end

# custom print function
import Base.summary
function Base.summary(io::IO, basis::SPBasis{SPMSBasisState{BS}}) where {BS}
    print(io, string(length(basis))*"-element SP multi-site basis for states of type "*string(BS))
end
function get_sites(basis::SPBasis{BS}) where {SPSSBS, BS <: SPMSBasisState{SPSSBS}}
    return unique([s.site for s in states(basis)])
end