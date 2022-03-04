struct TetramerBasisState{B<:AbstractSPSSBasisState} <: AbstractSPBasisState
    # the SPSS basis state itself
    state :: B
    # bonding type
    bonding_type :: Symbol
    # coefficients
#     coeff :: Vector{Complex{Float64}}
    # sites
    site1 :: Int64
    site2 :: Int64
    site3 :: Int64
    site4 :: Int64
end


# custom print function

import Base.show
# function Base.show(io::IO, state::TetramerBasisState{BS} where BS) 
    
#     # note: state.state must be of a type whose name is "BasisState"*s, where s is a string of arbitrary legth
    
#     bsstr =  haskey(io, :compact) ? "" : string(typeof(state.state))[11:length(string(typeof(state.state)))]
    
#     s=bsstr*summary(state.state, ["|#"*string(state.site1)*",","⟩ "])
#     if state.bonding_type == :A0a || state.bonding_type == :A0b     
#         s=s*"+"
#     else
#         s=s*"-"
#     end
#     s=s*" "*bsstr*summary(state.state, ["|#"*string(state.site2)*",","⟩ "])
#     if state.bonding_type == :A0a || state.bonding_type == :A0c   
#         s=s*"+"
#     else
#         s=s*"-"
#     end
#     s=s*" "*bsstr*summary(state.state, ["|#"*string(state.site3)*",","⟩ "])
#     if state.bonding_type == :A0a || state.bonding_type == :A0d     
#         s=s*"+"
#     else
#         s=s*"-"
#     end
#     s=s*" "*bsstr*summary(state.state, ["|#"*string(state.site4)*",","⟩ "])
    
#     print(io, s)
# end

# version 2
function Base.show(io::IO, state::TetramerBasisState{BS} where BS) 
    
    bsstr =  haskey(io, :compact) ? "" : string(typeof(state.state))[11:length(string(typeof(state.state)))]
    s=bsstr*""
    
    if state.bonding_type == :Ap
        print(io,bsstr*"|"*string(state.state.orbital)*",B+,"*(state.state.ms>0 ? '↑' : '↓')*"⟩ ")
    elseif state.bonding_type == :Am
        print(io,bsstr*"|"*string(state.state.orbital)*",B-,"*(state.state.ms>0 ? '↑' : '↓')*"⟩ ")
    elseif state.bonding_type == :Ap
        print(io,bsstr*"|"*string(state.state.orbital)*",A+,"*(state.state.ms>0 ? '↑' : '↓')*"⟩ ")
    else # :Am
        print(io,bsstr*"|"*string(state.state.orbital)*",A-,"*(state.state.ms>0 ? '↑' : '↓')*"⟩ ")
    end
    
end



# custom summary function
function summary(bs::TetramerBasisState, brackets="()")
    return brackets[1]*"$(bs.bonding_type)"*summary(bs.state,[',', brackets[2]])
end


# export the basis type
export TetramerBasisState



#########################################
#   Pre-implemented tetramer basis   #
#########################################



"""
    getTetramerBasis(basis :: SPBasis{BS}, site1::Int64, site2::Int64, site3::Int64, site4::Int64) :: SPBasis{TetramerBasisState{BS}} where {BS<:AbstractSPSSBasisState}

This function provides the pre-implemented single particle - four site delocalized basis built on the SPSS basis states `BS`.
"""
function getTetramerBasis(basis :: SPBasis{BS}, site1::Int64, site2::Int64, site3::Int64, site4::Int64) :: SPBasis{TetramerBasisState{BS}} where {BS<:AbstractSPSSBasisState}
    # make a list of multisite states
    tetramer_states = TetramerBasisState{BS}[]
    # push all states for all sites
    for b in states(basis)
        push!(tetramer_states, TetramerBasisState{BS}(b, :Bp, site1, site2, site3, site4))
        push!(tetramer_states, TetramerBasisState{BS}(b, :Bm, site1, site2, site3, site4))
        push!(tetramer_states, TetramerBasisState{BS}(b, :Ap, site1, site2, site3, site4))
        push!(tetramer_states, TetramerBasisState{BS}(b, :Am, site1, site2, site3, site4))
    end
    # return the multisite basis
    return SPBasis{TetramerBasisState{BS}}(tetramer_states)
end

export getTetramerBasis