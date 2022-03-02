struct TdSymBasisState <: AbstractSPBasisState
    # symmetry type
    symmetry_type :: Symbol
    # spin
    ms :: Rational{Int64}
    # hopping parameters
    A :: Float64
    B :: Float64
    C :: Float64
    #beware though: when studying wrt hopping, have to build the basis every time...
    # sites
    site1 :: Int64
    site2 :: Int64
    site3 :: Int64
    site4 :: Int64
end

# custom show function
import Base.show
function Base.show(io::IO, state::TdSymBasisState) 
    
    bsstr =  haskey(io, :compact) ? "" : "Td"
        
    print(io, bsstr*"|"*string(state.symmetry_type)*","*(state.ms>0 ? '↑' : '↓')*"⟩")
end

# custom summary function
function summary(bs::TdSymBasisState, brackets="()")
    return brackets[1]*"$(bs.symmetry_type), $(bs.ms>0 ? '↑' : '↓')"*brackets[2]
end

export TdSymBasisState

# α, β functions for the hopping parameters A,B,C
α(A,B,C)= (A-B)/C - 1/2 +(1/(2*C))*sqrt(4*(A-B)*(A-B-C)+9*C^2)
β(A,B,C)= (A-B)/C - 1/2 -(1/(2*C))*sqrt(4*(A-B)*(A-B-C)+9*C^2)

export α,β

"""
    getTdSymBasis(basis :: SPBasis{TetramerBasisState{BS}} where {BS<:AbstractSPSSBasisState}, 
        site1::Int64, site2::Int64, site3::Int64, site4::Int64) :: SPBasis{TdSymBasisState{BS}} 
        where {BS<:AbstractSPSSBasisState}

This function provides the pre-implemented single particle - four site delocalized basis in the tetramer symmetry.
"""
function getTdSymBasis(A::Float64, B::Float64, C::Float64,
        site1::Int64, site2::Int64, site3::Int64, site4::Int64) :: SPBasis{TdSymBasisState}
    # make a list of multisite states
    TdSym_states = TdSymBasisState[]
    # push all states for all sites
    
    push!(TdSym_states, TdSymBasisState(:a_1,      1//2, A, B, C, site1, site2, site3, site4))
    push!(TdSym_states, TdSymBasisState(:a_1,     -1//2, A, B, C, site1, site2, site3, site4))
    
    push!(TdSym_states, TdSymBasisState(:e_alpha,  1//2, A, B, C, site1, site2, site3, site4))
    push!(TdSym_states, TdSymBasisState(:e_alpha, -1//2, A, B, C, site1, site2, site3, site4))
    
    push!(TdSym_states, TdSymBasisState(:e_beta,   1//2, A, B, C, site1, site2, site3, site4))
    push!(TdSym_states, TdSymBasisState(:e_beta,  -1//2, A, B, C, site1, site2, site3, site4))
    
    push!(TdSym_states, TdSymBasisState(:t_2alpha, 1//2, A, B, C, site1, site2, site3, site4))
    push!(TdSym_states, TdSymBasisState(:t_2alpha,-1//2, A, B, C, site1, site2, site3, site4))
    
    push!(TdSym_states, TdSymBasisState(:t_2beta,  1//2, A, B, C, site1, site2, site3, site4))
    push!(TdSym_states, TdSymBasisState(:t_2beta, -1//2, A, B, C, site1, site2, site3, site4))
    
    push!(TdSym_states, TdSymBasisState(:t_2gamma, 1//2, A, B, C, site1, site2, site3, site4))
    push!(TdSym_states, TdSymBasisState(:t_2gamma,-1//2, A, B, C, site1, site2, site3, site4))
    
    push!(TdSym_states, TdSymBasisState(:t_balpha, 1//2, A, B, C, site1, site2, site3, site4))
    push!(TdSym_states, TdSymBasisState(:t_balpha,-1//2, A, B, C, site1, site2, site3, site4))
    
    push!(TdSym_states, TdSymBasisState(:t_bbeta,  1//2, A, B, C, site1, site2, site3, site4))
    push!(TdSym_states, TdSymBasisState(:t_bbeta, -1//2, A, B, C, site1, site2, site3, site4))

    push!(TdSym_states, TdSymBasisState(:t_bgamma, 1//2, A, B, C, site1, site2, site3, site4))
    push!(TdSym_states, TdSymBasisState(:t_bgamma,-1//2, A, B, C, site1, site2, site3, site4))
    
    push!(TdSym_states, TdSymBasisState(:t_aalpha, 1//2, A, B, C, site1, site2, site3, site4))
    push!(TdSym_states, TdSymBasisState(:t_aalpha,-1//2, A, B, C, site1, site2, site3, site4))
    
    push!(TdSym_states, TdSymBasisState(:t_abeta,  1//2, A, B, C, site1, site2, site3, site4))
    push!(TdSym_states, TdSymBasisState(:t_abeta, -1//2, A, B, C, site1, site2, site3, site4))
    
    push!(TdSym_states, TdSymBasisState(:t_agamma, 1//2, A, B, C, site1, site2, site3, site4))
    push!(TdSym_states, TdSymBasisState(:t_agamma,-1//2, A, B, C, site1, site2, site3, site4))
        
    # return the multisite basis
    return SPBasis{TdSymBasisState}(TdSym_states)
end

export getTdSymBasis