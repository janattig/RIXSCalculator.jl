#################################################################################
#   SINGLE PARTICLE SINGLE SITE basis state LS BASIS
#   |l,ml, s,ms>, represented as {Int64, Int64, Rational{Int64}, Rational{Int64}}
#################################################################################
"""
    BasisStateLS <: AbstractSPSSBasisState

``\\left| l,m_l,s,m_s \\right>`` represented as `{Int64, Int64, Rational{Int64}, Rational{Int64}}`.
"""
struct BasisStateLS <: AbstractSPSSBasisState
    # l part
    l   :: Int64
    ml  :: Int64
    # s part
    s   :: Rational{Int64}
    ms  :: Rational{Int64}
end

# short hand notation for l=1, s=1/2
function BasisStateLS(ml :: Integer, ms :: Rational)
    @assert ml in [-1,0,1] && ms in [-1//2,1//2] "ml=$(ml) and ms=$(ms) still consistent with l=1, s=1//2"
    return BasisStateLS(1,ml,1//2,ms)
end



# Custom show function
import Base.show
function Base.show(io::IO, bs::BasisStateLS)
    bsstr = haskey(io, :compact) ? "" : "LS"
    if bs.l == 1 && bs.s == 1//2
        print(io, bsstr*"|"*(bs.ml>0 ? "+$(bs.ml)" : (bs.ml==0 ? " 0" : "$(bs.ml)"))*","*(bs.ms>0 ? '↑' : '↓')*">")
    else
        print(io, bsstr*"|"*string(bs.l)*","*string(bs.ml)*", "*string(bs.s)*","*string(bs.ms)*">")
    end
end

# custom summary function
function summary(bs::BasisStateLS, brackets="()")
    return brackets[1]*"ml=$(bs.ml),$(bs.ms>0 ? '↑' : '↓')"*brackets[2]
end


# export the basis type
export BasisStateLS


#################################
#   Pre-implemented t2g bases   #
#################################
"""
    getT2GBasisLS() :: SPBasis{BasisStateLS}

This function provides the pre-implemented single particle - single site LS basis for the t2g.
"""
function getT2GBasisLS() :: SPBasis{BasisStateLS}
    # construct the list of states
    states = BasisStateLS[
        BasisStateLS(1,+1, 1//2,+1//2),
        BasisStateLS(1,+1, 1//2,-1//2),
        BasisStateLS(1, 0, 1//2,+1//2),
        BasisStateLS(1, 0, 1//2,-1//2),
        BasisStateLS(1,-1, 1//2,+1//2),
        BasisStateLS(1,-1, 1//2,-1//2)
    ]
    # return the basis
    return SPBasis{BasisStateLS}(states)
end
export getT2GBasisLS

