##############################################################
#   SINGLE PARTICLE SINGLE SITE basis state J BASIS
#   J BASIS for l=1, s=1//2,
#   |j,mj>, represented as {Rational, Rational}
##############################################################
"""
    BasisStateJ <: AbstractSPSSBasisState

``\\left| j,m_j \\right>`` represented as `{Rational{Int64}, Rational{Int64}}`.

The J basis is defined through the LS basis as follows:

``\\left| \\frac{1}{2},-\\frac{1}{2} \\right>=\\frac{1}{\\sqrt{3}}  \\left|0,\\downarrow\\right> - \\sqrt{\\frac{2}{3}} \\left| -1,\\uparrow\\right> ``

``\\left| \\frac{1}{2},+\\frac{1}{2} \\right>=\\frac{1}{\\sqrt{3}}  \\left|0,\\uparrow\\right> - \\sqrt{\\frac{2}{3}} \\left| +1,\\downarrow\\right> ``

``\\left| \\frac{3}{2},-\\frac{3}{2} \\right>=\\left| -1,\\downarrow \\right>``

``\\left| \\frac{3}{2},-\\frac{1}{2} \\right>=\\sqrt{\\frac{2}{3}} \\left|0,\\downarrow\\right> + \\sqrt{\\frac{1}{3}} \\left| -1,\\uparrow\\right>``

``\\left| \\frac{3}{2},+\\frac{1}{2} \\right>=\\sqrt{\\frac{2}{3}} \\left|0,\\uparrow\\right> + \\sqrt{\\frac{1}{3}} \\left| +1,\\downarrow\\right>``

``\\left| \\frac{3}{2},+\\frac{3}{2} \\right>=\\left| +1, \\uparrow \\right>``
"""
struct BasisStateJ <: AbstractSPSSBasisState
    # j part
    j   :: Rational{Int64}
    mj  :: Rational{Int64}
end

# Custom show function
import Base.show
function Base.show(io::IO, bs::BasisStateJ)
    bsstr = haskey(io, :compact) ? "" : "J"
    print(io, bsstr*"|"*string(bs.j.num)*"/"*string(bs.j.den)*","*(bs.mj==0 ? " " : (bs.mj>0 ? "+" : "-"))*string(abs(bs.mj.num))*"/"*string(abs(bs.mj.den))*"⟩")
end

# custom summary function
function summary(bs::BasisStateJ, brackets="()")
    return brackets[1]*"j=$(bs.j.num)/$(bs.j.den),mj=$(bs.mj>0 ? '+' : '-')$(abs(bs.mj.num))/$(bs.mj.den)"*brackets[2]
end


# export the basis type
export BasisStateJ


#################################
#   Pre-implemented t2g basis   #
#################################
"""
    getT2GBasisJ() :: SPBasis{BasisStateJ}

This function provides the pre-implemented single particle - single site J basis for the t2g.
"""
function getT2GBasisJ() :: SPBasis{BasisStateJ}
    # reset the basis object in the RIXSSite object
    states = BasisStateJ[
        BasisStateJ(1//2, +1//2),
        BasisStateJ(1//2, -1//2),
        BasisStateJ(3//2, +3//2),
        BasisStateJ(3//2, +1//2),
        BasisStateJ(3//2, -1//2),
        BasisStateJ(3//2, -3//2)
    ]
    # return the basis
    return SPBasis{BasisStateJ}(states)
end
export getT2GBasisJ