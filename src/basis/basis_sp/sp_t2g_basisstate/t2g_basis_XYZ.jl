##############################################################
#   SINGLE PARTICLE SINGLE SITE basis state XYZ BASIS
#   XYZ BASIS (only for s=1/2, l=1)
#   |XYZ, ms>, represented as {Symbol, Rational}
#   with definition of XYZ basis in terms of LS basis
#   --> |x> = -1/sqrt(2) ( |1> - |-1> )
#   --> |y> = im/sqrt(2) ( |1> + |-1> )
#   --> |z> = |0>
##############################################################
"""
    BasisStateXYZ <: AbstractSPSSBasisState

``\\left| orbital,m_s \\right>`` represented as `{Symbol, Rational{Int64}}`.

The XYZ basis is defined in terms of the LS basis as follows:

`` \\left| x \\right> = -\\frac{1}{\\sqrt{2}} ( \\left| 1 \\right> - \\left| -1 \\right> ) ``

`` \\left| y \\right> = +\\frac{i}{\\sqrt{2}} ( \\left| 1 \\right> + \\left| -1 \\right> ) ``

`` \\left| z \\right> = \\left| 0 \\right> ``

consistent with L. Ament and G. Khaliullin (Phys. Rev. B 84, 020403(R) (2011)).
"""
struct BasisStateXYZ <: AbstractSPSSBasisState
    # orbital part
    orbital :: Symbol
    # spin part
    ms      :: Rational{Int64}
end

# Custom show function
import Base.show
function Base.show(io::IO, bs::BasisStateXYZ)
    bsstr = haskey(io, :compact) ? "" : "XYZ"
    print(io, bsstr*"|"*string(bs.orbital)*","*(bs.ms>0 ? '↑' : '↓')*">")
end
# custom summary function
function summary(bs::BasisStateXYZ, brackets="()")
    return brackets[1]*"$(bs.orbital),$(bs.ms>0 ? '↑' : '↓')"*brackets[2]
end


# export the basis type
export BasisStateXYZ


#################################
#   Pre-implemented t2g basis   #
#################################
"""
    getT2GBasisXYZ() :: SPBasis{BasisStateXYZ}

This function provides the pre-implemented single particle - single site XYZ basis for the t2g.
"""
function getT2GBasisXYZ() :: SPBasis{BasisStateXYZ}
    # reset the basis object in the RIXSSite object
    states = BasisStateXYZ[
        # add all terms into the basis
        BasisStateXYZ(:x, +1//2),
        BasisStateXYZ(:x, -1//2),
        BasisStateXYZ(:y, +1//2),
        BasisStateXYZ(:y, -1//2),
        BasisStateXYZ(:z, +1//2),
        BasisStateXYZ(:z, -1//2)
    ]
    # return the basis
    return SPBasis{BasisStateXYZ}(states)
end
export getT2GBasisXYZ











# t2g LS - XYZ conversion

# The t2g LS basis is related to the t2g XYZ basis as follows:

# ``\\left| +1 \\right>=-\\frac{1}{\\sqrt{2}} (\\left|x\\right> + i \\left| y\\right>)``

# ``\\left| -1 \\right>=+\\frac{1}{\\sqrt{2}} (\\left|x\\right> - i \\left| y\\right>)``

# ``\\left| 0 \\right>=\\left| z\\right>``

# consistent with L. Ament and G. Khaliullin (Phys. Rev. B 84, 020403(R) (2011)).