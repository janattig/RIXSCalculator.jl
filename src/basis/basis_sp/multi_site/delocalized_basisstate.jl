struct DelocalizedBasisStateXYZ <: AbstractSPBasisState
    # type of bond
    bonding_type :: Symbol
    # orbital
    orbital :: Symbol
    # ms
    ms :: Rational{Int64}
end

# Custom show function
import Base.show
function Base.show(io::IO, bs::DelocalizedBasisStateXYZ)
    bsstr = haskey(io, :compact) ? "" : "Delocalized"
    print(io, bsstr*"|"*string(bs.bonding_type)*","*string(bs.orbital)*","*(bs.ms>0 ? '↑' : '↓')*">")
end
# custom summary function
function summary(bs::DelocalizedBasisStateXYZ, brackets="()")
    return brackets[1]*"$(bs.bonding_type),$(bs.orbital),$(bs.ms>0 ? '↑' : '↓')"*brackets[2]
end


# export the basis type
export DelocalizedBasisStateXYZ



#########################################
#   Pre-implemented delocalized basis   #
#########################################


"""
    getT2GDelocalizedBasisDelocalizedBasisStateXYZ() :: SPBasis{DelocalizedBasisStateXYZ}

This function provides the pre-implemented single particle - two site bonding and antibonding basis for the t2g.
"""
function getT2GDelocalizedBasis() :: SPBasis{DelocalizedBasisStateXYZ}
    # reset the basis object in the RIXSSite object
    states = DelocalizedBasisStateXYZ[
        # add all terms into the basis
        DelocalizedBasisStateXYZ(:bonding,:x, +1//2),
        DelocalizedBasisStateXYZ(:antibonding,:x, +1//2),
        DelocalizedBasisStateXYZ(:bonding,:y, +1//2),
        DelocalizedBasisStateXYZ(:antibonding,:y, +1//2),
        DelocalizedBasisStateXYZ(:bonding,:z, +1//2),
        DelocalizedBasisStateXYZ(:antibonding,:z, +1//2),
        
        DelocalizedBasisStateXYZ(:bonding,:x, -1//2),
        DelocalizedBasisStateXYZ(:antibonding,:x, -1//2),
        DelocalizedBasisStateXYZ(:bonding,:y, -1//2),
        DelocalizedBasisStateXYZ(:antibonding,:y, -1//2),
        DelocalizedBasisStateXYZ(:bonding,:z, -1//2),
        DelocalizedBasisStateXYZ(:antibonding,:z, -1//2),
    ]
    # return the basis
    return SPBasis{DelocalizedBasisStateXYZ}(states)
end

export getT2GDelocalizedBasis