##############################################################
#   SINGLE PARTICLE SINGLE SITE basis state A1G BASIS
#   A1G BASIS (only for s=1/2, l=1)
#   |orbital, ms>, represented as {Symbol, Rational}
#   with definition of a1g basis in terms of XYZ basis
#   --> |a1g> = 1/sqrt(3) ( |x> + |y> + |z> )
#   --> |eg+> = 1/sqrt(3) ( e^(-im 2pi/3)|x> + e^(+im 2pi/3)|y> + |z> )
#   --> |eg-> =-1/sqrt(3) ( e^(+im 2pi/3)|x> + e^(-im 2pi/3)|y> + |z> )
##############################################################
"""
    BasisStateA1G <: AbstractSPSSBasisState

``\\left| orbital,m_s \\right>`` represented as `{Symbol, Rational{Int64}}`.

With the electron A1G basis defined in terms of the XYZ basis as follows:

``\\left| a_{1g} \\right>= \frac{1}{\\sqrt{3}} ( \\left| x \\right> + \\left| y \\right> + \\left| z \\right> ) ``

``\\left| e_{g+}^\\pi \\right>=+\frac{1}{\\sqrt{3}} ( e^{-i 2\\pi/3}\\left| x \\right> + e^{+i 2\\pi/3}\\left| y \\right>  + \\left| z \\right>) ``

``\\left| e_{g-}^\\pi \\right>=-\frac{1}{\\sqrt{3}} ( e^{+i 2\\pi/3}\\left| x \\right> + e^{-i 2\\pi/3}\\left| y \\right>  + \\left| z \\right>)``

See also D. Khomskii "Transition Metal Compounds", Equation (3.13), and transform to the hole picture (cc).

``\\left| e_{g1} \\right>=-\\frac{i}{\\sqrt{2}} ( \\left| e_{g+} \\right> + \\left| e_{g-} \\right> )``

``\\left| e_{g2} \\right>=+\\frac{1}{\\sqrt{2}} ( \\left| e_{g+} \\right> - \\left| e_{g-} \\right> )``

see again Khomskii's Book, Equation (3.14, 3.15).
"""
struct BasisStateA1G <: AbstractSPSSBasisState
    # orbital part
    orbital :: Symbol
    # spin part
    ms      :: Rational{Int64}
end

function orb_string_a1g(orbital :: Symbol)
    if orbital == :a1g
        return "a1g"
    elseif orbital == :egp
        return "eg+"
    elseif orbital == :egm
        return "eg-"
    elseif orbital == :eg1
        return "eg1"
    elseif orbital == :eg2
        return "eg2"
    else
        return "UNKNOWN A1G ORBITAL $(orbital)"
    end
end

# Custom show function
function Base.show(io::IO, bs::BasisStateA1G)
    bsstr = haskey(io, :compact) ? "" : "A1G"
    print(io, bsstr*"|"*orb_string_a1g(bs.orbital)*","*(bs.ms>0 ? '↑' : '↓')*">")
end
# custom summary function
function summary(bs::BasisStateA1G, brackets="()")
    return brackets[1]*"$(orb_string_a1g(bs.orbital)),$(bs.ms>0 ? '↑' : '↓')"*brackets[2]
end

# export the basis type
export BasisStateA1G



#################################
#   Pre-implemented t2g basis   #
#################################
"""
    getT2GBasisA1G(real_orbitals::Bool=false) :: SPBasis{BasisStateA1G}

This function provides the pre-implemented single particle - single site A1G basis for the t2g.

Note that the function does not return real orbitals by default. In order to get them, enter `true` as the input parameter.
"""
function getT2GBasisA1G(real_orbitals::Bool=false) :: SPBasis{BasisStateA1G}
    # check wether real orbitals are requested
    if real_orbitals
        # reset the basis object in the RIXSSite object
        states = BasisStateA1G[
            # add all terms into the basis
            BasisStateA1G(:a1g, +1//2),
            BasisStateA1G(:a1g, -1//2),
            BasisStateA1G(:eg1, +1//2),
            BasisStateA1G(:eg1, -1//2),
            BasisStateA1G(:eg2, +1//2),
            BasisStateA1G(:eg2, -1//2)
        ]
        # return the basis
        return SPBasis{BasisStateA1G}(states)
    else
        # reset the basis object in the RIXSSite object
        states = BasisStateA1G[
            # add all terms into the basis
            BasisStateA1G(:a1g, +1//2),
            BasisStateA1G(:a1g, -1//2),
            BasisStateA1G(:egp, +1//2),
            BasisStateA1G(:egp, -1//2),
            BasisStateA1G(:egm, +1//2),
            BasisStateA1G(:egm, -1//2)
        ]
        # return the basis
        return SPBasis{BasisStateA1G}(states)
    end
end
export getT2GBasisA1G