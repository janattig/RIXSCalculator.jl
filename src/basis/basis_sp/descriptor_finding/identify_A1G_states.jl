#################################
# DESCRIPTOR FINDING A1G STATES #
#################################

function getSPBasisState(
            :: Type{BasisStateA1G},
            descriptor :: Tuple{Union{String,Symbol}, Real}
        )

    # find out the orbital
    orb = "NONE"
    orb_readin = lowercase(string(descriptor[1]))
    if orb_readin in ["a1g", "egp", "egm", "eg1", "eg2"]
        orb = orb_readin
    elseif orb_readin == "eg+"
        orb = "egp"
    elseif orb_readin == "eg-"
        orb = "egm"
    end

    # find out the spin
    spn = 0//1
    if abs(descriptor[2]) == 0.5
        spn = sign(descriptor[2]) * 1//2
    end

    # check if finding the orbital and spin was successful
    if orb == "NONE" || spn == 0//1
        error("Could not process input for BasisStateA1G: $(descriptor)")
    end

    # return the state
    return BasisStateA1G(Symbol(orb),spn)
end
function getSPBasisState(
            :: Type{BasisStateA1G},
            descriptor :: Tuple{Union{String,Symbol}, Union{String,Symbol}}
        )

    # find out the spin
    if lowercase(string(descriptor[2])) == "down"
        return getSPBasisState(BasisStateA1G, (descriptor[1], -1/2))
    elseif lowercase(string(descriptor[2])) == "up"
        return getSPBasisState(BasisStateA1G, (descriptor[1], 1/2))
    end
end
function getSPBasisState(
            :: Type{BasisStateA1G},
            descriptor :: Tuple{Any, Any}
        )
    error("Could not process input for BasisStateA1G: $(descriptor)")
end
export getSPBasisState