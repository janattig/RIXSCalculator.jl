#################################
# DESCRIPTOR FINDING XYZ STATES #
#################################

function getSPBasisState(
            :: Type{BasisStateXYZ},
            descriptor :: Tuple{Union{String,Symbol}, Real}
        )

    # find out the orbital
    orb = "NONE"
    orb_readin = lowercase(string(descriptor[1]))
    if orb_readin in ["x","y","z"]
        orb = orb_readin
    end

    # find out the spin
    spn = 0//1
    if abs(descriptor[2]) == 0.5
        spn = sign(descriptor[2]) * 1//2
    end

    # check if finding the orbital and spin was successful
    if orb == "NONE" || spn == 0//1
        error("Could not process input for BasisStateXYZ: $(descriptor)")
    end

    # return the state
    return BasisStateXYZ(Symbol(orb),spn)
end
function getSPBasisState(
            :: Type{BasisStateXYZ},
            descriptor :: Tuple{Union{String,Symbol}, Union{String,Symbol}}
        )

    # find out the spin
    if lowercase(string(descriptor[2])) == "down"
        return getSPBasisState(BasisStateXYZ, (descriptor[1], -1/2))
    elseif lowercase(string(descriptor[2])) == "up"
        return getSPBasisState(BasisStateXYZ, (descriptor[1], 1/2))
    end
end
function getSPBasisState(
            :: Type{BasisStateXYZ},
            descriptor :: Tuple{Any, Any}
        )
    error("Could not process input for BasisStateXYZ: $(descriptor)")
end
export getSPBasisState