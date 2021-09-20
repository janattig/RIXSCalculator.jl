################################
# DESCRIPTOR FINDING LS STATES #
################################
function getSPBasisState(
            :: Type{BasisStateLS},
            descriptor :: Tuple{Real, Real}
        )

    # find out the ml value
    ml = round(Int64, descriptor[1])

    # find out the spin
    spn = 0//1
    if abs(descriptor[2]) == 0.5
        spn = 1//2 * sign(descriptor[2])
    end

    # check if finding the orbital and spin was successful
    if spn == 0//1 || !(ml in [-1,0,1])
        error("Could not process input for BasisStateLS: $(descriptor)")
    end

    # return the state
    return BasisStateLS(1, ml, 1//2, spn)
end
function getSPBasisState(
            :: Type{BasisStateLS},
            descriptor :: Tuple{Real, Union{String,Symbol}}
        )

    # find out the spin
    if lowercase(string(descriptor[2])) == "down"
        return getSPBasisState(BasisStateLS, (descriptor[1], -1/2))
    elseif lowercase(string(descriptor[2])) == "up"
        return getSPBasisState(BasisStateLS, (descriptor[1], 1/2))
    end
end
function getSPBasisState(
            :: Type{BasisStateLS},
            descriptor :: Tuple{Any, Any}
        )
    error("Could not process input for BasisStateLS: $(descriptor)")
end
export getSPBasisState