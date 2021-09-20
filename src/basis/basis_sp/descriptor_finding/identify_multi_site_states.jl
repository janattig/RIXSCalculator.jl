################################
# DESCRIPTOR FINDING MS STATES #
################################

function getSPBasisState(
            :: Union{Type{SPMSBasisState{SPSS}},Type{SPSS}},
            descriptor :: Union{
                    Tuple{Integer, String, Real},
                    Tuple{Integer, Symbol, Real},
                    Tuple{Integer, Real, Real},
                    Tuple{Integer, String, Symbol},
                    Tuple{Integer, Symbol, Symbol},
                    Tuple{Integer, Real, Symbol},
                    Tuple{Integer, String, String},
                    Tuple{Integer, Symbol, String},
                    Tuple{Integer, Real, String}
                }
        ) where {SPSS <: AbstractSPSSBasisState}
    return SPMSBasisState{SPSS}(getSPBasisState(SPSS, (descriptor[2], descriptor[3])), descriptor[1])
end
export getSPBasisState