################################################################################
#
#   DEFINITION OF SPECTRUM ABSTRACT TYPES
#
################################################################################


# spectrum
abstract type AbstractSpectrum{T<:AbstractTransition} end

# export the abstract type
export AbstractSpectrum



################################################################################
#
#   INTERFACING ABSTRACT TYPES
#
################################################################################

# get the intensity at a certain energy
function intensity(spectrum :: S, energy :: Real) :: Float64 where {S<:AbstractSpectrum{T} where T<:AbstractTransition}
    error("not implemented function 'intensity' for spectrum type " * S)
end


# export the interface
export intensity