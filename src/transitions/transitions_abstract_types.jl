################################################################################
#
#   DEFINITION OF SPECTRUM ABSTRACT TYPES
#
################################################################################

# transition
abstract type AbstractTransition end

# spectrum
abstract type AbstractSpectrum{T<:AbstractTransition} end


# export the abstract types
export AbstractTransition, AbstractSpectrum



################################################################################
#
#   INTERFACING ABSTRACT TYPES
#
################################################################################

# TRANSITION

# get the transition weight
function weight(transition :: T) :: Float64 where {T<:AbstractTransition}
    error("not implemented function 'weight' for transition type " * T)
end

# get the transition linewidth
function linewidth(transition :: T) :: Float64 where {T<:AbstractTransition}
    error("not implemented function 'linewidth' for transition type " * T)
end

# get the transition frequency
function frequency(transition :: T) :: Float64 where {T<:AbstractTransition}
    error("not implemented function 'frequency' for transition type " * T)
end


# export the interface
export weight, linewidth, frequency




# SPECTRUM

# get the intensity at a certain energy
function intensity(spectrum :: S, energy :: Real) :: Float64 where {S<:AbstractSpectrum{T} where T<:AbstractTransition}
    error("not implemented function 'intensity' for spectrum type " * S)
end


# export the interface
export intensity