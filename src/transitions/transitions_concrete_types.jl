################################################################################
#
#   DEFINITION OF SPECTRUM CONCRETE TYPE
#
################################################################################


# TYPE DEFINITION OF A SPECTRUM
"""
    Spectrum{T} <: AbstractSpectrum{T}

This object is a mutable struct defined by the list of transitions of the spectrum,  `transitions :: Vector{T}`.
"""
mutable struct Spectrum{T} <: AbstractSpectrum{T}

    # list of transitions
    transitions :: Vector{T}

end

# export the spectrum type
export Spectrum



# IMPLEMENTATION OF INTERFACE

# intensity
function intensity(spectrum :: Spectrum{T}, energy :: Real) :: Float64 where {T<:AbstractTransition}
    # add up weights of all transitions
    intensity = 0.0 :: Float64
    # iterate over all transitions
    for t in spectrum.transitions
        intensity += abs(energy - frequency(t)) < linewidth(t)*15 ? weight(t) / (sqrt(2*pi) * linewidth(t)) * exp(-(energy - frequency(t))^2 / (2*linewidth(t)^2)) : 0.0
    end
    # return value
    return intensity
end


# Function to add two spectra
function +(spectrum_1::Spectrum{T}, spectrum_2::Spectrum{T}) :: Spectrum{T} where {T<:AbstractTransition}
    # return a new spectrum object
    return Spectrum{T}(Base.vcat(spectrum_1.transitions, spectrum_2.transitions))
end