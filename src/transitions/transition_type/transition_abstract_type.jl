################################################################################
#
#   DEFINITION OF TRANSITION ABSTRACT TYPES
#
################################################################################

# transition
abstract type AbstractTransition end


# export the abstract types
export AbstractTransition



################################################################################
#
#   INTERFACING ABSTRACT TYPES
#
################################################################################

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