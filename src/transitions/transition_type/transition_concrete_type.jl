################################################################################
#
#   DEFINITION OF TRANSITION CONCRETE TYPE
#
################################################################################

# DEFINITION OF STRUCT
"""
    mutable struct Transition <: AbstractTransition

This object defines a Transition.

# Fields

- `weight :: Float64`
- `linewidth :: Float64`
- `frequency :: Float64`

"""
mutable struct Transition <: AbstractTransition

    # weight
    weight      :: Float64

    # linewidth
    linewidth   :: Float64

    # frequency
    frequency   :: Float64

end

# export the struct
export Transition




# IMPLEMENTATION OF INTERFACE

# get the transition weight
function weight(transition :: Transition) :: Float64
    return transition.weight
end

# get the transition linewidth
function linewidth(transition :: Transition) :: Float64
    return transition.linewidth
end

# get the transition frequency
function frequency(transition :: Transition) :: Float64
    return transition.frequency
end