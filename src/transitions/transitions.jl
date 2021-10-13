##############################################################
#
#   TRANSITIONS
#
##############################################################


# DEFINITION OF TRANSITION
include("transition_type/transition_abstract_type.jl")
include("transition_type/transition_concrete_type.jl")

# DEFINITION OF SPECTRUM
include("spectrum_type/spectrum_abstract_type.jl")
include("spectrum_type/spectrum_concrete_type.jl")

# INCLUDE VARIOUS FUNCTIONS
include("transitions_calculations.jl")
include("transitions_auxiliary_functions.jl")