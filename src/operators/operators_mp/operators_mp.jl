##############################################################
#
#   MULTI PARTICLE OPERATOR DEFINITIONS
#   FILE STRUCTURE
#
#   1) MP operators of SP operators
#   - type definition
#   - interface functions
#   - convenience functions
#
#   2) MP operators for density density interactions
#   - type definition
#   - interface functions
#   - convenience functions
#
#   3) MP operators for two particle scattering
#   - type definition
#   - interface functions
#   - convenience functions
#
##############################################################





# ABSTRACT TYPE DEFINITION
include("operators_mp_abstract_type.jl")





##############################################################
#
#   1) MP operators of SP operators (Single particle)
#   - abstract type definition
#   - type definition
#   - interface functions
#   - convenience functions
#
##############################################################


# everything else included in subfile
include("operators_mp_1p_general_sp.jl")




##############################################################
#
#   2) MP operators for density density interactions (2 particle)
#   - type definition
#   - interface functions
#   - convenience functions
#
##############################################################


# included in subfile
include("operators_mp_2p_density.jl")




##############################################################
#
#   3) MP operators for two particle scattering (2 particle)
#   - type definition
#   - interface functions
#   - convenience functions
#
##############################################################

# included in subfile
include("operators_mp_2p_scattering.jl")






##############################################################
#
#   4) Interaction Hamiltonians (on-site)
#   - abstract type definition
#   - type definition Perkins Woelfle Hamiltonian
#   - interface functions
#   - convenience functions
#
##############################################################


# included in subfile
include("operators_mp_interaction_ham.jl")
