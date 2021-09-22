##############################################################
#
#   SINGLE PARTICLE OPERATORS
#   FILE STRUCTURE
#
#   1) basic operators for SP SS (in bare numbers)
#
#   2) Definition of concrete SP SS operators
#   - spin orbit LS operator
#   - distortion Ln^2 operator
#   - magnetic field BS operator
#
#   3) basic operators for SP MS
#
##############################################################


##############################################################
#
#   1) basic operators for SP SS (in bare numbers)
#   (Definition of delta function as well as Lx, Ly, Lz, etc.)
#
##############################################################

# included in subfile
include("operators_sp_basic.jl")



##############################################################
#
#   2) Definition of concrete SP SS operators
#   - spin orbit LS operator
#   - distortion Ln^2 operator
#   - magnetic field BS operator
#
##############################################################

# spin orbit operator in subfile
include("../operators_specific/singleparticle_singlesite/operators_sp_ss_spin_orbit.jl")

# distortion operator in subfile
include("../operators_specific/singleparticle_singlesite/operators_sp_ss_distortion.jl")

# magnetic field operator in subfile
include("../operators_specific/singleparticle_singlesite/operators_sp_ss_magnetic_field.jl")




##############################################################
#
#   3) basic operators for SP MS
#
##############################################################

# included in subfile
include("../operators_specific/singleparticle_multisite/operators_sp_ms.jl")





##############################################################
#
#   4) basic operators for SP hopping
#
##############################################################

# included in subfile
include("../operators_specific/operators_sp_hopping.jl")
