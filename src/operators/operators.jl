# using packages
using LinearAlgebra
using Combinatorics



# abstract type
include("operator_abstract_type.jl")



# single particle operators
include("operators_sp/operators_sp.jl")

# multi particle operators
include("operators_mp/operators_mp.jl")



# projectors
include("projector_operators/operators_projectors.jl")



# specific operators
include("specific_operators/fundamental_operators.jl")

include("specific_operators/spin_orbit_coupling.jl")
include("specific_operators/magnetic_field.jl")
include("specific_operators/crystal_field_distortion.jl")

include("specific_operators/hopping.jl")

include("specific_operators/hubbard_interaction.jl")

include("specific_operators/dipole_operator.jl")

include("specific_operators/Jz_operator.jl")




# math functions for operators
include("operator_math/operator_math.jl")



# eigensystem of an operator
include("operators_eigensystem.jl")
