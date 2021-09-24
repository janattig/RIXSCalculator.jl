# abstract type
include("operators_mp_abstract_type.jl")


# 1 particle
include("1p_abstract_type.jl")
include("1p_generalized_sp_type_definition.jl")


# 2 particles
include("2p_abstract_type.jl")
include("2p_densitydensity_type_definition.jl")
include("2p_densitydensity_generate_terms.jl")
include("2p_scattering_type_definition.jl")
include("2p_scattering_generate_terms.jl")


# N particles
include("Np_interaction_hamiltonian_abstract_type.jl")
