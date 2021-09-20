# abstract SP basis state types
include("sp_basisstate_abstract_type.jl")
# abstract SP basis types
include("sp_basis_type.jl")


# include all t2g bases definitions
include("sp_t2g_basisstate/t2g_basis_LS.jl")
include("sp_t2g_basisstate/t2g_basis_J.jl")
include("sp_t2g_basisstate/t2g_basis_XYZ.jl")
include("sp_t2g_basisstate/t2g_basis_A1G.jl")



# include single particle - multi site description
include("sp_multi_site/multi_site_basisstate.jl")
include("sp_multi_site/multi_site_functions.jl")

# include composite basis state definitions and functions
include("sp_composite_basisstate/single_site_composite_basisstate.jl")
include("sp_composite_basisstate/multi_site_composite_basisstate.jl")


# include overlap definitions
include("sp_overlaps.jl")


# include the descriptor finding functions
include("descriptor_finding/identify_LS_states.jl")
include("descriptor_finding/identify_XYZ_states.jl")
include("descriptor_finding/identify_A1G_states.jl")

include("descriptor_finding/identify_composite_states.jl")
include("descriptor_finding/identify_multi_site_states.jl")

include("descriptor_finding/functions.jl")
