# abstract SP basis state types
include("basisstate_abstract_type_sp.jl")
include("basisstate_abstract_type_sp_single_site.jl")
# abstract SP basis types
include("basis_sp_type_definition.jl")


# include all t2g bases definitions
include("basisstate_types_for_t2g/t2g_basis_LS.jl")
include("basisstate_types_for_t2g/t2g_basis_J.jl")
include("basisstate_types_for_t2g/t2g_basis_XYZ.jl")
include("basisstate_types_for_t2g/t2g_basis_A1G.jl")



# include single particle - multi site description
include("multi_site/multi_site_basisstate.jl")
include("multi_site/multi_site_functions.jl")

# include composite basis state definitions and functions
include("composite_basisstate/single_site_composite_basisstate.jl")
include("composite_basisstate/multi_site_composite_basisstate.jl")


# include overlap definitions
include("sp_overlaps.jl")


# include the descriptor finding functions
include("descriptor_finding/identify_LS_states.jl")
include("descriptor_finding/identify_XYZ_states.jl")
include("descriptor_finding/identify_A1G_states.jl")

include("descriptor_finding/identify_composite_states.jl")
include("descriptor_finding/identify_multi_site_states.jl")

include("descriptor_finding/functions.jl")
