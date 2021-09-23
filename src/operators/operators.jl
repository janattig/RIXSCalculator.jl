# IMPORTS FOR LATER USE
# function for inheritance of AbstractArray
# import Base.size
# import Base.getindex
# import Base.setindex
# import Base.IndexStyle


# using Linear algebra package
using LinearAlgebra
using Combinatorics


include("operators_abstract_type.jl")

################################################################################
#
#   DEFINITION OF SP OPERATORS AND FUNCTIONS
#
################################################################################


# included in subfile
include("operators_sp/operators_sp.jl")





################################################################################
#
#   DEFINITION OF MP OPERATORS AND FUNCTIONS
#
################################################################################

# included in subfile
include("operators_mp/operators_mp.jl")







################################################################################
#
#   DEFINITION OF SUMMING TWO OPERATORS
#
################################################################################

# summing two single site operators, included in subfile
include("operators_sum.jl")

# scalar factor in front of operator
include("operators_factor.jl")

# zero element of operators
include("operators_zero.jl")




################################################################################
#
#   DEFINITION OF SP OPERATORS THAT TRANSFORM TO DIFFERENT BASES
#
################################################################################

# included in subfile
include("operators_projectors.jl")









################################################################################
#
#   DEFINITION OF UTILITY FUNCTIONS FOR EIGENSYSTEMS
#
################################################################################

# included in subfile
include("operators_eigensystem.jl")
