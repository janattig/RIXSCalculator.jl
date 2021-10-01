################################################################################
#
#   GENERAL MODULE FOR RIXS CALCULATIONS
#   AUTHOR Jan Attig
#          Luca Peterlini
#   JULIA v.1+
#
################################################################################



# start of module
module RIXSCalculator


    # GENERAL USINGS
    using LinearAlgebra



    # INCLUDE EVERTHING FROM SUBFILES
    include("coordinate_frames/coordinate_frames.jl")
    include("basis/basis.jl")
    include("operators/operators.jl")
    #include("transitions/transitions.jl")


# end of module
end
