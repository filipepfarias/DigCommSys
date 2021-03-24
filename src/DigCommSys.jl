# __precompile__()

module DigCommSys
    import Statistics
    import Plots
    import SpecialFunctions

    include("GrayCode.jl")
    export GrayCode

    include("QAM.jl")
    export MQAM

    include("bitModulation.jl")
    export bitModulation

    include("bitDemodulation.jl")
    export bitDemodulation
end
