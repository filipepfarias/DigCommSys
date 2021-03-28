# __precompile__()

module DigCommSys
    import Statistics
    import Plots
    import SpecialFunctions
    import Random

    include("GrayCode.jl")
    export GrayCode

    include("QAM.jl")
    export MQAM

    include("PSK.jl")
    export MPSK

    include("bitModulation.jl")
    export bitModulation

    include("bitDemodulation.jl")
    export bitDemodulation

    include("signalClassifier.jl")
    export absClassifier

    include("AWGN.jl")
    export AWGN

    include("BERSER.jl")
    export BER, SER, TER
end
