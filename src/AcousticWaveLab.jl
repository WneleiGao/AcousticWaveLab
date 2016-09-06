module AcousticWaveLab
    using Requires
    include("FDMtx/FDMtx.jl")
    include("plot/plotting.jl")
    include("imaging/imaging.jl")
    include("sourceLocation/sourceLocation.jl")
    @require PyPlot include("plot/plotting.jl")
end
