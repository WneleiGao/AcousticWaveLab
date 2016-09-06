using PyCall
@pyimport matplotlib.pyplot as plt
@pyimport matplotlib.lines as lines
@pyimport matplotlib.animation as anim

export sptPlot,
       SeisPlot,
       waveAnim

include("sptPlot.jl")
include("SeisPlot.jl")
include("waveAnim.jl")
