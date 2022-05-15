"""
# This file is used to plot the L1 and L2 PO that is used to generate the tubes.
# Coded By: KK
# Last Updated: 05/15/2022
"""
##-----------------------------
## First load necessary files and Pacakages
@time include("Functions/LoadPackages.jl")
@time include("Functions/TubeFunctionTBP.jl")
@time include("Functions/PlotSettings.jl")

# Load the PO data
μ=9.537e-4;EngTarget=-1.5175;

@save "Results/L1PO_μ_$(μ)_Eng_$(EngTarget)_Tube.jld2" μ par x_inc_L1 PO1_Time PO1_Data EqP1 EqP2 EqP3 EngP1 EngP2 EngP3
@save "Results/L2PO_μ_$(μ)_Eng_$(EngTarget)_Tube.jld2" μ par x_inc_L2 PO2_Time PO2_Data EqP1 EqP2 EqP3 EngP1 EngP2 EngP3

## Plot the PO
# Define the size of marker and linewidth
mks1=2.5
ls1=3

pp1=Plots.scatter([EqP1,EqP2],[0,0],markersize=2,label="",size=resol1,tickfontsize=tickfontsize,
    marker = (:circle, mks1, 0.99, :blue, stroke(1.0, 0.9, :black, :solid)))
pp1=Plots.plot!(PO1_Data[1,:],PO1_Data[2,:],linewidth=ls1,color=col_p,label="",
    size=resol1,tickfontsize=tickfontsize,framestyle=:box)
pp1=Plots.plot!(PO2_Data[1,:],PO2_Data[2,:],linewidth=ls1,color=col_p,label="",
    size=resol1,tickfontsize=tickfontsize,framestyle=:box)

# Save it
#savefig(pp1,"Figures/TBP_TubePOs.pdf")

display(pp1)


