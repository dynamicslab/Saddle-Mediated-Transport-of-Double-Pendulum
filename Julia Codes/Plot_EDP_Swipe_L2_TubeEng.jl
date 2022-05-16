""" 
# This file will be sued to plot the L2 Homoclinic Poincare cut at different energy level
# Coded By: KK
# Last Updated: 05/15/2022
"""
##-----------------------------
## First load necessary files and Pacakages
@time include("Functions/LoadPackages.jl")
@time include("Functions/TubeFunctionDP.jl")
@time include("Functions/PlotSettings.jl")

## Next load the necessary data for plotting
@load "Results/EDP_Swipe_L2_Tube_Eng.jld2" StoredCutDataL2 EngListL2

EngListL2=round.(EngListL2[EngListL2.!=0],digits=3)

## Start plotting
for i=1:size(StoredCutDataL2,1)
    # Extract data
    Data=StoredCutDataL2[i,1]
    # Plot the result
    PoincarPlot=Plots.plot(size=(650,400),tickfontsize=20)
    PoincarPlot=Plots.scatter!(PoincarPlot,Data[2,:],Data[4,:],color=col_p,
        label="",markersize=1.5,markerstrokewidth=0,markerstrokealpha=0,framestyle=:box)
    PoincarPlot=Plots.scatter!(PoincarPlot,-Data[2,:],Data[4,:],color=col_o,
        label="",markersize=1.5,markerstrokewidth=0,markerstrokealpha=0,framestyle=:box)
    #display(PoincarPlot)

    ## Save the result
    savefig(PoincarPlot,"Figures/EDP_Swipe_L2_Tube_EngIndex_$(i).pdf")

end

## Plot all of them together
PoincarPlot=Plots.plot(size=(650,400),tickfontsize=20)

for i=1:size(StoredCutDataL2,1)
    # Extract data
    Data=StoredCutDataL2[i,1]
    # Plot the result
    PoincarPlot=Plots.scatter!(PoincarPlot,Data[2,:],Data[4,:],color=col_p,
        label="",markersize=1.5,markerstrokewidth=0,markerstrokealpha=0,framestyle=:box)
    PoincarPlot=Plots.scatter!(PoincarPlot,-Data[2,:],Data[4,:],color=col_o,
        label="",markersize=1.5,markerstrokewidth=0,markerstrokealpha=0,framestyle=:box)
    #display(PoincarPlot)

end

# Save the result
savefig(PoincarPlot,"Figures/EDP_Swipe_L2_Tube_All.pdf")

## Plot all of them together (Keep onlt the half)
PoincarPlot=Plots.plot(size=(650,400),tickfontsize=20)

for i=1:size(StoredCutDataL2,1)
    # Extract data
    Data=StoredCutDataL2[i,1]
    # Plot the result
    PoincarPlot=Plots.scatter!(PoincarPlot,Data[2,:],Data[4,:],color=col_p,
        label="",markersize=1,markerstrokewidth=0,markerstrokealpha=0,framestyle=:box)
    #display(PoincarPlot)

end

# Save the result
savefig(PoincarPlot,"Figures/EDP_Swipe_L2_Tube_Half_All.pdf")


