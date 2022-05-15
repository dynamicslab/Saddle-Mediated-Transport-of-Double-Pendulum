"""
# This file is used to plot the L1 Homoclinic orbits
# Coded By: KK
# Last Updated: 05/15/2022
"""
##-----------------------------
## First load necessary files and Pacakages
@time include("Functions/LoadPackages.jl")
@time include("Functions/TubeFunctionTBP.jl")
@time include("Functions/PlotSettings.jl")

## Load the Results
μ=9.537e-4;Gap=10;EngTarget=-1.5175;

@load "Results/TBP_Manifold_$(μ)_Gap_$(Gap)_Eng_$(EngTarget).jld2" PODataPO1 PODataPO2 μ par Traj1PO1 Traj2PO1 Traj3PO1 Traj4PO1 Traj1PO2 Traj2PO2 Traj3PO2 Traj4PO2

## Get the Poincare cut data
PDataL1_unstable_pos=GetPoincareData(Traj1PO1)
PDataL1_unstable_neg=GetPoincareData(Traj2PO1)
PDataL1_stable_pos=GetPoincareData(Traj3PO1)
PDataL1_stable_neg=GetPoincareData(Traj4PO1)

## Calculate the index of trajectory
L1IndexUnstable=findmax(Int64.(abs.(PDataL1_unstable_neg[3,:]).<=2e-4))
L1IndexStable=findmax(Int64.(abs.(PDataL1_stable_pos[3,:]).<=2e-4))

# Get the corresponding trajectory
HomoL1Data1=Array(Traj2PO1[L1IndexUnstable[1]])
HomoL1Data2=Array(Traj3PO1[L1IndexStable[1]])

# Now plot the symetric Homiclinic orbits
# HomoPlotL1=Scene(resolution=(1000,1000))
# HomoPlotL1=Makie.lines!(HomoPlotL1,HomoL1Data1[1,:],HomoL1Data1[2,:],color=col_p,linewidth=4)
# HomoPlotL1=Makie.lines!(HomoPlotL1,HomoL1Data2[1,:],HomoL1Data2[2,:],color=col_o,linewidth=4)
ls=0.5

HomoPlotL1=Plots.plot(size=(400,400),tickfontsize=tickfontsize,framestyle=:box)
HomoPlotL1=Plots.plot!(HomoL1Data1[1,:],HomoL1Data1[2,:],color=col_p,linewidth=3*ls,label="")
HomoPlotL1=Plots.plot!(HomoL1Data2[1,:],HomoL1Data2[2,:],color=col_o,linewidth=3*ls,label="")

# Now plot the unsymetric L1 Homoclinic orbits
Dummy=findmax(Int64.(abs.(PDataL1_unstable_neg[3,:]-PDataL1_stable_pos[3,:]).<1e-2))

HomoL1Data1=Array(Traj2PO1[Dummy[2]])
HomoL1Data2=Array(Traj3PO1[Dummy[2]])

# HomoPlotL1=Scene(resolution=(1000,1000))
# HomoPlotL1=Makie.lines!(HomoPlotL1,HomoL1Data1[1,:],HomoL1Data1[2,:],color=col_p,linewidth=4,linestyle=:dashdot)
# HomoPlotL1=Makie.lines!(HomoPlotL1,HomoL1Data2[1,:],HomoL1Data2[2,:],color=col_o,linewidth=4,linestyle=:dashdot)
# HomoPlotL1=Makie.lines!(HomoPlotL1,sin.(collect(0:0.01:2*pi)),cos.(collect(0:0.01:2*pi)),linestyle=:dash,linewidth=2)
HomoPlotL1=Plots.plot!(HomoL1Data1[1,:],HomoL1Data1[2,:],color=col_p,linewidth=3*ls,label="",linestyle=:dashdot)
HomoPlotL1=Plots.plot!(HomoL1Data2[1,:],HomoL1Data2[2,:],color=col_o,linewidth=3*ls,label="",linestyle=:dashdot)
HomoPlotL1=Plots.scatter!(HomoPlotL1,[par[1]],[0],label="",
    marker = (:circle, mks1, 0.99, :red, stroke(1.0, 0.9, :black, :solid)))
HomoPlotL1=Plots.plot!(sin.(collect(0:0.01:2*pi)),cos.(collect(0:0.01:2*pi)),
    color=:black,linewidth=2*ls,linestyle=:dash,label="")
savefig(HomoPlotL1,"Figures/TBP_HomoclinicL1.pdf")
display(HomoPlotL1)


