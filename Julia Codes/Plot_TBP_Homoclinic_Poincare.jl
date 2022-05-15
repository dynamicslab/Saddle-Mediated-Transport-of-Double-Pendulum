"""
# This file is used to plot the L1 and L2 homoclinic poincare cut.
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
#
PDataL2_unstable_pos=GetPoincareData(Traj1PO2)
PDataL2_unstable_neg=GetPoincareData(Traj2PO2)
PDataL2_stable_pos=GetPoincareData(Traj3PO2)
PDataL2_stable_neg=GetPoincareData(Traj4PO2)

## Plot the L1 manifold tubes, that generates the homoclinic orbits and orbits inside the Hill's region
ls=0.5

HomoPoincareL1=Plots.plot(size=(400,400),tickfontsize=tickfontsize,framestyle=:box)
HomoPoincareL1=Plots.plot!(PDataL1_unstable_neg[1,:],PDataL1_unstable_neg[3,:],
    color=col_p,label="")
HomoPoincareL1=Plots.plot!(PDataL1_stable_pos[1,:],PDataL1_stable_pos[3,:],
    color=col_o,label="",xticks=[-0.72,-0.66,-0.6])
# HomoPoincareL1=Plots.plot!(PDataL1_unstable_pos[1,:],PDataL1_unstable_pos[3,:],
#     color=col_p,label="")
# HomoPoincareL1=Plots.plot!(PDataL1_stable_neg[1,:],PDataL1_stable_neg[3,:],
#     color=col_o,label="",xticks=[-0.8,-0.7,-0.6])
HomoPoincareL1=Plots.hline!([0],color=:black,label="",ls=:dash)

savefig(HomoPoincareL1,"Figures/TBP_HomoPoincareL1.pdf")
display(HomoPoincareL1)

## Plot the L2 manifold tubes, that generates the homoclinic orbits and orbits outside the Hill's region
HomoPoincareL2=Plots.plot(size=(400,400),tickfontsize=tickfontsize,framestyle=:box)
HomoPoincareL2=Plots.plot!(PDataL2_unstable_pos[1,:],PDataL2_unstable_pos[3,:],
    color=col_p,label="")
HomoPoincareL2=Plots.plot!(PDataL2_stable_pos[1,:],PDataL2_stable_pos[3,:],
    color=col_o,label="")
# HomoPoincareL2=Plots.plot!(PDataL2_unstable_neg[1,:],PDataL2_unstable_neg[3,:],
#     color=col_p,label="")
# HomoPoincareL2=Plots.plot!(PDataL2_stable_neg[1,:],PDataL2_stable_neg[3,:],
#     color=col_o,label="",xticks=[-2.2,-1.8,-1.4])
HomoPoincareL2=Plots.hline!([0],color=:black,label="",ls=:dash)

savefig(HomoPoincareL2,"Figures/TBP_HomoPoincareL2.pdf")
display(HomoPoincareL2)

## Plot the heteroclinic poincare cut
HeteroPoincare=Plots.plot(size=(400,400),tickfontsize=tickfontsize,framestyle=:box)
HeteroPoincare=Plots.plot!(PDataL1_unstable_pos[2,:],PDataL1_unstable_pos[4,:],
    color=col_p,label="",markersize=0.01)
HeteroPoincare=Plots.plot!(PDataL2_stable_neg[2,:],PDataL2_stable_neg[4,:],
    color=col_o,label="",markersize=0.01)
# HeteroPoincare=Plots.scatter!(PDataL1_unstable_neg[2,:],PDataL1_unstable_neg[4,:],
#     color=col_p,label="",markersize=0.01)
# HeteroPoincare=Plots.scatter!(PDataL2_stable_pos[2,:],PDataL2_stable_pos[4,:],
#     color=col_o,label="",markersize=0.01)
HeteroPoincare=Plots.hline!([0],color=:black,label="",ls=:dash)

savefig(HeteroPoincare,"Figures/TBP_HeteroPoincare.pdf")
display(HeteroPoincare)



