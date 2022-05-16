"""
# This file is used to plot the heteroclinic Poincare cut of double pendulum
# Coded By: KK
# Last Updated: 05/15/2022
"""
##-----------------------------
## First load necessary files and Pacakages
@time include("Functions/LoadPackages.jl")
@time include("Functions/TubeFunctionDP.jl")
@time include("Functions/PlotSettings.jl")

# Define the simulation parameters
T=3.0;tspan=(0.0,T);dt=0.001;Tsim=50;

# Calculate the energy at different equlibrium points
EngEq1=HamiltonianEDP([0,0,0,0],p_nf_edp)
EngEq2=HamiltonianEDP([0,pi,0,0],p_nf_edp)
EngEq3=HamiltonianEDP([pi,0,0,0],p_nf_edp)
EngEq4=HamiltonianEDP([pi,pi,0,0],p_nf_edp)

println("The energy of down-down is: ",EngEq1)
println("The energy of down-up is: ",EngEq2)
println("The energy of up-down is: ",EngEq3)
println("The energy of up-up is: ",EngEq4)

## Load the desired file
EngTarget=0.2;

@load "Results/EDP_L1PO_Eng_$(EngTarget)_Tube.jld2" p_nf_edp x_inc_L1 PO1_Time PO1_Data EngEq1 EngEq2 EngEq3 EngEq4
@load "Results/EDP_L2PO_Eng_$(EngTarget)_Tube.jld2" p_nf_edp x_inc_L2 PO2_Time PO2_Data EngEq1 EngEq2 EngEq3 EngEq4

## Load the data
Tsim=3
N_Inc=9448;
@time @load "Results/EDP_Eng_$(EngTarget)_Tsim_$(Tsim)_POIncNum_$(N_Inc)_F_P_UN.jld2" TrajVarInc1_pos_unstable_NoFric Inc1_1
@time @load "Results/EDP_Eng_$(EngTarget)_Tsim_$(Tsim)_POIncNum_$(N_Inc)_F_N_UN.jld2" TrajVarInc1_neg_unstable_NoFric Inc1_2
@time @load "Results/EDP_Eng_$(EngTarget)_Tsim_$(Tsim)_POIncNum_$(N_Inc)_B_P_S.jld2" TrajVarInc1_pos_stable_NoFric Inc1_3
@time @load "Results/EDP_Eng_$(EngTarget)_Tsim_$(Tsim)_POIncNum_$(N_Inc)_B_N_S.jld2" TrajVarInc1_neg_stable_NoFric Inc1_4
##
Tsim=2
N_Inc=82;
@time @load "Results/EDP_Hetero_Eng_$(EngTarget)_Tsim_$(Tsim)_POIncNum_$(N_Inc)_F_P_UN_PO2.jld2" TrajVarInc2_pos_unstable_NoFric Inc2_1
@time @load "Results/EDP_Hetero_Eng_$(EngTarget)_Tsim_$(Tsim)_POIncNum_$(N_Inc)_F_N_UN_PO2.jld2" TrajVarInc2_neg_unstable_NoFric Inc2_2
@time @load "Results/EDP_Hetero_Eng_$(EngTarget)_Tsim_$(Tsim)_POIncNum_$(N_Inc)_B_P_S_PO2.jld2" TrajVarInc2_pos_stable_NoFric Inc2_3
@time @load "Results/EDP_Hetero_Eng_$(EngTarget)_Tsim_$(Tsim)_POIncNum_$(N_Inc)_B_N_S_PO2.jld2" TrajVarInc2_neg_stable_NoFric Inc2_4

## Plot the Poincare cut
Pos=1;AxisNum=1

PoincarPlot=Plots.plot(size=(650,400),tickfontsize=20,
    xticks = ([-π,0,π], ["-\\pi", "0", "\\pi"]),xlim=(-π-0.1,π+0.1))

CutData1,PoincarPlot=PlotFinalPointPlots(PoincarPlot,
    TrajVarInc1_pos_unstable_NoFric,col_p,AxisNum,LinePlot=0)

CutData2,PoincarPlot=PlotFinalPointPlots(PoincarPlot,
    TrajVarInc2_neg_stable_NoFric,col_o,AxisNum,LinePlot=0)
#savefig(PoincarPlot,"Figures/EDP_Hetero_PoincareCut_L1PosL2Neg_Eng_$(EngTarget).pdf") 
display(PoincarPlot)

## Plot the zoom in plot
PoincarPlot=Plots.plot!(PoincarPlot,xlims=(-0.5,0),ylims=(1.25,1.75))
display(PoincarPlot)
#savefig(PoincarPlot,"Figures/EDP_Hetero_PoincareCut_L1PosL2Neg_Eng_$(EngTarget)_ZoomIn.pdf") 

## Plot the zoom in plot
PoincarPlot=Plots.plot!(PoincarPlot,xlims=(-0.5,0),ylims=(1.25,1.75))
display(PoincarPlot)
#savefig(PoincarPlot,"Figures/EDP_Hetero_PoincareCut_L1PosL2Neg_Eng_$(EngTarget)_ZoomIn.pdf") 

## Plot the zoom in plot with the special points
SpecialPO=[-0.3036270751949797,-0.3036270751949797,20.104866956769918,-21.443747916415752]


## 
PoincarPlot=Plots.plot(size=(650,400),tickfontsize=20,
    xticks = ([-π,0,π], ["-\\pi", "0", "\\pi"]),xlim=(-π-0.1,π+0.1))

CutData1,PoincarPlot=PlotFinalPointPlots(PoincarPlot,
    TrajVarInc1_neg_stable_NoFric,col_o,AxisNum)

CutData2,PoincarPlot=PlotFinalPointPlots(PoincarPlot,
    TrajVarInc2_pos_unstable_NoFric,col_p,AxisNum)

#display(PoincarPlot)
savefig(PoincarPlot,"Figures/EDP_Hetero_PoincareCut_L1NegL2Pos_Eng_$(EngTarget).pdf") 

##
PoincarPlot=Plots.plot!(PoincarPlot,xlims=(-0.5,0),ylims=(-1.75,-1.25))
#display(PoincarPlot)
savefig(PoincarPlot,"Figures/EDP_Hetero_PoincareCut_L1NegL2Pos_Eng_$(EngTarget)_ZoomIn.pdf") 
