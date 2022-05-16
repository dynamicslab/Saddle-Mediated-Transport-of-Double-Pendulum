"""
# This file is used to plot tubes of double pendulum's L2 poincare cut (Using experimental parameters)
# Coded By: KK
# Last Updated: 05/15/2022
"""
##-----------------------------
## First load necessary files and Pacakages
@time include("Functions/LoadPackages.jl")
@time include("Functions/TubeFunctionDP.jl")
@time include("Functions/PlotSettings.jl")

## Define the loading parameters
EngTarget=0.2;Tsim=3;SizeInc=8144;

## Save the calculation result
@time @load "Results/EDP_Eng_$(EngTarget)_Tsim_$(Tsim)_POIncNum_$(SizeInc)_F_P_UN_PO2.jld2" TrajVarInc2_pos_unstable_NoFric Inc2_1
@time @load "Results/EDP_Eng_$(EngTarget)_Tsim_$(Tsim)_POIncNum_$(SizeInc)_F_N_UN_PO2.jld2" TrajVarInc2_neg_unstable_NoFric Inc2_2
@time @load "Results/EDP_Eng_$(EngTarget)_Tsim_$(Tsim)_POIncNum_$(SizeInc)_B_P_S_PO2.jld2" TrajVarInc2_pos_stable_NoFric Inc2_3
@time @load "Results/EDP_Eng_$(EngTarget)_Tsim_$(Tsim)_POIncNum_$(SizeInc)_B_N_S_PO2.jld2" TrajVarInc2_neg_stable_NoFric Inc2_4

## Plot the positiveeps plot
Pos=1;AxisNum=2

PoincarPlot=Plots.plot(size=(650,400),tickfontsize=20)
    #xticks = ([0,π,2*π], ["0","\\pi","2\\pi"]),xlim=(-0.1,2*π+0.1))

CutData1,PoincarPlot=PlotFinalPointPlots(PoincarPlot,
    TrajVarInc2_pos_stable_NoFric,col_o,AxisNum)

CutData2,PoincarPlot=PlotFinalPointPlots(PoincarPlot,
    TrajVarInc2_neg_unstable_NoFric,col_p,AxisNum)

savefig(PoincarPlot,"Figures/EDP_PO2_PoincareCut_PosEps_Eng_$(EngTarget).pdf") 

## Plot the negative plot
PoincarPlot=Plots.plot(size=(650,400),tickfontsize=20)
    #xticks = ([0,π,2*π], ["0","\\pi","2\\pi"]),xlim=(-0.1,2*π+0.1))

CutData1,PoincarPlot=PlotFinalPointPlots(PoincarPlot,
    TrajVarInc2_neg_stable_NoFric,col_o,AxisNum)

CutData2,PoincarPlot=PlotFinalPointPlots(PoincarPlot,
    TrajVarInc2_pos_unstable_NoFric,col_p,AxisNum)

savefig(PoincarPlot,"Figures/EDP_PO2_PoincareCut_NegEps_Eng_$(EngTarget).pdf") 


