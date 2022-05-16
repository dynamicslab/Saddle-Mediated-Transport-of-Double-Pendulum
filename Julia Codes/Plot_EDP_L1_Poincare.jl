"""
# This file is used to plot tubes of double pendulum's L1 poincare cut (using experimental parameter)
# Coded By: KK
# Last Updated: 05/15/2022
"""
##-----------------------------
## First load necessary files and Pacakages
@time include("Functions/LoadPackages.jl")
@time include("Functions/TubeFunctionDP.jl")
@time include("Functions/PlotSettings.jl")

## Define the loading parameters
EngTarget=0.2;Tsim=3;SizeInc=9448;

## Load the POData
@load "Results/EDP_L1PO_Eng_$(EngTarget)_Tube.jld2" p_nf_edp x_inc_L1 PO1_Time PO1_Data EngEq1 EngEq2 EngEq3 EngEq4
@load "Results/EDP_L2PO_Eng_$(EngTarget)_Tube.jld2" p_nf_edp x_inc_L2 PO2_Time PO2_Data EngEq1 EngEq2 EngEq3 EngEq4

## Save the calculation result
@time @load "Results/EDP_Eng_$(EngTarget)_Tsim_$(Tsim)_POIncNum_$(SizeInc)_F_P_UN.jld2" TrajVarInc1_pos_unstable_NoFric Inc1_1
@time @load "Results/EDP_Eng_$(EngTarget)_Tsim_$(Tsim)_POIncNum_$(SizeInc)_F_N_UN.jld2" TrajVarInc1_neg_unstable_NoFric Inc1_2
@time @load "Results/EDP_Eng_$(EngTarget)_Tsim_$(Tsim)_POIncNum_$(SizeInc)_B_P_S.jld2" TrajVarInc1_pos_stable_NoFric Inc1_3
@time @load "Results/EDP_Eng_$(EngTarget)_Tsim_$(Tsim)_POIncNum_$(SizeInc)_B_N_S.jld2" TrajVarInc1_neg_stable_NoFric Inc1_4

## Plot the positiveeps plot
Pos=1;AxisNum=1
PlotTheta2=0;LinePlot=0;

PoincarPlot=Plots.plot(size=(650,400),tickfontsize=20,
    xticks = ([-π,0,π], ["-\\pi", "0", "\\pi"]),xlim=(-π-0.1,π+0.1))

CutData1,PoincarPlot=PlotFinalPointPlots(PoincarPlot,
    TrajVarInc1_pos_stable_NoFric,col_o,AxisNum,PlotTheta2=PlotTheta2,LinePlot=LinePlot)

CutData2,PoincarPlot=PlotFinalPointPlots(PoincarPlot,
    TrajVarInc1_pos_unstable_NoFric,col_p,AxisNum,PlotTheta2=PlotTheta2,LinePlot=LinePlot)
savefig(PoincarPlot,"Figures/EDP_PO1_PoincareCut_PosEps_Eng_$(EngTarget).pdf") 
display(PoincarPlot)

## Plot the negative plot
PoincarPlot=Plots.plot(size=(650,400),tickfontsize=20,
    xticks = ([-π,0,π], ["-\\pi", "0", "\\pi"]),xlim=(-π-0.1,π+0.1))

CutData1,PoincarPlot=PlotFinalPointPlots(PoincarPlot,
    TrajVarInc1_neg_stable_NoFric,col_o,AxisNum,PlotTheta2=PlotTheta2,LinePlot=LinePlot)

CutData2,PoincarPlot=PlotFinalPointPlots(PoincarPlot,
    TrajVarInc1_neg_unstable_NoFric,col_p,AxisNum,PlotTheta2=PlotTheta2,LinePlot=LinePlot)
savefig(PoincarPlot,"Figures/EDP_PO1_PoincareCut_NegEps_Eng_$(EngTarget).pdf") 
#display(PoincarPlot)

## Plot the positiveeps plot (θ1 vs dθ2)
Pos=1;AxisNum=1
PlotTheta2=1;LinePlot=0;

PoincarPlot=Plots.plot(size=(650,400),tickfontsize=20,
    xticks = ([-π,0,π], ["-\\pi", "0", "\\pi"]),xlim=(-π-0.1,π+0.1))

CutData1,PoincarPlot=PlotFinalPointPlots(PoincarPlot,
    TrajVarInc1_pos_stable_NoFric,col_o,AxisNum,PlotTheta2=PlotTheta2,LinePlot=LinePlot)

CutData2,PoincarPlot=PlotFinalPointPlots(PoincarPlot,
    TrajVarInc1_pos_unstable_NoFric,col_p,AxisNum,PlotTheta2=PlotTheta2,LinePlot=LinePlot)

#display(PoincarPlot)
savefig(PoincarPlot,"Figures/EDP_PO1_PoincareCut_PosEps_Eng_$(EngTarget)_1_vs_4.pdf") 

## Plot the negative plot
PoincarPlot=Plots.plot(size=(650,400),tickfontsize=20,
    xticks = ([-π,0,π], ["-\\pi", "0", "\\pi"]),xlim=(-π-0.1,π+0.1))

CutData1,PoincarPlot=PlotFinalPointPlots(PoincarPlot,
    TrajVarInc1_neg_stable_NoFric,col_o,AxisNum,PlotTheta2=PlotTheta2,LinePlot=LinePlot)

CutData2,PoincarPlot=PlotFinalPointPlots(PoincarPlot,
    TrajVarInc1_neg_unstable_NoFric,col_p,AxisNum,PlotTheta2=PlotTheta2,LinePlot=LinePlot)

#display(PoincarPlot)
savefig(PoincarPlot,"Figures/DP_PO1_PoincareCut_NegEps_Eng_$(EngTarget)_1_vs_4.pdf") 

## Plot the positiveeps plot (θ1 vs dθ2) with friction
Pos=1;AxisNum=1
PlotTheta2=1;LinePlot=0;

PoincarPlot=Plots.plot(size=(650,400),tickfontsize=20,
    xticks = ([-π,0,π], ["-\\pi", "0", "\\pi"]),xlim=(-π-0.1,π+0.1))

CutData1,PoincarPlot=PlotFinalPointPlots(PoincarPlot,
    TrajVarInc1_pos_stable_NoFric,col_o,AxisNum,PlotTheta2=PlotTheta2,LinePlot=LinePlot)

CutData1WithFric,PoincarPlot=PlotFinalPointPlots(PoincarPlot,
    TrajVarInc1_pos_stable_WithFric,col_o,AxisNum,PlotTheta2=PlotTheta2,LinePlot=LinePlot)

CutData2,PoincarPlot=PlotFinalPointPlots(PoincarPlot,
    TrajVarInc1_pos_unstable_NoFric,col_p,AxisNum,PlotTheta2=PlotTheta2,LinePlot=LinePlot)
    
CutData2WithFric,PoincarPlot=PlotFinalPointPlots(PoincarPlot,
    TrajVarInc1_pos_unstable_WithFric,col_p,AxisNum,PlotTheta2=PlotTheta2,LinePlot=LinePlot)

display(PoincarPlot)
#savefig(PoincarPlot,"Figures/DP_PO1_PoincareCut_PosEps_Eng_$(EngTarget)_1_vs_4.pdf") 







