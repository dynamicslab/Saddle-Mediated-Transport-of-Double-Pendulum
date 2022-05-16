"""
# This file is used to plot tubes of double pendulum's L2 saddle points 
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

## Load the POData
@load "Results/EDP_L1PO_Eng_$(EngTarget)_Tube.jld2" p_nf_edp x_inc_L1 PO1_Time PO1_Data EngEq1 EngEq2 EngEq3 EngEq4
@load "Results/EDP_L2PO_Eng_$(EngTarget)_Tube.jld2" p_nf_edp x_inc_L2 PO2_Time PO2_Data EngEq1 EngEq2 EngEq3 EngEq4

## Save the calculation result
@time @load "Results/EDP_Eng_$(EngTarget)_Tsim_$(Tsim)_POIncNum_$(SizeInc)_F_P_UN_PO2.jld2" TrajVarInc2_pos_unstable_NoFric Inc2_1
@time @load "Results/EDP_Eng_$(EngTarget)_Tsim_$(Tsim)_POIncNum_$(SizeInc)_F_N_UN_PO2.jld2" TrajVarInc2_neg_unstable_NoFric Inc2_2
@time @load "Results/EDP_Eng_$(EngTarget)_Tsim_$(Tsim)_POIncNum_$(SizeInc)_B_P_S_PO2.jld2" TrajVarInc2_pos_stable_NoFric Inc2_3
@time @load "Results/EDP_Eng_$(EngTarget)_Tsim_$(Tsim)_POIncNum_$(SizeInc)_B_N_S_PO2.jld2" TrajVarInc2_neg_stable_NoFric Inc2_4

## Now plot the 3D tubes
ThreeD=1;Axis1=1;Axis2=2;Axis3=3
color1=col_p;color2=col_o;
Mod=0;Transform=0
mks=1;lineplot=1;lw=1;lw2=3;camview=(2,8);gap_line=250;gap_point=1;

##
TrajPlot=Plots.plot(size=(650,400),tickfontsize=20,
    xticks = ([-π,0,π], ["-\\pi", "0", "\\pi"]),yticks = ([0,π,2*π], ["0", "\\pi", "2\\pi"]))
TrajPlot=PlotsPlotMono(TrajPlot,PO2_Data,TrajVarInc2_pos_unstable_NoFric,color1,
    ThreeD,Axis1,Axis2,Axis3,Mod,Transform,mks,lineplot,lw,lw2,camview,gap_point,gap_line=gap_line)
TrajPlot=PlotsPlotMono(TrajPlot,PO2_Data,TrajVarInc2_neg_stable_NoFric,color2,
    ThreeD,Axis1,Axis2,Axis3,Mod,Transform,mks,lineplot,lw,lw2,camview,gap_point,gap_line=gap_line)
#savefig(TrajPlot,"Figures/EDP_PO2_3D_Positive_Eng_$(EngTarget).pdf") 
# display(TrajPlot)

## Make this plot larger
TrajPlot=Plots.plot(size=(1300,400),tickfontsize=20,
    xticks = ([-π,0,π], ["-\\pi", "0", "\\pi"]),yticks = ([0,π,2*π], ["0", "\\pi", "2\\pi"]))
TrajPlot=PlotsPlotMono(TrajPlot,PO2_Data,TrajVarInc2_pos_unstable_NoFric,color1,
    ThreeD,Axis1,Axis2,Axis3,Mod,Transform,mks,lineplot,lw,lw2,camview,gap_point,gap_line=gap_line)
TrajPlot=PlotsPlotMono(TrajPlot,PO2_Data,TrajVarInc2_neg_stable_NoFric,color2,
    ThreeD,Axis1,Axis2,Axis3,Mod,Transform,mks,lineplot,lw,lw2,camview,gap_point,gap_line=gap_line)
#savefig(TrajPlot,"Figures/EDP_PO2_3D_Positive_Eng_$(EngTarget)_Large.pdf") 
#display(TrajPlot)

## Plot the negative parts
TrajPlot=Plots.plot(size=(650,400),tickfontsize=20,
    xticks = ([-π,0,π], ["-\\pi", "0", "\\pi"]),yticks = ([0,π,2*π], ["0", "\\pi", "2\\pi"]))
TrajPlot=PlotsPlotMono(TrajPlot,PO2_Data,TrajVarInc2_neg_unstable_NoFric,color1,
    ThreeD,Axis1,Axis2,Axis3,Mod,Transform,mks,lineplot,lw,lw2,camview,gap_point,gap_line=gap_line)
TrajPlot=PlotsPlotMono(TrajPlot,PO2_Data,TrajVarInc2_pos_stable_NoFric,color2,
    ThreeD,Axis1,Axis2,Axis3,Mod,Transform,mks,lineplot,lw,lw2,camview,gap_point,gap_line=gap_line)
savefig(TrajPlot,"Figures/EDP_PO2_3D_Negative_Eng_$(EngTarget).pdf") 
#display(TrajPlot)

## Make this plot larger
TrajPlot=Plots.plot(size=(1300,400),tickfontsize=20,
    xticks = ([-π,0,π], ["-\\pi", "0", "\\pi"]),yticks = ([0,π,2*π], ["0", "\\pi", "2\\pi"]))
TrajPlot=PlotsPlotMono(TrajPlot,PO2_Data,TrajVarInc2_neg_unstable_NoFric,color1,
    ThreeD,Axis1,Axis2,Axis3,Mod,Transform,mks,lineplot,lw,lw2,camview,gap_point,gap_line=gap_line)
TrajPlot=PlotsPlotMono(TrajPlot,PO2_Data,TrajVarInc2_pos_stable_NoFric,color2,
    ThreeD,Axis1,Axis2,Axis3,Mod,Transform,mks,lineplot,lw,lw2,camview,gap_point,gap_line=gap_line)
savefig(TrajPlot,"Figures/EDP_PO2_3D_Negative_Eng_$(EngTarget)_Large.pdf") 
#display(TrajPlot)

## Plot the pos and neg together
TrajPlot=Plots.plot(size=(1300,400),tickfontsize=20,
    xticks = ([-π,0,π], ["-\\pi", "0", "\\pi"]),yticks = ([0,π,2*π], ["0", "\\pi", "2\\pi"]))
TrajPlot=PlotsPlotMono(TrajPlot,PO2_Data,TrajVarInc2_pos_unstable_NoFric,color1,
    ThreeD,Axis1,Axis2,Axis3,Mod,Transform,mks,lineplot,lw,lw2,camview,gap_point,gap_line=gap_line)
TrajPlot=PlotsPlotMono(TrajPlot,PO2_Data,TrajVarInc2_pos_stable_NoFric,color2,
    ThreeD,Axis1,Axis2,Axis3,Mod,Transform,mks,lineplot,lw,lw2,camview,gap_point,gap_line=gap_line)
TrajPlot=PlotsPlotMono(TrajPlot,PO2_Data,TrajVarInc2_neg_unstable_NoFric,color1,
    ThreeD,Axis1,Axis2,Axis3,Mod,Transform,mks,lineplot,lw,lw2,camview,gap_point,gap_line=gap_line)
TrajPlot=PlotsPlotMono(TrajPlot,PO2_Data,TrajVarInc2_neg_stable_NoFric,color2,
    ThreeD,Axis1,Axis2,Axis3,Mod,Transform,mks,lineplot,lw,lw2,camview,gap_point,gap_line=gap_line)
savefig(TrajPlot,"Figures/EDP_PO2_3D_Full_Eng_$(EngTarget)_Large.pdf") 
#display(TrajPlot)

## Now plot the 2D version of the plot
ThreeD=0
TrajPlot=Plots.plot(size=(650,400),tickfontsize=20,
    xticks = ([-π,0,π], ["-\\pi", "0", "\\pi"]),yticks = ([0,π,2*π], ["0", "\\pi", "2\\pi"]))
TrajPlot=PlotsPlotMono(TrajPlot,PO2_Data,TrajVarInc2_pos_unstable_NoFric,color1,
    ThreeD,Axis1,Axis2,Axis3,Mod,Transform,mks,lineplot,lw,lw2,camview,gap_point,gap_line=gap_line)
TrajPlot=PlotsPlotMono(TrajPlot,PO2_Data,TrajVarInc2_neg_stable_NoFric,color2,
    ThreeD,Axis1,Axis2,Axis3,Mod,Transform,mks,lineplot,lw,lw2,camview,gap_point,gap_line=gap_line)
#savefig(TrajPlot,"Figures/EDP_PO2_2D_Positive_Eng_$(EngTarget).pdf") 

## 
TrajPlot=Plots.plot(size=(650,400),tickfontsize=20,
    xticks = ([-π,0,π], ["-\\pi", "0", "\\pi"]),yticks = ([0,π,2*π], ["0", "\\pi", "2\\pi"]))
TrajPlot=PlotsPlotMono(TrajPlot,PO2_Data,TrajVarInc2_neg_unstable_NoFric,color1,
    ThreeD,Axis1,Axis2,Axis3,Mod,Transform,mks,lineplot,lw,lw2,camview,gap_point,gap_line=gap_line)
TrajPlot=PlotsPlotMono(TrajPlot,PO2_Data,TrajVarInc2_pos_stable_NoFric,color2,
    ThreeD,Axis1,Axis2,Axis3,Mod,Transform,mks,lineplot,lw,lw2,camview,gap_point,gap_line=gap_line)
savefig(TrajPlot,"Figures/EDP_PO2_2D_Negative_Eng_$(EngTarget).pdf") 


