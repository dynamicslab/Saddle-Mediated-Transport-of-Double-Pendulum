"""
# This file is used to plot tubes of double pendulum's L1 homoclinic trajectory
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
AxisNum1=1

## Plot the positive plot
PlotTheta2=0;LinePlot=0;

PoincarPlot=Figure()
ax1=Makie.Axis(PoincarPlot[1,1])
CutData1,PoincarPlot=MakiePlotFinalPoint(ax1,
    TrajVarInc1_pos_stable_NoFric,col_o,AxisNum1,PlotTheta2=PlotTheta2,LinePlot=LinePlot)
CutData2,PoincarPlot=MakiePlotFinalPoint(ax1,
    TrajVarInc1_pos_unstable_NoFric,col_p,AxisNum1,PlotTheta2=PlotTheta2,LinePlot=LinePlot)
#Makie.scale!(PoincarPlot, 1, 0.1)
ax1.aspect=AxisAspect(1/0.5)
current_figure()

## Get the symetric points 1 to 3
Points1_1=CutData1[:,(CutData1[1,:].<0.0).*(CutData1[1,:].>-0.002).*(CutData1[3,:].>9.592).*(CutData1[3,:].<9.594)]
Points2_1=CutData2[:,(CutData2[1,:].<-0.0).*(CutData2[2,:].>-0.001).*(CutData2[3,:].>9.594).*(CutData2[3,:].<9.596)]
#
Points1_2=CutData1[:,(CutData1[1,:].<0.0).*(CutData1[1,:].>-0.0001).*(CutData1[3,:].>-9.22)]
Points2_2=CutData2[:,(CutData2[1,:].<-0.0001).*(CutData2[1,:].>-0.00015).*(CutData2[3,:].>-9.2185).*(CutData2[3,:].<-9.218)]
#
Points1_3=CutData1[:,(CutData1[1,:].<0.001).*(CutData1[1,:].>0).*(CutData1[3,:].<-20.205)]
Points2_3=CutData2[:,(CutData2[1,:].<0.0005).*(CutData2[1,:].>0).*(CutData2[3,:].<-20.208)]
#
Points1=[Points1_1 Points1_2 Points1_3]
Points2=[Points2_1 Points2_2 Points2_3]

#PoincarPlot=Scene()
if PlotTheta2==1
    Makie.scatter!(PoincarPlot,Points1[1,:],Points1[4,:],color="red")
    Makie.scatter!(PoincarPlot,Points2[1,:],Points2[4,:],color="green")
else
    Makie.scatter!(PoincarPlot,Points1[1,:],Points1[3,:],color="red")
    Makie.scatter!(PoincarPlot,Points2[1,:],Points2[3,:],color="green")
end
display(PoincarPlot)
current_figure()

## Get the unsymetric points 4 to 5
Points3=CutData1[:,(CutData1[1,:].<-0.297).*(CutData1[1,:].>-0.299).*(CutData1[3,:].>8.315).*(CutData1[3,:].<8.325)]
Points4=CutData2[:,(CutData2[1,:].<-0.297).*(CutData2[1,:].>-0.299).*(CutData2[3,:].>8.315).*(CutData2[3,:].<8.325)]

#PoincarPlot=Scene()
Makie.scatter!(PoincarPlot,Points3[1,:],Points3[3,:],color="red")
Makie.scatter!(PoincarPlot,Points4[1,:],Points4[3,:],color="green")

## Get the unsymetric points 6
Points5=CutData1[:,(CutData1[3,:].<8.335).*(CutData1[3,:].>8.315).*(CutData1[1,:].>0.2960).*(CutData1[1,:].<0.2975)]
Points6=CutData2[:,(CutData2[3,:].<8.335).*(CutData2[3,:].>8.315).*(CutData2[1,:].>0.2960).*(CutData2[1,:].<0.2975)]

#PoincarPlot=Scene()
Makie.scatter!(PoincarPlot,Points5[1,:],Points5[3,:],color="red")
Makie.scatter!(PoincarPlot,Points6[1,:],Points6[3,:],color="green")

## Get the unsymetric points 7
Points7=CutData1[:,(CutData1[3,:].<-3.125).*(CutData1[3,:].>-3.135).*(CutData1[1,:].>-1.692).*(CutData1[1,:].<-1.691)]
Points8=CutData2[:,(CutData2[3,:].<-3.125).*(CutData2[3,:].>-3.135).*(CutData2[1,:].>-1.692).*(CutData2[1,:].<-1.691)]

#PoincarPlot=Scene()
Makie.scatter!(PoincarPlot,Points7[1,:],Points7[3,:],color="red")
Makie.scatter!(PoincarPlot,Points8[1,:],Points8[3,:],color="green")

## Get the intersection
Points9=CutData1[:,(CutData1[3,:].<-3.127).*(CutData1[4,:].>-3.132).*(CutData1[1,:].>1.6916).*(CutData1[1,:].<1.6918)]
Points10=CutData2[:,(CutData2[3,:].<-3.127).*(CutData2[4,:].>-3.132).*(CutData2[1,:].>1.6916).*(CutData2[1,:].<1.6918)]

#PoincarPlot=Scene()
Makie.scatter!(PoincarPlot,Points9[1,:],Points9[3,:],color="red")
Makie.scatter!(PoincarPlot,Points10[1,:],Points10[3,:],color="green")

## Store all the points
L1HomoStablePoints=hcat(Points1,Points3,Points5,Points7,Points9)
L1HomoUnstablePoints=hcat(Points2,Points4,Points6,Points8,Points10)

## Get the index
L1HomoStableIndex=zeros(7,1)
L1HomoUnstableIndex=zeros(7,1)

for i=1:size(L1HomoStablePoints,2)
    L1HomoStableIndex[i,1]=Int(findmax(sum((L1HomoStablePoints[:,i].==CutData1),dims=1))[2][2])
    L1HomoUnstableIndex[i,1]=Int(findmax(sum((L1HomoUnstablePoints[:,i].==CutData2),dims=1))[2][2])
end

L1HomoStableIndex=Int64.(L1HomoStableIndex)
L1HomoUnstableIndex=Int64.(L1HomoUnstableIndex)

## Now plot the trajectory 
AxisNum1=1;AxisNum2=3
mks1=2.5;lw=3;lw2=3;gap_plot=1
fig_resol=(400,400)

for i=1:7
    # Extract the trajectory
    DataTrajStable=Array(TrajVarInc1_pos_stable_NoFric[L1HomoStableIndex[i,1]])
    DataTrajUnstable=Array(TrajVarInc1_pos_unstable_NoFric[L1HomoUnstableIndex[i,1]])

    # Plot the result
    HomoL1Polt=[]
    if AxisNum1==1 && (AxisNum2==3 || AxisNum2==4)
        HomoL1Polt=Plots.plot(size=fig_resol,tickfontsize=tickfontsize,
            xticks = ([-1.5,-0.5,0.5,1.5], ["-1.5", "-0.5", "0.5", "1.5"]),framestyle=:box)
        HomoL1Polt=Plots.scatter!([0],[0],markersize=2,label="",size=resol1,tickfontsize=tickfontsize,
            marker = (:circle, mks1, 0.99, :blue, stroke(1.0, 0.9, :black, :solid)))
    elseif AxisNum1==2 && (AxisNum2==3 || AxisNum2==4)
        HomoL1Polt=Plots.plot(size=fig_resol,tickfontsize=tickfontsize,
            xticks = ([0,π,2π], ["0","\\pi","2\\pi"]))
        HomoL1Polt=Plots.scatter!([π],[0],markersize=2,label="",size=resol1,tickfontsize=tickfontsize,
            marker = (:circle, mks1, 0.99, :blue, stroke(1.0, 0.9, :black, :solid)),framestyle=:box)
    end
    # 
    HomoL1Polt=Plots.plot!(HomoL1Polt,PO1_Data[AxisNum1,1:gap_plot:end],PO1_Data[AxisNum2,1:gap_plot:end],
        color=:black,linewidth=lw2,label="")

    HomoL1Polt=Plots.plot!(HomoL1Polt,DataTrajStable[AxisNum1,:],DataTrajStable[AxisNum2,:],
        color=col_o,linewidth=lw,label="")
    HomoL1Polt=Plots.plot!(HomoL1Polt,DataTrajUnstable[AxisNum1,:],DataTrajUnstable[AxisNum2,:],
        color=col_p,linewidth=lw,label="")
    
    display(HomoL1Polt)
    savefig("Figures/EDP_L1_Homo_Eng_$(EngTarget)_x_$(AxisNum1)_y_$(AxisNum2)_Index_$(i).pdf")
end

## Generate the initial conditions
DataTrajStableInc=zeros(4,7)
DataTrajUnstableInc=zeros(4,7)

for i=1:7
    # Extract the trajectory
    DataTrajStable=Array(TrajVarInc1_pos_stable_NoFric[L1HomoStableIndex[i,1]])
    DataTrajUnstable=Array(TrajVarInc1_pos_unstable_NoFric[L1HomoUnstableIndex[i,1]])

    DataTrajStableInc[:,i]=DataTrajStable[:,1]
    DataTrajUnstableInc[:,i]=DataTrajUnstable[:,1]

end

## Generate animation of the Homoclinic orbits
## Define the parameters for animation
gap=10;color1=col_p;color2=col_o;color3=:green;

# Define the size of the pendulum arm
mks1=38;mks2=38;ls1=20;ls2=2;Ratio=1;

# Define the frame rate of the animation
framerate=60

# Define the boundary of the pendulum animation
xlims1=-0.25;xlims2=0.25;ylims1=-0.25;ylims2=0.25

## Define the data to use
for i=1:7
    Data_HomoclinicStable=Array(TrajVarInc1_pos_stable_NoFric[L1HomoStableIndex[i,1]])
    Data_HomoclinicUnstable=Array(TrajVarInc1_pos_unstable_NoFric[L1HomoUnstableIndex[i,1]])
    Data_Homoclinic=hcat(Data_HomoclinicUnstable,reverse(Data_HomoclinicStable,dims=2))

    ## Generate the animation for L1 point
    DpPlot=Figure(resolution=(1920,1920))
    ax1=Makie.Axis(DpPlot[1,1])

    # Old animation function
    # DpPlotL1,streamsL1=DP_Animate(DpPlotL1,Data_Homoclinic,p_nf_edp,gap,
    #    color1,color2,color3,cmap,mks1,mks2,ls1,Ratio,framerate,xlims1,xlims2,ylims1,ylims2)
    
    path="Figures/EDP_L1_Homo_Eng_$(EngTarget)_Index_$(i).mp4";
    p_animation=[p_nf_edp[1:2],0.17272,0.2286,p_nf_edp[end]]
    DP_Animate(DpPlot,ax1,Data_Homoclinic,p_nf_edp,gap,color1,color2,color3,cmap,mks1,mks2,ls1,Ratio,path,framerate,xlims1,xlims2,ylims1,ylims2)
    
end

