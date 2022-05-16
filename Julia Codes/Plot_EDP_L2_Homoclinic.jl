"""
# This file is used to plot tubes of double pendulum's L2 homoclinic trajectory
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
@load "Results/EDP_Eng_$(EngTarget)_Tsim_$(Tsim)_POIncNum_$(SizeInc)_F_P_UN_PO2.jld2" TrajVarInc2_pos_unstable_NoFric Inc2_1
@load "Results/EDP_Eng_$(EngTarget)_Tsim_$(Tsim)_POIncNum_$(SizeInc)_F_N_UN_PO2.jld2" TrajVarInc2_neg_unstable_NoFric Inc2_2
@load "Results/EDP_Eng_$(EngTarget)_Tsim_$(Tsim)_POIncNum_$(SizeInc)_B_P_S_PO2.jld2" TrajVarInc2_pos_stable_NoFric Inc2_3
@load "Results/EDP_Eng_$(EngTarget)_Tsim_$(Tsim)_POIncNum_$(SizeInc)_B_N_S_PO2.jld2" TrajVarInc2_neg_stable_NoFric Inc2_4

## Plot the positiveeps plot
Pos=1;AxisNum=2

PoincarPlot=Figure()
ax1=Makie.Axis(PoincarPlot[1,1])
CutData1,PoincarPlot=MakiePlotFinalPoint(ax1,
    TrajVarInc2_pos_stable_NoFric,col_o,AxisNum,LinePlot=0)
CutData2,PoincarPlot=MakiePlotFinalPoint(ax1,
    TrajVarInc2_neg_unstable_NoFric,col_p,AxisNum,LinePlot=0)
#Makie.scale!(PoincarPlot, 1, 0.1)
ax1.aspect=AxisAspect(1/0.5)
current_figure()

## Get the symetric points 1
Points1=CutData1[:,(CutData1[2,:].<0.0005).*(CutData1[2,:].>0.0).*(CutData1[4,:].>21.052).*(CutData1[4,:].<21.056)]
Points2=CutData2[:,(CutData2[2,:].<0.0005).*(CutData2[2,:].>0.0).*(CutData2[4,:].>21.052).*(CutData2[4,:].<21.056)]

#PoincarPlot=Scene()
Makie.scatter!(PoincarPlot,Points1[2,:],Points1[4,:],color="green")
Makie.scatter!(PoincarPlot,Points2[2,:],Points2[4,:],color="red")

## Get the symetric points 2
Points3=CutData1[:,(CutData1[2,:].<0).*(CutData1[2,:].>-0.0002).*(CutData1[4,:].<18.728).*(CutData1[4,:].>18.724)]
Points4=CutData2[:,(CutData2[2,:].<0).*(CutData2[2,:].>-0.0002).*(CutData2[4,:].<18.728).*(CutData2[4,:].>18.724)]

#PoincarPlot=Scene()
Makie.scatter!(PoincarPlot,Points3[2,:],Points3[4,:],color="green")
Makie.scatter!(PoincarPlot,Points4[2,:],Points4[4,:],color="red")

## Get the symetric points 3
Points5=CutData1[:,(CutData1[2,:].<0.4366).*(CutData1[2,:].>0.4362).*(CutData1[4,:].<17.87).*(CutData1[4,:].>17.865)]
Points6=CutData2[:,(CutData2[2,:].<0.4366).*(CutData2[2,:].>0.4362).*(CutData2[4,:].<17.87).*(CutData2[4,:].>17.865)]

#PoincarPlot=Scene()
Makie.scatter!(PoincarPlot,Points5[2,:],Points5[4,:],color="green")
Makie.scatter!(PoincarPlot,Points6[2,:],Points6[4,:],color="red")

## Get the symetric points 4
Points7=CutData1[:,(CutData1[2,:].<-0.4366).*(CutData1[2,:].>-0.4368).*(CutData1[4,:].<17.866).*(CutData1[4,:].>17.864)]
Points8=CutData2[:,(CutData2[2,:].<-0.4366).*(CutData2[2,:].>-0.4368).*(CutData2[4,:].<17.866).*(CutData2[4,:].>17.864)]

#PoincarPlot=Scene()
Makie.scatter!(PoincarPlot,Points7[2,:],Points7[4,:],color="green")
Makie.scatter!(PoincarPlot,Points8[2,:],Points8[4,:],color="red")

## Store all the homoclinic points
L2HomoStablePoints=hcat(Points1,Points3,Points5,Points7)
L2HomoUnstablePoints=hcat(Points2,Points4,Points6,Points8)

## Get the index
L2HomoStableIndex=zeros(4,1)
L2HomoUnstableIndex=zeros(4,1)

for i=1:size(L2HomoStablePoints,2)
    L2HomoStableIndex[i,1]=Int(findmax(sum((L2HomoStablePoints[:,i].==CutData1),dims=1))[2][2])
    L2HomoUnstableIndex[i,1]=Int(findmax(sum((L2HomoUnstablePoints[:,i].==CutData2),dims=1))[2][2])
end

L2HomoStableIndex=Int64.(L2HomoStableIndex)
L2HomoUnstableIndex=Int64.(L2HomoUnstableIndex)

## Now plot the trajectory 
Axis1=2;Axis2=3
mks1=2.5;lw=2;lw2=2;gap_plot=1
fig_resol=(400,400)

for i=1:4
    # Extract the trajectory
    DataTrajStable=Array(TrajVarInc2_pos_stable_NoFric[L2HomoStableIndex[i,1]])
    DataTrajUnstable=Array(TrajVarInc2_neg_unstable_NoFric[L2HomoUnstableIndex[i,1]])

    # Plot the result
    HomoL2Polt=[]
    if Axis1==2 && Axis2==4
        HomoL2Polt=Plots.plot(size=fig_resol,tickfontsize=tickfontsize,framestyle=:box)#,
            #xticks = ([-1.5,-0.5,0.5,1.5], ["-1.5", "-0.5", "0.5", "1.5"]))
        HomoL2Polt=Plots.scatter!([0],[0],markersize=2,label="",size=resol1,tickfontsize=tickfontsize,
            marker = (:circle, mks1, 0.99, :blue, stroke(1.0, 0.9, :black, :solid)))
    elseif Axis1==2 && Axis2==3
        HomoL2Polt=Plots.plot(size=fig_resol,tickfontsize=tickfontsize,framestyle=:box)#,
            #xticks = ([0,π,2π], ["0","\\pi","2\\pi"]))
        HomoL2Polt=Plots.scatter!([0],[0],markersize=2,label="",size=resol1,tickfontsize=tickfontsize,
            marker = (:circle, mks1, 0.99, :blue, stroke(1.0, 0.9, :black, :solid)))
    end
    # 
    HomoL2Polt=Plots.plot!(HomoL2Polt,PO2_Data[Axis1,1:gap_plot:end],PO2_Data[Axis2,1:gap_plot:end],
        color=:black,linewidth=lw2,label="")

    HomoL2Polt=Plots.plot!(HomoL2Polt,DataTrajStable[Axis1,:],DataTrajStable[Axis2,:],
        color=col_o,linewidth=lw,label="")
    HomoL2Polt=Plots.plot!(HomoL2Polt,DataTrajUnstable[Axis1,:],DataTrajUnstable[Axis2,:],
        color=col_p,linewidth=lw,label="")
    #display(HomoL2Polt)
    savefig("Figures/EDP_L2_Homo_Eng_$(EngTarget)_x_$(Axis1)_y_$(Axis2)_Index_$(i).pdf")
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
for i=1:4
    Data_HomoclinicStable=Array(TrajVarInc2_pos_stable_NoFric[L2HomoStableIndex[i,1]])
    Data_HomoclinicUnstable=Array(TrajVarInc2_neg_unstable_NoFric[L2HomoUnstableIndex[i,1]])
    Data_Homoclinic=hcat(Data_HomoclinicUnstable,reverse(Data_HomoclinicStable,dims=2))

    ## Generate the animation for L1 point
    # Old animation command
    #DpPlotL2=Scene(resolution=(1920,1920))
    #DpPlotL2,streamsL2=DP_Animate(DpPlotL2,Data_Homoclinic,p_nf,gap,
    #    color1,color2,color3,cmap,mks1,mks2,ls1,Ratio,framerate,xlims1,xlims2,ylims1,ylims2)
    #save("Figures/DP_L2_Homo_Eng_$(EngTarget)_Index_$(i).mp4",streamsL2)

    ## Generate the animation for L1 point
    DpPlot=Figure(resolution=(1920,1920))
    ax1=Makie.Axis(DpPlot[1,1])

    path="Figures/EDP_L2_Homo_Eng_$(EngTarget)_Index_$(i).mp4";
    p_animation=[p_nf_edp[1:2],0.17272,0.2286,p_nf_edp[end]]
    DP_Animate(DpPlot,ax1,Data_Homoclinic,p_nf_edp,gap,color1,color2,color3,cmap,mks1,mks2,ls1,Ratio,path,framerate,xlims1,xlims2,ylims1,ylims2)

end

