"""
# This file is used to identify the Heteroclinic points of double pendulum (using experimental parameters)
# Coded By: KK
# Last Updated: 05/15/2022
"""
##-----------------------------
## First load necessary files and Pacakages
@time include("Functions/LoadPackages.jl")
@time include("Functions/TubeFunctionDP.jl")
@time include("Functions/PlotSettings.jl")

# Define the simulation parameters
T=3.0;tspan=(0.0,T);dt=0.001;

# Calculate the energy at different equlibrium points
EngEq1=HamiltonianEDP([0,0,0,0],p_nf_edp)
EngEq2=HamiltonianEDP([0,pi,0,0],p_nf_edp)
EngEq3=HamiltonianEDP([pi,0,0,0],p_nf_edp)
EngEq4=HamiltonianEDP([pi,pi,0,0],p_nf_edp)
#
println("The energy of down-down is: ",EngEq1)
println("The energy of down-up is: ",EngEq2)
println("The energy of up-down is: ",EngEq3)
println("The energy of up-up is: ",EngEq4)

## Load the desired file
EngTarget=0.2;

@load "Results/EDP_L1PO_Eng_$(EngTarget)_Tube.jld2" p_nf_edp x_inc_L1 PO1_Time PO1_Data EngEq1 EngEq2 EngEq3 EngEq4
@load "Results/EDP_L2PO_Eng_$(EngTarget)_Tube.jld2" p_nf_edp x_inc_L2 PO2_Time PO2_Data EngEq1 EngEq2 EngEq3 EngEq4

HamiltonianEDP(x_inc_L1,p_nf_edp)
HamiltonianEDP(x_inc_L2,p_nf_edp)

## Load the data
N_Inc=12651;
@time @load "Results/EDP_Eng_$(EngTarget)_POIncNum_$(N_Inc)_F_P_UN_HeteroInc.jld2"  Inc1_1
@time @load "Results/EDP_Eng_$(EngTarget)_POIncNum_$(N_Inc)_F_N_UN_HeteroInc.jld2"  Inc1_2
@time @load "Results/EDP_Eng_$(EngTarget)_POIncNum_$(N_Inc)_B_P_S_HeteroInc.jld2"  Inc1_3
@time @load "Results/EDP_Eng_$(EngTarget)_POIncNum_$(N_Inc)_B_N_S_HeteroInc.jld2"  Inc1_4
##
N_Inc=9810;
@time @load "Results/EDP_Eng_$(EngTarget)_POIncNum_$(N_Inc)_F_P_UN_PO2_HeteroInc.jld2"  Inc2_1
@time @load "Results/EDP_Eng_$(EngTarget)_POIncNum_$(N_Inc)_F_N_UN_PO2_HeteroInc.jld2"  Inc2_2
@time @load "Results/EDP_Eng_$(EngTarget)_POIncNum_$(N_Inc)_B_P_S_PO2_HeteroInc.jld2"  Inc2_3
@time @load "Results/EDP_Eng_$(EngTarget)_POIncNum_$(N_Inc)_B_N_S_PO2_HeteroInc.jld2"  Inc2_4

## Define the Poincare Cut
function affect_cut!(integrator) 
    terminate!(integrator)
end

function cut_function1!(u,t,integrator) 
    # Define when to stop when theta1 hits the target
    return u[1]-u[2]
end

function cut_function2!(u,t,integrator) 
    # Define when to stop when theta1 hits the target
    return u[1]+u[2]
end

function cut_function3!(u,t,integrator) 
    # Define when to stop when theta1 hits the target
    return u[2]+pi
end

function cut_function4!(u,t,integrator) 
    # Define when to stop when theta2 hits the target
    return u[1]+pi
end

function cut_function5!(u,t,integrator) 
    # Define when to stop when theta1 hits the target
    return u[2]
end

function cut_function6!(u,t,integrator) 
    # Define when to stop when theta2 hits the target
    return u[1]
end

cb_cut1=ContinuousCallback(cut_function1!,affect_cut!,affect_cut!,save_positions = (true,false))
cb_cut2=ContinuousCallback(cut_function1!,affect_cut!,affect_cut!,save_positions = (true,false))
cb_cut3=ContinuousCallback(cut_function3!,affect_cut!,affect_cut!,save_positions = (true,false))
cb_cut4=ContinuousCallback(cut_function4!,affect_cut!,affect_cut!,save_positions = (true,false))
cb_cut5=ContinuousCallback(cut_function5!,affect_cut!,affect_cut!,save_positions = (true,false))
cb_cut6=ContinuousCallback(cut_function6!,affect_cut!,affect_cut!,save_positions = (true,false))

## Define the plotting parameters
# Now plot the 3D tubes
ThreeD=1;AxisNum1=1;AxisNum2=2;AxisNum3=3
color1=col_p;color2=col_o;
Mod=0;Transform=0
mks=1;lineplot=1;lw=1;lw2=3;camview=(100,10);

## Now simulate stable point forward in time (perform continuation)
gap=1;
backward=0;Tsim=3;everystepsave=false;
UnstableSimL1=SimulateManifold(Inc1_2[:,1:gap:end],p_nf_edp,backward,Tsim,everystepsave,CallbackSet(cb_cut1,cb_cut3),experimental=true)

##
backward=1;Tsim=3;everystepsave=false;
StableSimL1=SimulateManifold(Inc1_4[:,1:gap:end],p_nf_edp,backward,Tsim,everystepsave,CallbackSet(cb_cut1,cb_cut3),experimental=true)

##
ThreeD=1;gap_plot=250;lineplot=1
TrajPlot=Scene()
TrajPlot=PlotPeriod(TrajPlot,PO1_Data,p_nf_edp,Transform,AxisNum1,AxisNum2,AxisNum3,ThreeD)
TrajPlot=MakiePlotMono(TrajPlot,PO1_Data,UnstableSimL1,color1,ThreeD,AxisNum1,AxisNum2,AxisNum3,Mod,Transform,mks,lineplot,lw,gap=gap_plot)
TrajPlot=MakiePlotMono(TrajPlot,PO1_Data,StableSimL1,color2,ThreeD,AxisNum1,AxisNum2,AxisNum3,Mod,Transform,mks,lineplot,lw,gap=gap_plot)
Makie.scale!(TrajPlot,1,1,0.1)
display(TrajPlot)

## Now simulate unstable point backward in time (perform continuation)
backward=1;Tsim=3
StableSimL2=SimulateManifold(Inc2_4[:,1:gap:end],p_nf_edp,backward,Tsim,everystepsave,CallbackSet(cb_cut2,cb_cut4),experimental=true)

##
backward=0;Tsim=3
UnstableSimL2=SimulateManifold(Inc2_1[:,1:gap:end],p_nf_edp,backward,Tsim,everystepsave,CallbackSet(cb_cut2,cb_cut4),experimental=true)

##
gap_plot=250;lineplot=1
TrajPlot=Scene()
TrajPlot=PlotPeriod(TrajPlot,PO2_Data,p_nf_edp,Transform,AxisNum1,AxisNum2,AxisNum3,ThreeD)
TrajPlot=MakiePlotMono(TrajPlot,PO2_Data,StableSimL2,color2,ThreeD,AxisNum1,AxisNum2,AxisNum3,Mod,Transform,mks,lineplot,lw,gap=gap_plot)
TrajPlot=MakiePlotMono(TrajPlot,PO2_Data,UnstableSimL2,color1,ThreeD,AxisNum1,AxisNum2,AxisNum3,Mod,Transform,mks,lineplot,lw,gap=gap_plot)
Makie.scale!(TrajPlot,1,1,0.1)
display(TrajPlot)


## Simulate the 2Pi PO of up-down and down-up saddle
TwoPiPO_UpDown=[0,0,20.704284411720227,-20.03102661145768]
TwoPiPO_DownUp=[0,pi,-9.65638354688606,-20.5152473222608]
#
HamiltonianEDP(TwoPiPO_UpDown,p_nf_edp)
HamiltonianEDP(TwoPiPO_DownUp,p_nf_edp)
#
HamiltonianEDP(Inc2_1[:,4],p_nf_edp)
#
Prob_TwoPiPO_UpDown=ODEProblem(DPE_NoFric!,TwoPiPO_UpDown,tspan,p_nf_edp)
Prob_TwoPiPO_DownUp=ODEProblem(DPE_NoFric!,TwoPiPO_DownUp,tspan,p_nf_edp)
#
Sol_TwoPiPO_UpDown=DifferentialEquations.solve(Prob_TwoPiPO_UpDown,Vern9(),save_start=true,save_end=true,
    saveat=dt,abstol=1e-16,reltol=1e-16,callback=cb_cut5)
Sol_TwoPiPO_DownUp=DifferentialEquations.solve(Prob_TwoPiPO_DownUp,Vern9(),save_start=true,save_end=true,
    saveat=dt,abstol=1e-16,reltol=1e-16,callback=cb_cut6)
#
TwoPiPO_UpDown_Data=Array(Sol_TwoPiPO_UpDown)
TwoPiPO_DownUp_Data=Array(Sol_TwoPiPO_DownUp)

## Plot the tubes together, and save the plots
gap_plot=250;lineplot=1
TrajPlot=Scene(resolution=(1920,1680))
TrajPlot=PlotPeriod(TrajPlot,PO1_Data,p_nf_edp,Transform,AxisNum1,AxisNum2,AxisNum3,ThreeD)
TrajPlot=PlotPeriod(TrajPlot,PO2_Data,p_nf_edp,Transform,AxisNum1,AxisNum2,AxisNum3,ThreeD)
TrajPlot=MakiePlotMono(TrajPlot,PO2_Data,StableSimL2,color2,ThreeD,AxisNum1,AxisNum2,AxisNum3,Mod,Transform,mks,lineplot,lw,gap=gap_plot)
#TrajPlot=MakiePlotMono(TrajPlot,PO2_Data,UnstableSimL2,color1,ThreeD,AxisNum1,AxisNum2,AxisNum3,Mod,Transform,mks,lineplot,lw,gap=gap_plot)
TrajPlot=MakiePlotMono(TrajPlot,PO1_Data,UnstableSimL1,color1,ThreeD,AxisNum1,AxisNum2,AxisNum3,Mod,Transform,mks,lineplot,lw,gap=gap_plot)
#TrajPlot=MakiePlotMono(TrajPlot,PO1_Data,StableSimL1,color2,ThreeD,AxisNum1,AxisNum2,AxisNum3,Mod,Transform,mks,lineplot,lw,gap=gap_plot)
Makie.lines!(TrajPlot,TwoPiPO_UpDown_Data[AxisNum1,:],TwoPiPO_UpDown_Data[AxisNum2,:],TwoPiPO_UpDown_Data[AxisNum3,:],linewidth=10,color=:green)
Makie.lines!(TrajPlot,TwoPiPO_DownUp_Data[AxisNum1,:],TwoPiPO_DownUp_Data[AxisNum2,:],TwoPiPO_DownUp_Data[AxisNum3,:],linewidth=10,color=:red)
#
Makie.scale!(TrajPlot,1,1,0.1)
Makie.xlabel!("theta1")
Makie.ylabel!("theta2")
Makie.zlabel!("dtheta2")
#save("Figures/HeteroclinicTubes3D.png",TrajPlot)
current_figure()
display(TrajPlot)

## Generate the 2D plot of it
 ThreeD=0;gap_point=5;gap_line=250
TrajPlot=Plots.plot(size=(650,400),tickfontsize=20,
    xticks = ([-π,0,π], ["-\\pi", "0", "\\pi"]),yticks = ([0,π,2*π], ["0", "\\pi", "2\\pi"]))
TrajPlot=PlotsPlotMono(TrajPlot,PO1_Data,UnstableSimL1,color1,
    ThreeD,AxisNum1,AxisNum2,AxisNum3,Mod,Transform,mks,lineplot,lw,lw2,camview,gap_point,gap_line=gap_line)
TrajPlot=PlotsPlotMono(TrajPlot,PO2_Data,StableSimL2,color2,
    ThreeD,AxisNum1,AxisNum2,AxisNum3,Mod,Transform,mks,lineplot,lw,lw2,camview,gap_point,gap_line=gap_line)
#savefig(TrajPlot,"Figures/EDP_Heteroclinic_L1Pos_L2Neg_Eng_$(EngTarget).pdf") 
display(TrajPlot)

## 
TrajPlot=Plots.plot(size=(650,400),tickfontsize=20,
    xticks = ([-π,0,π], ["-\\pi", "0", "\\pi"]),yticks = ([0,π,2*π], ["0", "\\pi", "2\\pi"]))
TrajPlot=PlotsPlotMono(TrajPlot,PO1_Data,StableSimL1,color2,
    ThreeD,AxisNum1,AxisNum2,AxisNum3,Mod,Transform,mks,lineplot,lw,lw2,camview,gap_point,gap_line=gap_line)
TrajPlot=PlotsPlotMono(TrajPlot,PO2_Data,UnstableSimL2,color1,
    ThreeD,AxisNum1,AxisNum2,AxisNum3,Mod,Transform,mks,lineplot,lw,lw2,camview,gap_point,gap_line=gap_line)
##savefig(TrajPlot,"Figures/EDP_Heteroclinic_L1Neg_L2Pos_Eng_$(EngTarget).pdf") 
display(TrajPlot)

## Plot the Poincare cut
Pos=1;AxisNum=1
PoincarPlot=Plots.plot(size=(650,400),tickfontsize=20,
    xticks = ([-π,0,π], ["-\\pi", "0", "\\pi"]),xlim=(-π-0.1,π+0.1))
CutData1,PoincarPlot=PlotFinalPointPlots(PoincarPlot,
    UnstableSimL1,col_p,AxisNum,LinePlot=0)
CutData2,PoincarPlot=PlotFinalPointPlots(PoincarPlot,
    StableSimL2,col_o,AxisNum,LinePlot=1)
#savefig(PoincarPlot,"Figures/EDP_Hetero_PoincareCut_L1PosL2Neg_Eng_$(EngTarget).pdf") 
display(PoincarPlot)

## Plot the zoom in plot
PoincarPlot=Plots.plot!(PoincarPlot,xlims=(-0.8,0),ylims=(18,22))
#savefig(PoincarPlot,"Figures/EDP_Hetero_PoincareCut_L1PosL2Neg_Eng_$(EngTarget)_ZoomIn.pdf") 

## Plot the zoom in plot
PoincarPlot=Plots.plot!(PoincarPlot,xlims=(-0.675,-0.25),ylims=(18.5,20.5),
    xticks = ([-0.6,-0.5,-0.4,-0.3], ["-0.6", "-0.5", "-0.4", "-0.3"]),xlim=(-0.7,-0.25))
pointN=[-0.3036270751949797,-0.3036270751949797,0,20.104866956769918,-21.443747916415752,0]
PoincarPlot=Plots.scatter!(PoincarPlot,[pointN[1]],[pointN[4]],label="")
savefig(PoincarPlot,"Figures/EDP_Hetero_PoincareCut_L1PosL2Neg_Eng_WithSpecialPoint_$(EngTarget)_ZoomIn.pdf") 

## 
PoincarPlot=Plots.plot(size=(650,400),tickfontsize=20,
    xticks = ([-π,0,π], ["-\\pi", "0", "\\pi"]),xlim=(-π-0.1,π+0.1))
CutData1,PoincarPlot=PlotFinalPointPlots(PoincarPlot,
    StableSimL1,col_o,AxisNum,LinePlot=0)
CutData2,PoincarPlot=PlotFinalPointPlots(PoincarPlot,
    UnstableSimL2,col_p,AxisNum,LinePlot=1)
#savefig(PoincarPlot,"Figures/EDP_Hetero_PoincareCut_L1NegL2Pos_Eng_$(EngTarget).pdf") 
display(PoincarPlot)

##
PoincarPlot=Plots.plot!(PoincarPlot,xlims=(-0.8,0),ylims=(-22,-18))
#savefig(PoincarPlot,"Figures/EDP_Hetero_PoincareCut_L1NegL2Pos_Eng_$(EngTarget)_ZoomIn.pdf") 

## Plot the Poincare section
CutData1=ExtractIncFromEndPoint(UnstableSimL1)
CutData2=ExtractIncFromEndPoint(StableSimL2)

FinalPoint1=CutData1[:,abs.(CutData1[1,:]-CutData1[2,:]).<=1e-12]
FinalPoint2=CutData2[:,abs.(CutData2[1,:]-CutData2[2,:]).<=1e-12]

## Plot it
mks=0.5;lineplot=1;gap_plot=100;
TrajPlot=Scene()
Makie.scatter!(TrajPlot,FinalPoint1[1,:],FinalPoint1[3,:],color=color1,markersize=5*mks)
Makie.scatter!(TrajPlot,FinalPoint2[1,:],FinalPoint2[3,:],color=color2,markersize=5*mks)
#TrajPlot=PlotPeriod(TrajPlot,PO1_Data,p_nf_edp,Transform,AxisNum1,AxisNum2,AxisNum3,ThreeD)
#TrajPlot=PlotPeriod(TrajPlot,PO2_Data,p_nf_edp,Transform,AxisNum1,AxisNum2,AxisNum3,ThreeD)
#TrajPlot=MakiePlotMono(TrajPlot,PO1_Data,UnstableSim,color1,ThreeD,AxisNum1,AxisNum2,AxisNum3,Mod,Transform,mks,lineplot,lw,gap=gap_plot)
#TrajPlot=MakiePlotMono(TrajPlot,PO2_Data,StableSim,color2,ThreeD,AxisNum1,AxisNum2,AxisNum3,Mod,Transform,mks,lineplot,lw,gap=gap_plot)
Makie.scale!(TrajPlot,1,0.1)
display(TrajPlot)

## Now find out the intersection points, those points will be our heterclinic connection point
# We say two points are connected if the L2 norm difference between the two points is less than certain threshold
MinValDistancePoint=zeros(size(FinalPoint1,2))
MinIndexDistancePoint=Array{Any}(undef,(1,size(FinalPoint1,2)))
@time for i=1:size(FinalPoint1,2)
    PointDistance=(FinalPoint2.-FinalPoint1[:,i]).^2
    NormPointDistance=sum(PointDistance,dims=1)
    minval,minindex=findmin(NormPointDistance)
    MinValDistancePoint[i]=minval
    MinIndexDistancePoint[i]=minindex[2]
end

## Now sort the minimum points
MinNum=2

SortedIndex=sortperm(MinValDistancePoint)

MinValDistancePoint[SortedIndex[1:MinNum]]

HeteroPoints1_Index=SortedIndex[1:MinNum]
HeteroPoints1=FinalPoint1[:,HeteroPoints1_Index]

HeteroPoints2_Index=MinIndexDistancePoint[SortedIndex[1:MinNum]]
HeteroPoints2=FinalPoint2[:,HeteroPoints2_Index]

# Plot the heteroclinic point
mks=5;lineplot=1;gap_plot=100;
TrajPlot=Figure()
ax1=Axis(TrajPlot[1,1])
Makie.scatter!(ax1,FinalPoint1[1,:],FinalPoint1[3,:],color=color1,markersize=mks)
Makie.scatter!(ax1,FinalPoint2[1,:],FinalPoint2[3,:],color=color2,markersize=mks)
Makie.scatter!(ax1,HeteroPoints1[1,:],HeteroPoints1[3,:],color="red",markersize=5*mks)
Makie.scatter!(ax1,HeteroPoints2[1,:],HeteroPoints2[3,:],color="green",markersize=5*mks)
#TrajPlot=MakiePlotMono(TrajPlot,PO1_Data,UnstableSim,color1,ThreeD,AxisNum1,AxisNum2,AxisNum3,Mod,Transform,mks,lineplot,lw,gap=gap_plot)
#TrajPlot=MakiePlotMono(TrajPlot,PO2_Data,StableSim,color2,ThreeD,AxisNum1,AxisNum2,AxisNum3,Mod,Transform,mks,lineplot,lw,gap=gap_plot)
#Makie.scale!(TrajPlot,1,0.1)
display(TrajPlot)

## Now, get the trajectories of the stable and unstable tube, once we connect them, we can have heteroclinic orbits!
TrajPlot=Scene()
AxisNum1=1;AxisNum2=2;AxisNum3=3

HeterclinicOrbits=[]

for i=1:MinNum
    # Plot the stable trajectories
    ~,minIndex=findmin(abs.(sum(CutData2.-HeteroPoints2[:,i],dims=1)))
    HeteroStableData=Array(StableSimL2[minIndex[2]])
    Makie.lines!(TrajPlot,HeteroStableData[AxisNum1,:],HeteroStableData[AxisNum2,:],HeteroStableData[AxisNum3,:],color=color2)

    # Plot the unstable trajectories
    ~,minIndex=findmin(abs.(sum(CutData1.-HeteroPoints1[:,i],dims=1)))
    HeteroUnstableData=Array(UnstableSimL1[minIndex[2]])
    Makie.lines!(TrajPlot,HeteroUnstableData[AxisNum1,:],HeteroUnstableData[AxisNum2,:],HeteroUnstableData[AxisNum3,:],color=color1)

    # Get the heteroclinic orbits for generating animation
    push!(HeterclinicOrbits,hcat(HeteroUnstableData,reverse(HeteroStableData,dims=2)))
end

## Save the Heteroclinic points for later use
HeterTrajUnstable1=UnstableSimL1[HeteroPoints1_Index[1]]
HeterTrajUnstable2=UnstableSimL1[HeteroPoints1_Index[2]]
#
HeterTrajStable1=StableSimL2[HeteroPoints2_Index[1]]
HeterTrajStable2=StableSimL2[HeteroPoints2_Index[2]]
#
#@save "Results/EDP_HeteroResult_$(EngTarget).jld2" p_nf_edp x_inc_L1 x_inc_L2 PO1_Time PO1_Data PO2_Time PO2_Data EngEq1 EngEq2 EngEq3 EngEq4 HeteroPoints1_Index HeteroPoints2_Index HeteroPoints1 HeteroPoints2 HeterclinicOrbits HeterTrajUnstable1 HeterTrajUnstable2 HeterTrajStable1 HeterTrajStable2

##
mks=20
ThreeD=1
TrajPlot=PlotPeriod(TrajPlot,PO1_Data,p_nf_edp,Transform,AxisNum1,AxisNum2,AxisNum3,ThreeD)
TrajPlot=PlotPeriod(TrajPlot,PO2_Data,p_nf_edp,Transform,AxisNum1,AxisNum2,AxisNum3,ThreeD)
#Makie.scatter!(TrajPlot,FinalPoint1[AxisNum1,:],FinalPoint1[AxisNum2,:],FinalPoint1[AxisNum3,:],color=color1,markersize=mks)
#Makie.scatter!(TrajPlot,FinalPoint2[AxisNum1,:],FinalPoint2[AxisNum2,:],FinalPoint2[AxisNum3,:],color=color2,markersize=mks)
Makie.scatter!(TrajPlot,HeteroPoints1[AxisNum1,:],HeteroPoints1[AxisNum2,:],HeteroPoints1[AxisNum3,:],color="red",markersize=5*mks)
Makie.scatter!(TrajPlot,HeteroPoints2[AxisNum1,:],HeteroPoints2[AxisNum2,:],HeteroPoints2[AxisNum3,:],color="green",markersize=5*mks)
Makie.xlabel!("")
Makie.ylabel!("")
if AxisNum3==3
    Makie.zlabel!("")
else
    Makie.zlabel!("")
end
Makie.scale!(TrajPlot,1,1,0.1)
#save("Figures/HeteroclinicOrbits3D.png",TrajPlot)
display(TrajPlot)

## Plot the heteroclinic trajectory in 2D plane
TrajPlot=Scene()
AxisNum1=1;AxisNum2=2;

for i=1:MinNum
    # Plot the stable trajectories
    ~,minIndex=findmin(abs.(sum(CutData2.-HeteroPoints2[:,i],dims=1)))
    HeteroStableData=Array(StableSimL2[minIndex[2]])
    Makie.lines!(TrajPlot,HeteroStableData[AxisNum1,:],HeteroStableData[AxisNum2,:],color=color2)

    # Plot the unstable trajectories
    ~,minIndex=findmin(abs.(sum(CutData1.-HeteroPoints1[:,i],dims=1)))
    HeteroUnstableData=Array(UnstableSimL1[minIndex[2]])
    Makie.lines!(TrajPlot,HeteroUnstableData[AxisNum1,:],HeteroUnstableData[AxisNum2,:],color=color1)
end

display(TrajPlot)

## Plot all the heteroclinic orbits and save it
AxisNum1=1;AxisNum2=2;ThreeD=0;lw=2
pgfplotsx()
for i=1:MinNum
    # Plot the stable trajectories
    ~,minIndex=findmin(abs.(sum(CutData2.-HeteroPoints2[:,i],dims=1)))
    HeteroStableData=Array(StableSimL2[minIndex[2]])
    
    # Plot the unstable trajectories
    ~,minIndex=findmin(abs.(sum(CutData1.-HeteroPoints1[:,i],dims=1)))
    HeteroUnstableData=Array(UnstableSimL1[minIndex[2]])


    mks=5;lineplot=1;gap_plot=100;
    TrajPlot=Plots.plot(size=(650,400),tickfontsize=20,
    xticks = ([-π,0,π], ["-\\pi", "0", "\\pi"]),yticks = ([-π,0,π],["-\\pi", "0", "\\pi"]))

    # Plot the periodic orbits
    TrajPlot=Plots.plot!(TrajPlot,PO1_Data[AxisNum1,:],PO1_Data[AxisNum2,:],color="black",label="")
    TrajPlot=Plots.plot!(TrajPlot,PO2_Data[AxisNum1,:],PO2_Data[AxisNum2,:],color="black",label="")

    # Plot the saddle point
    TrajPlot=Plots.scatter!(TrajPlot,[PO1_Data[AxisNum1,1]],[PO1_Data[AxisNum2,1]],label="",
        marker = (:circle, 6, 0.99, :blue, stroke(2, 0.8, :black, :solid)),framestyle=:box)
    TrajPlot=Plots.scatter!(TrajPlot,[PO2_Data[AxisNum1,1]],[PO2_Data[AxisNum2,1]],label="",
        marker = (:circle, 6, 0.99, :blue, stroke(2, 0.8, :black, :solid)),framestyle=:box)
    
    # Now plot the Heteroclinic points
    TrajPlot=Plots.plot!(TrajPlot,HeteroUnstableData[AxisNum1,:],HeteroUnstableData[AxisNum2,:],color=col_p,
        label="",linewidth=lw,framestyle=:box)
    TrajPlot=Plots.plot!(TrajPlot,HeteroStableData[AxisNum1,:],HeteroStableData[AxisNum2,:],color=col_o,
        label="",linewidth=lw,framestyle=:box)

    #savefig(TrajPlot,"Figures/EDP_Heteroclinic_Trajectory_$(i)_x_$(AxisNum1)_y_$(AxisNum2)_Eng_$(EngTarget).pdf") 
end

## Gnerate the animation for Heteroclinic points
# Define the colors of the pendulum
color1=col_p;color2=col_o;color3=:green;

# Define the size of the pendulum arm
mks1=38;mks2=38;ls1=20;ls2=2;Ratio=1;

# Define the frame rate of the animation
framerate=60;gap=10;

# Define the boundary of the pendulum animation
xlims1=-0.25;xlims2=0.25;ylims1=-0.25;ylims2=0.25

for i=1:MinNum
    println(i)
    DpPlot=Figure(resolution=(1920,1920))
    ax1=Axis(DpPlot[1,1])
    path="Figures/EDoublePendulum_HeteroclinicTrajectories_$(i).mp4";

    DP_Animate(DpPlot,ax1,HeterclinicOrbits[i],p_nf_edp,gap,color1,color2,color3,cmap,mks1,mks2,ls1,Ratio,path,framerate,xlims1,xlims2,ylims1,ylims2)    
end




