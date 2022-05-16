"""
# This file is used to identify the Heteroclinic second cut points of double pendulum (using experimental parameters)
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

function cut_function7!(u,t,integrator) 
    # Define when to stop when theta2 hits the target
    return u[1]-pi
end

function cut_function8!(u,t,integrator) 
    # Define when to stop when theta2 hits the target
    return u[2]-pi
end

cb_cut1=ContinuousCallback(cut_function1!,affect_cut!,affect_cut!,save_positions = (true,false))
cb_cut2=ContinuousCallback(cut_function1!,affect_cut!,affect_cut!,save_positions = (true,false))
cb_cut3=ContinuousCallback(cut_function3!,affect_cut!,affect_cut!,save_positions = (true,false))
cb_cut4=ContinuousCallback(cut_function4!,affect_cut!,affect_cut!,save_positions = (true,false))
cb_cut5=ContinuousCallback(cut_function5!,affect_cut!,affect_cut!,save_positions = (true,false))
cb_cut6=ContinuousCallback(cut_function6!,affect_cut!,affect_cut!,save_positions = (true,false))
#
cb_cut7=ContinuousCallback(cut_function1!,nothing,affect_cut!,save_positions = (true,false))
cb_cut8=ContinuousCallback(cut_function7!,affect_cut!,affect_cut!,save_positions = (true,false))
cb_cut9=ContinuousCallback(cut_function8!,affect_cut!,affect_cut!,save_positions = (true,false))
cb_cut10=ContinuousCallback(cut_function1!,affect_cut!,nothing,save_positions = (true,false))

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

## We also want to simulate those points until they flow back 
IncSecondCut=ExtractIncFromEndPoint(UnstableSimL1)
UnstableSimL1SecondCut=SimulateManifold(Inc1_2[:,1:gap:end],p_nf_edp,backward,Tsim,everystepsave,CallbackSet(cb_cut3,cb_cut7,cb_cut8),experimental=true)

## 
backward=1;Tsim=3;everystepsave=false;
StableSimL1=SimulateManifold(Inc1_4[:,1:gap:end],p_nf_edp,backward,Tsim,everystepsave,CallbackSet(cb_cut1,cb_cut3),experimental=true)

## Plot the tubes
ThreeD=1;gap_plot=250;lineplot=1
TrajPlot=Scene()
TrajPlot=PlotPeriod(TrajPlot,PO1_Data,p_nf_edp,Transform,AxisNum1,AxisNum2,AxisNum3,ThreeD)
#TrajPlot=MakiePlotMono(TrajPlot,PO1_Data,UnstableSimL1,color1,ThreeD,AxisNum1,AxisNum2,AxisNum3,Mod,Transform,mks,lineplot,lw,gap=gap_plot)
TrajPlot=MakiePlotMono(TrajPlot,PO1_Data,StableSimL1,color2,ThreeD,AxisNum1,AxisNum2,AxisNum3,Mod,Transform,mks,lineplot,lw,gap=gap_plot)
#
TrajPlot=MakiePlotMono(TrajPlot,PO1_Data,UnstableSimL1SecondCut,color1,ThreeD,AxisNum1,AxisNum2,AxisNum3,Mod,Transform,mks,lineplot,lw,gap=gap_plot)
#
Makie.scale!(TrajPlot,1,1,0.1)
display(TrajPlot)

## Now, use the above information to plot the 3D plots showing how L1 unstable trajectories flows back to itself
ThreeD=1;gap_plot=50;lineplot=1
TrajPlot=Scene()
TrajPlot=PlotPeriod(TrajPlot,PO1_Data,p_nf_edp,Transform,AxisNum1,AxisNum2,AxisNum3,ThreeD)
TrajPlot=MakiePlotMono(TrajPlot,PO1_Data,StableSimL1,color2,ThreeD,AxisNum1,AxisNum2,AxisNum3,Mod,Transform,mks,lineplot,lw,gap=gap_plot)
# Run a for loop to determine which trajectories directly went to the cut plane
for i=1:gap_plot:size(UnstableSimL1SecondCut,3)
    # Extract the datas
    DummyData=Array(UnstableSimL1SecondCut[i])
    # Plot them one by one
    if abs(DummyData[1,end]-DummyData[2,end])<0.01
        Makie.lines!(TrajPlot,DummyData[AxisNum1,:],DummyData[AxisNum2,:],DummyData[AxisNum3,:],linewidth=lw,color=color1)
    end
end
#
Makie.scale!(TrajPlot,1,1,0.1)
display(TrajPlot)

## Now plot the Poincare plane of those tubes
DataL1DoubleCut=zeros(4,size(UnstableSimL1SecondCut,3))
pin=1
for i=1:size(UnstableSimL1SecondCut,3)
    DummyData=Array(UnstableSimL1SecondCut[i])
    if abs(DummyData[1,end]-DummyData[2,end])<0.01 
        DataL1DoubleCut[:,pin]=DummyData[:,end]
        pin=pin+1
    end
end
DataL1DoubleCut=DataL1DoubleCut[:,1:pin-1]

StableEndPionts=ExtractIncFromEndPoint(StableSimL1)

## Plot the Poincare cut using Plots.jl and save it
PoincarPlot=Plots.plot(size=(650,400),tickfontsize=20)
# Plot the result
PoincarPlot=Plots.scatter!(PoincarPlot,DataL1DoubleCut[1,:],DataL1DoubleCut[3,:],color=col_p,
    label="",markersize=1.5,markerstrokewidth=0,markerstrokealpha=0,framestyle=:box)
PoincarPlot=Plots.scatter!(PoincarPlot,StableEndPionts[1,:],StableEndPionts[3,:],color=col_o,
    label="",markersize=1.5,markerstrokewidth=0,markerstrokealpha=0,framestyle=:box)
display(PoincarPlot)

# Save the result
savefig(PoincarPlot,"Figures/EDP_L1_Hetero_SecondCut.pdf")

## Now simulate unstable point backward in time (perform continuation) (For L2)
backward=1;Tsim=3
StableSimL2=SimulateManifold(Inc2_4[:,1:gap:end],p_nf_edp,backward,Tsim,everystepsave,CallbackSet(cb_cut2,cb_cut4),experimental=true)

##
backward=0;Tsim=3
UnstableSimL2=SimulateManifold(Inc2_1[:,1:gap:end],p_nf_edp,backward,Tsim,everystepsave,CallbackSet(cb_cut2,cb_cut4),experimental=true)

## Simulate L2 tubes until they flow back to the θ1=θ2 plane
backward=0;Tsim=3
UnstableSimL2SecondCut=SimulateManifold(Inc2_1[:,1:gap:end],p_nf_edp,backward,Tsim,everystepsave,CallbackSet(cb_cut4,cb_cut9,cb_cut10),experimental=true)

## Plot the results in 3D
ThreeD=1;gap_plot=250;lineplot=1
TrajPlot=Scene()
TrajPlot=PlotPeriod(TrajPlot,PO2_Data,p_nf_edp,Transform,AxisNum1,AxisNum2,AxisNum3,ThreeD)
#TrajPlot=MakiePlotMono(TrajPlot,PO1_Data,UnstableSimL1,color1,ThreeD,AxisNum1,AxisNum2,AxisNum3,Mod,Transform,mks,lineplot,lw,gap=gap_plot)
TrajPlot=MakiePlotMono(TrajPlot,PO2_Data,StableSimL2,color2,ThreeD,AxisNum1,AxisNum2,AxisNum3,Mod,Transform,mks,lineplot,lw,gap=gap_plot)
#
TrajPlot=MakiePlotMono(TrajPlot,PO2_Data,UnstableSimL2SecondCut,color1,ThreeD,AxisNum1,AxisNum2,AxisNum3,Mod,Transform,mks,lineplot,lw,gap=gap_plot)
#
Makie.scale!(TrajPlot,1,1,0.1)
display(TrajPlot)

## Now only plot the trajectories that directly hits the θ1=θ2 plane
ThreeD=1;gap_plot=75;lineplot=1
TrajPlot=Scene()
TrajPlot=PlotPeriod(TrajPlot,PO2_Data,p_nf_edp,Transform,AxisNum1,AxisNum2,AxisNum3,ThreeD)
TrajPlot=MakiePlotMono(TrajPlot,PO2_Data,StableSimL2,color2,ThreeD,AxisNum1,AxisNum2,AxisNum3,Mod,Transform,mks,lineplot,lw,gap=gap_plot)
# Run a for loop to determine which trajectories directly went to the cut plane
for i=1:20:size(UnstableSimL2SecondCut,3)
    # Extract the datas
    DummyData=Array(UnstableSimL2SecondCut[i])
    # Plot them one by one
    if abs(DummyData[1,end]-DummyData[2,end])<0.01
        Makie.lines!(TrajPlot,DummyData[AxisNum1,:],DummyData[AxisNum2,:],DummyData[AxisNum3,:],linewidth=lw,color=color1)
    end
end
#
Makie.scale!(TrajPlot,1,1,0.1)
display(TrajPlot)

## Plot the Poincare section using Plots.jl
DataL2DoubleCut=zeros(4,size(UnstableSimL2SecondCut,3))
pin=1
for i=1:size(UnstableSimL2SecondCut,3)
    DummyData=Array(UnstableSimL2SecondCut[i])
    if abs(DummyData[1,end]-DummyData[2,end])<0.01 
        DataL2DoubleCut[:,pin]=DummyData[:,end]
        pin=pin+1
    end
end
DataL2DoubleCut=DataL2DoubleCut[:,1:pin-1]

StableEndPionts=ExtractIncFromEndPoint(StableSimL2)

## Plot the Poincare cut using Plots.jl and save it
PoincarPlot=Plots.plot(size=(650,400),tickfontsize=20)
# Plot the result
PoincarPlot=Plots.scatter!(PoincarPlot,DataL2DoubleCut[1,:],DataL2DoubleCut[3,:],color=col_p,
    label="",markersize=1.5,markerstrokewidth=0,markerstrokealpha=0,framestyle=:box)
PoincarPlot=Plots.scatter!(PoincarPlot,StableEndPionts[1,:],StableEndPionts[3,:],color=col_o,
    label="",markersize=1.5,markerstrokewidth=0,markerstrokealpha=0,framestyle=:box)
display(PoincarPlot)

# Save the result
savefig(PoincarPlot,"Figures/EDP_L2_Hetero_SecondCut.pdf")

## This plot will plot the L1 and L2 stable and unstable tubes together
gap_plot=250;lineplot=1
TrajPlot=Scene()
TrajPlot=PlotPeriod(TrajPlot,PO2_Data,p_nf_edp,Transform,AxisNum1,AxisNum2,AxisNum3,ThreeD)
TrajPlot=MakiePlotMono(TrajPlot,PO2_Data,StableSimL2,color2,ThreeD,AxisNum1,AxisNum2,AxisNum3,Mod,Transform,mks,lineplot,lw,gap=gap_plot)
TrajPlot=MakiePlotMono(TrajPlot,PO2_Data,UnstableSimL2,color1,ThreeD,AxisNum1,AxisNum2,AxisNum3,Mod,Transform,mks,lineplot,lw,gap=gap_plot)
TrajPlot=MakiePlotMono(TrajPlot,PO1_Data,UnstableSimL1,color1,ThreeD,AxisNum1,AxisNum2,AxisNum3,Mod,Transform,mks,lineplot,lw,gap=gap_plot)
TrajPlot=MakiePlotMono(TrajPlot,PO1_Data,StableSimL1,color2,ThreeD,AxisNum1,AxisNum2,AxisNum3,Mod,Transform,mks,lineplot,lw,gap=gap_plot)
Makie.scale!(TrajPlot,1,1,0.1)
display(TrajPlot)






































