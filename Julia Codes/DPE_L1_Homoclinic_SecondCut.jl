"""
# This file is used to plot the L1 homoclinic second cut (Using experimental parameters)
# Coded By: KK
# Last Updated: 05/15/2022
"""
##-----------------------------
## First load necessary files and Pacakages
@time include("Functions/LoadPackages.jl")
@time include("Functions/TubeFunctionDP.jl")
@time include("Functions/PlotSettings.jl")

# Define the simulation parameters
T=5.0;tspan=(0.0,T);dt=0.001;N_Inc=10;

## Load the result: Load the different P.O
@load "Results/EDP_L1PO_Family_Ninc_$(N_Inc).jld2"  PoincareSim_nf_L1 POL1 EngEq1 EngEq2 EngEq3 EngEq4

## Define the callback function to be used
# When the particle hits the x axis, we will let the simulation stop
affect!(integrator) = terminate!(integrator)

#Define a continuous callback function
function conditionSaddle1(u,t,integrator) 
    return u[1]
end

function cut_function1!(u,t,integrator) 
    # Define when to stop when theta1 hits the target
    return u[2]
end

function cut_function2!(u,t,integrator) 
    # Define when to stop when theta1 hits the target
    return u[2]-2*pi
end

function cut_function3!(u,t,integrator) 
    # Define when to stop when theta1 hits the target
    return u[2]-4*pi
end

function cut_function4!(u,t,integrator) 
    # Define when to stop when theta1 hits the target
    return u[2]-6*pi
end

# Now define the callback function
cb1 = ContinuousCallback(conditionSaddle1,affect!,nothing,save_positions = (true,true))
cb_cut1=ContinuousCallback(cut_function1!,affect!,affect!,save_positions = (true,true))
cb_cut2=ContinuousCallback(cut_function2!,affect!,affect!,save_positions = (true,true))
cb_cut3=ContinuousCallback(cut_function3!,affect!,affect!,save_positions = (true,true))
cb_cut4=ContinuousCallback(cut_function2!,nothing,affect!,save_positions = (true,true))
cb_cut5=ContinuousCallback(cut_function4!,affect!,affect!,save_positions = (true,true))
cb_cut6=ContinuousCallback(cut_function3!,nothing,affect!,save_positions = (true,true))


## Compute the monodromy matrix
nx=4

# Calculate the symbolic matrix
MonoMatrixDPE=CalSymMonoE(nx,p_nf_edp) 

## Now determine the energy level of PO you would like to use. Here we use the energy level -0.0698
LevelPin=4
println("The energy level we are using is: $(HamiltonianEDP(POL1[:,LevelPin],p_nf_edp))")

## Calculate the PO using the initial condition you chose
Data_Period1,POData1,eigMono1,eigVecMono1,MonoT1=GetValMonoE(POL1[:,LevelPin],nx,cb1,tspan,p_nf_edp)

# Show the eigenvector
display(eigMono1)

## Get the normalized stable and unstable direction
v_index1=1;v_index2=4;
vecS1,vecU1=GetSatbleVector(MonoT1,v_index1,v_index2)

# Determine the offset
δV=0.001;positive=1

# Determine the gap
gap=1

# Get initial condition for L1 orbits
Stable=0
Inc1_1=GetPeriodMono(Data_Period1[1:4,:],δV,positive,Stable,vecU1,vecS1,gap)
Stable=1
Inc1_2=GetPeriodMono(Data_Period1[1:4,:],δV,-positive,Stable,vecU1,vecS1,gap)

## Now simulate the L1 orbits
Tsim=3;everystepsave=false;

# Simulate the trajectory forward
backward=0
TrajVarInc1_pos_unstable_NoFric=SimulateManifold(Inc1_1,p_nf_edp,backward,Tsim,everystepsave,cb_cut2,experimental=true)
# Simulate the trajectory backward
backward=1
TrajVarInc1_neg_stable_NoFric=SimulateManifold(Inc1_2,p_nf_edp,backward,Tsim,everystepsave,cb_cut1,experimental=true)

## Plot the results
PoincarPlot=Scene()
Axis=1
CutDataL1_Unstable,PoincarPlot=MakiePlotFinalPoint(PoincarPlot,TrajVarInc1_pos_unstable_NoFric,col_p,Axis,LinePlot=1)
CutDataL1_Stable,PoincarPlot=MakiePlotFinalPoint(PoincarPlot,TrajVarInc1_neg_stable_NoFric,col_o,Axis,LinePlot=1)
#
Makie.scale!(PoincarPlot,1,0.1)
display(PoincarPlot)

## Plot the tubes in 3D
MonoPlot=Scene()
ThreeD=1;AxisNum1=1;AxisNum2=2;AxisNum3=3
color1=col_p;color2=col_o
Mod=0;Transform=0
mks=1;lw=2;lineplot=1;gap_line=200
MakiePlotMono(MonoPlot,POData1,TrajVarInc1_pos_unstable_NoFric,col_p,ThreeD,AxisNum1,AxisNum2,AxisNum3,Mod,Transform,mks,lineplot,lw,gap=gap_line)
MakiePlotMono(MonoPlot,POData1,TrajVarInc1_neg_stable_NoFric,col_o,ThreeD,AxisNum1,AxisNum2,AxisNum3,Mod,Transform,mks,lineplot,lw,gap=gap_line)
Makie.scale!(MonoPlot,1,1,0.1)
display(MonoPlot)

## Now simulate the PO until it hits the 4π poincare cut. The following simulation will stop the simulation of the trajectories that flows back
#backward=0
#TrajVarInc1_pos_unstable_NoFric_SecondCut=SimulateManifold(Inc1_1,p_nf_edp,backward,Tsim,everystepsave,CallbackSet(cb_cut3,cb_cut4),experimental=true)
#@save "Results/EDP_Eng_Minus0_7__Tsim_$(Tsim)_cb3_cb4.jld2" TrajVarInc1_pos_unstable_NoFric_SecondCut Inc1_1
@time @load "Results/EDP_Eng_Minus0_7__Tsim_$(Tsim)_cb3_cb4.jld2" TrajVarInc1_pos_unstable_NoFric_SecondCut Inc1_1

## Now simulate the PO until it hits the 4π poincare cut. The following simulation will NOT stop the simulation of the trajectories that flows back
#TrajVarInc1_pos_unstable_NoFric_SecondCut_NoStop=SimulateManifold(Inc1_1,p_nf_edp,backward,5*Tsim,everystepsave,cb_cut3,experimental=true)
#@save "Results/EDP_Eng_Minus0_7__Tsim_$(Tsim)_cb3.jld2" TrajVarInc1_pos_unstable_NoFric_SecondCut_NoStop Inc1_1
#@time @load "Results/EDP_Eng_Minus0_7__Tsim_$(Tsim)_cb3.jld2" TrajVarInc1_pos_unstable_NoFric_SecondCut_NoStop Inc1_1

## Now plot preview the tube
MonoPlot=Scene()
gap_line=50
MakiePlotMono(MonoPlot,POData1,TrajVarInc1_pos_unstable_NoFric_SecondCut,col_p,ThreeD,AxisNum1,AxisNum2,AxisNum3,Mod,Transform,mks,lineplot,lw,gap=gap_line)
Makie.scale!(MonoPlot,1,1,0.1)
display(MonoPlot)

## Plot all the trajectories that reach the 4π plane
MonoPlot=Scene()
#
for i=1:50:size(TrajVarInc1_pos_unstable_NoFric_SecondCut,3)
    DummyDataHolder=Array(TrajVarInc1_pos_unstable_NoFric_SecondCut[i])
    if DummyDataHolder[2,end]>4*π-0.01
        Makie.lines!(MonoPlot,DummyDataHolder[AxisNum1,:],DummyDataHolder[AxisNum2,:],DummyDataHolder[AxisNum3,:],
        linewidth=lw,color=col_p)
    end
end
#
Makie.lines!(MonoPlot,POData1[AxisNum1,:],POData1[AxisNum2,:],POData1[AxisNum3,:],linewidth=5*lw)
Makie.scale!(MonoPlot,1,1,0.1)
Makie.xlabel!("theta_1")
Makie.ylabel!("theta_2")
Makie.ylabel!("dtheta_1")
display(MonoPlot)

## Extract the final points
SecondCutFinalPoints=ExtractIncFromEndPoint(TrajVarInc1_pos_unstable_NoFric_SecondCut)
SecondCut4PiPoints=SecondCutFinalPoints[:,SecondCutFinalPoints[2,:].>4*π-0.01]

## Now plot the second cut final points with the stable tube
PoincarPlot=Scene()
Axis=1
CutDataL1_Stable,PoincarPlot=MakiePlotFinalPoint(PoincarPlot,TrajVarInc1_neg_stable_NoFric,col_o,Axis,LinePlot=1)
Makie.scatter!(PoincarPlot,SecondCut4PiPoints[1,:],SecondCut4PiPoints[3,:],color=col_p,markersize=3,strokewidth=0)
#
Makie.scale!(PoincarPlot,1,0.1)
display(PoincarPlot)

## Plot the Poincare cut using Plots.jl and save it
PoincarPlot=Plots.plot(size=(650,400),tickfontsize=20)
# Plot the result
PoincarPlot=Plots.scatter!(PoincarPlot,SecondCut4PiPoints[1,:],SecondCut4PiPoints[3,:],color=col_p,
    label="",markersize=1.5,markerstrokewidth=0,markerstrokealpha=0,framestyle=:box)
PoincarPlot=Plots.scatter!(PoincarPlot,CutDataL1_Stable[1,:],CutDataL1_Stable[3,:],color=col_o,
    label="",markersize=1.5,markerstrokewidth=0,markerstrokealpha=0,framestyle=:box)
display(PoincarPlot)

# Save the result
savefig(PoincarPlot,"Figures/EDP_L1_Tube_SecondCut.pdf")

## Now simulate all the points from 4π to 6π plane
TrajVarInc1_pos_unstable_NoFric_ThirdCut_NoStop=SimulateManifold(SecondCut4PiPoints,p_nf_edp,backward,Tsim,everystepsave,CallbackSet(cb_cut5,cb_cut6),experimental=true)

## Then plot those trajectory
MonoPlot=Scene()
gap_line=5
MakiePlotMono(MonoPlot,POData1,TrajVarInc1_pos_unstable_NoFric_ThirdCut_NoStop,col_p,ThreeD,AxisNum1,AxisNum2,AxisNum3,Mod,Transform,mks,lineplot,lw,gap=gap_line)
Makie.scale!(MonoPlot,1,1,0.1)
display(MonoPlot)

## Extract the final points
ThirdCutFinalPoints=ExtractIncFromEndPoint(TrajVarInc1_pos_unstable_NoFric_ThirdCut_NoStop)
ThirdCut6PiPoints=ThirdCutFinalPoints[:,ThirdCutFinalPoints[2,:].>6*π-0.01]

## Plot the Poincare cut
PoincarPlot=Scene()
Axis=1
CutDataL1_Stable,PoincarPlot=MakiePlotFinalPoint(PoincarPlot,TrajVarInc1_neg_stable_NoFric,col_o,Axis,LinePlot=1)
Makie.scatter!(PoincarPlot,ThirdCut6PiPoints[1,:],ThirdCut6PiPoints[3,:],color=col_p,markersize=3,strokewidth=0)
#
Makie.scale!(PoincarPlot,1,0.1)
display(PoincarPlot)

## Plot the third cut using Plots.jl
PoincarPlot=Plots.plot(size=(650,400),tickfontsize=20)
# Plot the result
PoincarPlot=Plots.scatter!(PoincarPlot,ThirdCut6PiPoints[1,:],ThirdCut6PiPoints[3,:],color=col_p,
    label="",markersize=1.5,markerstrokewidth=0,markerstrokealpha=0,framestyle=:box)
PoincarPlot=Plots.scatter!(PoincarPlot,CutDataL1_Stable[1,:],CutDataL1_Stable[3,:],color=col_o,
    label="",markersize=1.5,markerstrokewidth=0,markerstrokealpha=0,framestyle=:box)
display(PoincarPlot)

# Save the result
savefig(PoincarPlot,"Figures/EDP_L1_Tube_ThirdCut.pdf")
 


