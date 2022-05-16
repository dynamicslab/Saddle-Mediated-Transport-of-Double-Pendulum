"""
# This file is used to generate tubes of double pendulum's saddle points given PO, this calculation is used to generate heteroclinic orbits
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

println("The energy of down-down is: ",EngEq1)
println("The energy of down-up is: ",EngEq2)
println("The energy of up-down is: ",EngEq3)
println("The energy of up-up is: ",EngEq4)

## Load the desired file
EngTarget=0.2

@load "Results/EDP_L1PO_Eng_$(EngTarget)_Tube.jld2" p_nf_edp x_inc_L1 PO1_Time PO1_Data EngEq1 EngEq2 EngEq3 EngEq4
@load "Results/EDP_L2PO_Eng_$(EngTarget)_Tube.jld2" p_nf_edp x_inc_L2 PO2_Time PO2_Data EngEq1 EngEq2 EngEq3 EngEq4

## Define the callback function to be used
# When the particle hits the x axis, we will let the simulation stop
affect!(integrator) = terminate!(integrator)

#Define a continuous callback function
function conditionSaddle2(u,t,integrator) 
    return u[2]
end

# Now define the callback function
cb2 = ContinuousCallback(conditionSaddle2,affect!,affect!,save_positions = (true,true))
cb4 = ContinuousCallback(conditionSaddle2,affect!,nothing,save_positions = (true,true))

## Compute the monodromy matrix
nx=4

# Calculate the symbolic matrix
MonoMatrixDP=CalSymMonoE(nx,p_nf_edp) 

# Now define the new initial state
Data_Period2,POData2,eigMono2,eigVecMono2,MonoT2=GetValMonoE(x_inc_L2,nx,cb4,tspan,p_nf_edp)

# Plot the PO and check wether its correct
PO_Check=Scene()
Makie.lines!(PO_Check,Data_Period2[1,:],Data_Period2[3,:])
Makie.scale!(PO_Check,1,0.1)
display(PO_Check)

## Show the eigenvector
display(eigMono2)

## Get the normalized stable and unstable direction
v_index1=1;v_index2=4;
vecS2,vecU2=GetSatbleVector(MonoT2,v_index1,v_index2)

## Define the cut function
function affect_cut!(integrator) 
    terminate!(integrator)
end

function cut_function5!(u,t,integrator) 
    # Define when to stop when theta1 hits the target
    if u[1]<pi-1
        return u[2]
    else
        return false
    end
end

cb_cut5=ContinuousCallback(cut_function5!,affect_cut!,affect_cut!,save_positions = (true,true))

function cut_function6!(u,t,integrator) 
    # Define when to stop when theta1 hits the target
    if u[1]>2*pi-2
        return u[2]
    else
        return false
    end
end

cb_cut6=ContinuousCallback(cut_function6!,affect_cut!,affect_cut!,save_positions = (true,true));

## Now generate initial condition to simulate the tube
# Determine the offset
δV=0.001;positive=1

# Determine the gap
gap=100

# Get intiial condition for L2 orbits
Stable=0
Inc2_1=GetPeriodMono(Data_Period2[1:4,:],δV,positive,Stable,vecU2,vecS2,gap)
Inc2_2=GetPeriodMono(Data_Period2[1:4,:],δV,-positive,Stable,vecU2,vecS2,gap)
Stable=1
Inc2_3=GetPeriodMono(Data_Period2[1:4,:],δV,positive,Stable,vecU2,vecS2,gap)
Inc2_4=GetPeriodMono(Data_Period2[1:4,:],δV,-positive,Stable,vecU2,vecS2,gap)

## Now simulate the L1 orbits
Tsim=2;everystepsave=false;

## Next simulate the L2 orbits
backward=0
TrajVarInc2_pos_unstable_NoFric=SimulateManifold(Inc2_1,p_nf_edp,backward,Tsim,everystepsave,cb_cut6,experimental=true)
TrajVarInc2_neg_unstable_NoFric=SimulateManifold(Inc2_2,p_nf_edp,backward,Tsim,everystepsave,cb_cut5,experimental=true)
backward=1
TrajVarInc2_pos_stable_NoFric=SimulateManifold(Inc2_3,p_nf_edp,backward,Tsim,everystepsave,cb_cut6,experimental=true)
TrajVarInc2_neg_stable_NoFric=SimulateManifold(Inc2_4,p_nf_edp,backward,Tsim,everystepsave,cb_cut5,experimental=true)

## Test plot the tubes
ThreeD=1;AxisNum1=1;AxisNum2=2;AxisNum3=3
color1=col_p;color2=col_o
Mod=0;Transform=0
mks=1;lw=2;lineplot=1

TrajPlot=Scene()
TrajPlot=PlotPeriod(TrajPlot,POData2,p_nf_edp,Transform,AxisNum1,AxisNum2,AxisNum3,ThreeD)
TrajPlot=MakiePlotMono(TrajPlot,POData2,TrajVarInc2_pos_unstable_NoFric,color1,ThreeD,AxisNum1,AxisNum2,AxisNum3,Mod,Transform,mks,lineplot,lw)
TrajPlot=MakiePlotMono(TrajPlot,POData2,TrajVarInc2_neg_unstable_NoFric,color1,ThreeD,AxisNum1,AxisNum2,AxisNum3,Mod,Transform,mks,lineplot,lw)
TrajPlot=MakiePlotMono(TrajPlot,POData2,TrajVarInc2_pos_stable_NoFric,color2,ThreeD,AxisNum1,AxisNum2,AxisNum3,Mod,Transform,mks,lineplot,lw)
TrajPlot=MakiePlotMono(TrajPlot,POData2,TrajVarInc2_neg_stable_NoFric,color2,ThreeD,AxisNum1,AxisNum2,AxisNum3,Mod,Transform,mks,lineplot,lw)
Makie.scale!(TrajPlot,1,1,0.1)
display(TrajPlot)

## Save the calculation result
@save "Results/EDP_Hetero_Eng_$(EngTarget)_Tsim_$(Tsim)_POIncNum_$(size(Inc2_1,2))_F_P_UN_PO2.jld2" TrajVarInc2_pos_unstable_NoFric Inc2_1
@save "Results/EDP_Hetero_Eng_$(EngTarget)_Tsim_$(Tsim)_POIncNum_$(size(Inc2_2,2))_F_N_UN_PO2.jld2" TrajVarInc2_neg_unstable_NoFric Inc2_2
@save "Results/EDP_Hetero_Eng_$(EngTarget)_Tsim_$(Tsim)_POIncNum_$(size(Inc2_3,2))_B_P_S_PO2.jld2" TrajVarInc2_pos_stable_NoFric Inc2_3
@save "Results/EDP_Hetero_Eng_$(EngTarget)_Tsim_$(Tsim)_POIncNum_$(size(Inc2_4,2))_B_N_S_PO2.jld2" TrajVarInc2_neg_stable_NoFric Inc2_4




