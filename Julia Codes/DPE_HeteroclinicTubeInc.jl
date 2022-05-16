"""
# This file is used to get the initial condition from periodic orbits, this initial condition will then be used to generate the Heteroclinic tubes
# Coded By: KK
# Last Updated: 05/15/2022
"""
##-----------------------------
## First load necessary files and Pacakages
@time include("Functions/LoadPackages.jl")
@time include("Functions/TubeFunctionDP.jl")
@time include("Functions/PlotSettings.jl")

# Define the simulation parameters
T=3.0;tspan=(0.0,T);dt=0.0001;

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

# Define a continuous callback function
function conditionSaddle1(u,t,integrator) 
    return u[1]
end
function conditionSaddle2(u,t,integrator) 
    return u[2]
end

# Now define the callback function
cb1 = ContinuousCallback(conditionSaddle1,affect!,affect!,save_positions = (true,true))
cb2 = ContinuousCallback(conditionSaddle2,affect!,affect!,save_positions = (true,true))
cb3 = ContinuousCallback(conditionSaddle1,affect!,nothing,save_positions = (true,true))
cb4 = ContinuousCallback(conditionSaddle2,affect!,nothing,save_positions = (true,true))

## Compute the monodromy matrix
nx=4

# Calculate the symbolic matrix
MonoMatrixDP=CalSymMonoE(nx,p_nf_edp) 

# Now define the new initial state
Data_Period1,POData1,eigMono1,eigVecMono1,MonoT1=GetValMonoE(x_inc_L1,nx,cb3,tspan,p_nf_edp,dt=0.0001)
Data_Period2,POData2,eigMono2,eigVecMono2,MonoT2=GetValMonoE(x_inc_L2,nx,cb4,tspan,p_nf_edp,dt=0.0001)

## Plot the PO and check wether its correct
PO_Check=Scene()
Makie.lines!(PO_Check,Data_Period1[1,:],Data_Period1[3,:])
Makie.lines!(PO_Check,Data_Period2[1,:],Data_Period2[3,:])
display(PO_Check)

## Show the eigenvector
display(eigMono1)
display(eigMono2)

## Get the normalized stable and unstable direction
v_index1=1;v_index2=4;
vecS1,vecU1=GetSatbleVector(MonoT1,v_index1,v_index2)
vecS2,vecU2=GetSatbleVector(MonoT2,v_index1,v_index2)

## Define the cut function
function affect_cut!(integrator) 
    terminate!(integrator)
end

function cut_function1!(u,t,integrator) 
    # Define when to stop when theta1 hits the target
    return u[2]
end

cb_cut1=ContinuousCallback(cut_function1!,affect_cut!,affect_cut!,save_positions = (true,false))

function cut_function2!(u,t,integrator) 
    # Define when to stop when theta1 hits the target
    return u[2]-2*pi
end

cb_cut2=ContinuousCallback(cut_function2!,affect_cut!,affect_cut!,save_positions = (true,false));

function cut_function3!(u,t,integrator) 
    # Define when to stop when theta1 hits the target
    return u[1]
end

cb_cut3=ContinuousCallback(cut_function3!,affect_cut!,affect_cut!,save_positions = (true,true))

function cut_function4!(u,t,integrator) 
    # Define when to stop when theta1 hits the target
    return u[1]-2*pi
end

cb_cut4=ContinuousCallback(cut_function4!,affect_cut!,affect_cut!,save_positions = (true,true));

## Now generate initial condition to simulate the tube
# Determine the offset
δV=0.001;positive=1

# Determine the gap
gap=1

# Get initial condition for L1 orbits
Stable=0
Inc1_1=GetPeriodMono(Data_Period1[1:4,:],δV,positive,Stable,vecU1,vecS1,gap)
Inc1_2=GetPeriodMono(Data_Period1[1:4,:],δV,-positive,Stable,vecU1,vecS1,gap)
Stable=1
Inc1_3=GetPeriodMono(Data_Period1[1:4,:],δV,positive,Stable,vecU1,vecS1,gap)
Inc1_4=GetPeriodMono(Data_Period1[1:4,:],δV,-positive,Stable,vecU1,vecS1,gap)

## Get intiial condition for L2 orbits
Stable=0
Inc2_1=GetPeriodMono(Data_Period2[1:4,:],δV,positive,Stable,vecU2,vecS2,gap)
Inc2_2=GetPeriodMono(Data_Period2[1:4,:],δV,-positive,Stable,vecU2,vecS2,gap)
Stable=1
Inc2_3=GetPeriodMono(Data_Period2[1:4,:],δV,positive,Stable,vecU2,vecS2,gap)
Inc2_4=GetPeriodMono(Data_Period2[1:4,:],δV,-positive,Stable,vecU2,vecS2,gap)


## Save the calculation result
@save "Results/EDP_Eng_$(EngTarget)_POIncNum_$(size(Inc1_1,2))_F_P_UN_HeteroInc.jld2"  Inc1_1
@save "Results/EDP_Eng_$(EngTarget)_POIncNum_$(size(Inc1_2,2))_F_N_UN_HeteroInc.jld2"  Inc1_2
@save "Results/EDP_Eng_$(EngTarget)_POIncNum_$(size(Inc1_3,2))_B_P_S_HeteroInc.jld2"  Inc1_3
@save "Results/EDP_Eng_$(EngTarget)_POIncNum_$(size(Inc1_4,2))_B_N_S_HeteroInc.jld2"  Inc1_4
## 
@save "Results/EDP_Eng_$(EngTarget)_POIncNum_$(size(Inc2_1,2))_F_P_UN_PO2_HeteroInc.jld2"  Inc2_1
@save "Results/EDP_Eng_$(EngTarget)_POIncNum_$(size(Inc2_2,2))_F_N_UN_PO2_HeteroInc.jld2"  Inc2_2
@save "Results/EDP_Eng_$(EngTarget)_POIncNum_$(size(Inc2_3,2))_B_P_S_PO2_HeteroInc.jld2"  Inc2_3
@save "Results/EDP_Eng_$(EngTarget)_POIncNum_$(size(Inc2_4,2))_B_N_S_PO2_HeteroInc.jld2"  Inc2_4


