"""
# This file is used to swipe the L1 tube Poincare cut at different energy level (Using experimental parameters)
# Coded By: KK
# Last Updated: 05/15/2022
"""
##-----------------------------
## First load necessary files and Pacakages
@time include("Functions/LoadPackages.jl")
@time include("Functions/TubeFunctionDP.jl")
@time include("Functions/PlotSettings.jl")

# Define the simulation parameters
T=3.0;tspan=(0.0,T);dt=0.001;N_Inc=10;

# Load the result: Load the different P.O
@load "Results/EDP_L1PO_Family_Ninc_$(N_Inc).jld2"  PoincareSim_nf_L1 POL1 EngEq1 EngEq2 EngEq3 EngEq4

## Define the callback function to be used
# When the particle hits the x axis, we will let the simulation stop
affect!(integrator) = terminate!(integrator)

#Define a continuous callback function
function conditionSaddle1(u,t,integrator) 
    return u[1]
end

# Now define the callback function
cb3 = ContinuousCallback(conditionSaddle1,affect!,nothing,save_positions = (true,true))

function affect_cut!(integrator) 
    terminate!(integrator)
end

function cut_function1!(u,t,integrator) 
    # Define when to stop when theta1 hits the target
    return u[2]
end

cb_cut1=ContinuousCallback(cut_function1!,affect_cut!,affect_cut!,save_positions = (true,true))

## Compute the monodromy matrix
nx=4

# Calculate the symbolic matrix
MonoMatrixDP=CalSymMonoE(nx,p_nf_edp) 

## Now define the new initial state
PoincarPlot=Scene()
Axis=1
StoredCutDataL1=[];EngListL1=zeros(size(POL1,2),1)

for i in [1,4,5,6,8,9,10]
    # Calculate the PO
    Data_Period1,POData1,eigMono1,eigVecMono1,MonoT1=GetValMonoE(POL1[:,i],nx,cb3,tspan,p_nf_edp)

    # Show the eigenvector
    display(eigMono1)

    ## Get the normalized stable and unstable direction
    v_index1=1;v_index2=4;
    vecS1,vecU1=GetSatbleVector(MonoT1,v_index1,v_index2)

    # Determine the offset
    δV=0.001;positive=-1

    # Determine the gap
    gap=1

    # Get initial condition for L1 orbits
    Stable=0
    Inc1_1=GetPeriodMono(Data_Period1[1:4,:],δV,positive,Stable,vecU1,vecS1,gap)

    # Now simulate the L1 orbits
    Tsim=3;everystepsave=false;

    # Simulate the trajectory
    backward=0
    TrajVarInc1_pos_unstable_NoFric=SimulateManifold(Inc1_1,p_nf_edp,backward,Tsim,everystepsave,cb_cut1,experimental=true)

    # Plot the result
    CutDataL1,PoincarPlot=MakiePlotFinalPoint(PoincarPlot,TrajVarInc1_pos_unstable_NoFric,col_p,Axis)

    display(PoincarPlot)

    # Store the result
    push!(StoredCutDataL1,CutDataL1)

    # Store the energy level
    EngListL1[i,1]=HamiltonianEDP(POL1[:,i],p_nf_edp)
end

## Save the result
#@save "Results/EDP_Swipe_L1_Tube_Eng.jld2" StoredCutDataL1 EngListL1



