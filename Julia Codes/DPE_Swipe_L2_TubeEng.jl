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
@load "Results/EDP_L2PO_Family_Ninc_$(N_Inc).jld2"  PoincareSim_nf_L2 POL2 EngEq1 EngEq2 EngEq3 EngEq4

## Define the callback function to be used
# When the particle hits the x axis, we will let the simulation stop
affect!(integrator) = terminate!(integrator)

#Define a continuous callback function
function conditionSaddle2(u,t,integrator) 
    return u[2]
end

# Now define the callback function
cb4 = ContinuousCallback(conditionSaddle2,affect!,nothing,save_positions = (true,true))

function affect_cut!(integrator) 
    terminate!(integrator)
end

function cut_function3!(u,t,integrator) 
    # Define when to stop when theta1 hits the target
    return u[1]
end

cb_cut3=ContinuousCallback(cut_function3!,affect_cut!,affect_cut!,save_positions = (true,true))

## Compute the monodromy matrix
nx=4

# Calculate the symbolic matrix
MonoMatrixDP=CalSymMonoE(nx,p_nf_edp) 

## Now define the new initial state
PoincarPlot=Scene()
AxisNum=2
StoredCutDataL2=[];EngListL2=zeros(size(POL2,2),1)

for i in [1,4,5,6,8,9,10]
    # Calculate the PO
    Data_Period2,POData2,eigMono2,eigVecMono2,MonoT2=GetValMonoE(POL2[:,i],nx,cb4,tspan,p_nf_edp)

    # Show the eigenvector
    display(eigMono2)

    ## Get the normalized stable and unstable direction
    v_index1=1;v_index2=4;
    vecS1,vecU1=GetSatbleVector(MonoT2,v_index1,v_index2)

    # Determine the offset
    δV=0.001;positive=-1

    # Determine the gap
    gap=1

    # Get initial condition for L1 orbits
    Stable=0
    Inc1_1=GetPeriodMono(Data_Period2[1:4,:],δV,positive,Stable,vecU1,vecS1,gap)

    # Now simulate the L1 orbits
    Tsim=2;everystepsave=false;

    # Simulate the trajectory
    backward=0
    TrajVarInc2_pos_unstable_NoFric=SimulateManifold(Inc1_1,p_nf_edp,backward,Tsim,everystepsave,cb_cut3,experimental=true)

    # Plot the result
    CutDataL2,PoincarPlot=MakiePlotFinalPoint(PoincarPlot,TrajVarInc2_pos_unstable_NoFric,col_p,AxisNum,LinePlot=1)
    Makie.scale!(PoincarPlot,1,0.1)
    display(PoincarPlot)

    # Store the result
    push!(StoredCutDataL2,CutDataL2)

    # Store the energy level
    EngListL2[i,1]=HamiltonianEDP(POL2[:,i],p_nf_edp)
end

## Save the result
#@save "Results/EDP_Swipe_L2_Tube_Eng.jld2" StoredCutDataL2 EngListL2



