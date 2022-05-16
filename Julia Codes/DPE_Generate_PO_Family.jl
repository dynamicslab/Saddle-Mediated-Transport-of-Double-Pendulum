"""
# This file is used to generate the PO family of the double pendulum (for experimental parameters)
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

## Define the callback function to be used
# When the particle hits the x axis, we will let the simulation stop
affect!(integrator) = terminate!(integrator)

#Define a continuous callback function
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

## Generate a list of target enegy level
N_Inc=10
EngTargetList1=[EngEq2+(EngEq3-EngEq2)*i/N_Inc for i=1:N_Inc];
EngTargetList2=[EngEq3+(0.5-EngEq3)*i/N_Inc for i=1:N_Inc];

##
for i=1:N_Inc
    println(HamiltonianEDP(POL2[:,i],p_nf_edp))
end

## Calculate the a list of initial conditions
# Define a matrix to store the value
POL1=zeros(4,N_Inc);POL2=zeros(4,N_Inc);

## Define the optimization parameters
# Maximum iteration you would like to use
Niter=2000
# The tolerance you want
tol=1e-13
# Define the optimizer
opt=Flux.Optimise.ADAMW()

## Run and solve for the PO initial conditions
# Define which saddle you are tryying to use
index1=1
index2=2

## For L1
dtheta_opt=[3.9795526475085397,2.1940921449990083]

for i=1:N_Inc
    dtheta_opt=POL1[3:4,i]
    POL1[:,i]=DPEOptPO(dtheta_opt,Niter,tol,EngTargetList1[i],p_nf_edp,cb1,dt,T,tspan,index1,opt)
end

## Solve for the second P.O
dtheta_opt=[0.8682571420989703,3.330081525510946]

for i=1:N_Inc
    #if i!=1
    #    dtheta_opt=POL2[3:4,i-1]+[0.5,0.5]
    #end
    dtheta_opt=POL2[3:4,i]
    POL2[:,i]=DPEOptPO(dtheta_opt,Niter,tol,EngTargetList2[i],p_nf_edp,cb2,dt,T,tspan,index2,opt)
end

## Now define the ensembled ODE problem
# Now simulate all the initial conditions, and store the simulate values
# Define the ODE problem
DPProb=ODEProblem(DPE_NoFric!,POL1[:,1],tspan,p_nf_edp)

# Prepare for ensemble problem for L1 point
function prob_func_L1(prob,i,repeat)
    remake(prob,u0=POL1[:,i])
  end

# Prepare for ensemble problem for L2 point
function prob_func_L2(prob,i,repeat)
    remake(prob,u0=POL2[:,i])
end

# Define ensemble problem
ensemble_prob_nf_L1 = EnsembleProblem(DPProb,prob_func=prob_func_L1)
ensemble_prob_nf_L2 = EnsembleProblem(DPProb,prob_func=prob_func_L2)

## Solve for ODE
PoincareSim_nf_L1=[]
PoincareSim_nf_L1 = DifferentialEquations.solve(ensemble_prob_nf_L1,
        Vern9(),EnsembleThreads(),callback=cb3,
        saveat=0.001, save_start=true,save_end=true,
        abstol=1e-15,reltol=1e-15,trajectories=N_Inc);

PoincareSim_nf_L2=[]        
PoincareSim_nf_L2 = DifferentialEquations.solve(ensemble_prob_nf_L2,
        Vern9(),EnsembleThreads(),callback=cb4,
        saveat=0.001, save_start=true,save_end=true,
        abstol=1e-15,reltol=1e-15,trajectories=N_Inc);

## Save the result
@save "Results/EDP_L1PO_Family_Ninc_$(N_Inc).jld2"  PoincareSim_nf_L1 POL1 EngEq1 EngEq2 EngEq3 EngEq4
@save "Results/EDP_L2PO_Family_Ninc_$(N_Inc).jld2"  PoincareSim_nf_L2 POL2 EngEq1 EngEq2 EngEq3 EngEq4



