"""
# This file is used to generate the periodic orbits of CRTBP at different energy level.
# Coded By: KK
# Last Updated: 05/15/2022
"""
##-----------------------------
## First load necessary files and Pacakages
@time include("Functions/LoadPackages.jl")
@time include("Functions/TubeFunctionTBP.jl")

## Calculate the energy level of L1 to L3 points

# Define the time parameter to simulate the system
T=10.0;tspan=(0.0,T);dt=0.00001;

# Define the mass ratio of the system
μ=9.537e-4;par=[1-μ;μ];

# Calculate the L1 to L3 langrangian points
EqP1=FindL1(μ);EqP2=FindL2(μ);EqP3=FindL3(μ);

# Calculate the energy at each langrangian points
EngP1=Heng([EqP1,0,0,0],par);
EngP2=Heng([EqP2,0,0,0],par);
EngP3=Heng([EqP3,0,0,0],par);

# Generate a list of target enegy level
N_Inc=10
EngTargetEnd=-1.519
EngTargetList1=[EngP1+(EngTargetEnd-EngP1)*i/N_Inc for i=1:N_Inc];
EngTargetList2=[EngP2+(EngTargetEnd-EngP2)*i/N_Inc for i=1:N_Inc];

## Calculate the a list of initial conditions
# Define a matrix to store the value
POL1=zeros(4,N_Inc);POL2=zeros(4,N_Inc);

# Define the optimizer parameters
# When the particle hits the x axis, we will let the simulation stop
affect!(integrator) = terminate!(integrator)

#Define a continuous callback function
function condition(u,t,integrator) 
    # When u[2]=0, we will stop the integration
    return u[2]
end

# Now define the callback function
cb = ContinuousCallback(condition,affect!,affect!,save_positions = (true,true))
cb_pos = ContinuousCallback(condition,affect!,nothing,save_positions = (true,true))

# Maximum iteration you would like to use
Niter=2000
# The tolerance you want
tol=8e-14
# Updating step
δ=1e-3
# Define initial delta value
Ax=1e-4

## Start optimizing...
# Solve for L1...
for i=1:N_Inc
    # Define initial guess 
    x_bar=[EqP1+Ax,0,0,CalVoptEng(EngTargetList1[i],EqP1+Ax,par)]
    # Calculate the initial condition
    POL1[:,i]=OptPOTBP(x_bar,par,EngTargetList1[i],tspan,cb,Niter,tol,δ)
end

## Solve for L2
for i=1:N_Inc
    # Define initial guess 
    x_bar=[EqP2+Ax,0,0,CalVoptEng(EngTargetList2[i],EqP2+Ax,par)]
    # Calculate the initial condition
    POL2[:,i]=OptPOTBP(x_bar,par,EngTargetList2[i],tspan,cb,Niter,tol,δ)
end

## Now simulate all the initial conditions, and store the simulate values
# Define the ODE problem
TBProb=ODEProblem(TBPODE!,POL1[:,1],tspan,par)

# Prepare for ensemble problem for L1 point
function prob_func_L1(prob,i,repeat)
    remake(prob,u0=POL1[:,i])
  end

# Prepare for ensemble problem for L2 point
function prob_func_L2(prob,i,repeat)
    remake(prob,u0=POL2[:,i])
end

# Define ensemble problem
ensemble_prob_nf_L1 = EnsembleProblem(TBProb,prob_func=prob_func_L1)
ensemble_prob_nf_L2 = EnsembleProblem(TBProb,prob_func=prob_func_L2)

## Solve for ODE
PoincareSim_nf_L1=[]
PoincareSim_nf_L1 = DifferentialEquations.solve(ensemble_prob_nf_L1,
        Vern9(),EnsembleThreads(),callback=cb_pos,
        saveat=0.001, save_start=true,save_end=true,
        abstol=1e-15,reltol=1e-15,trajectories=N_Inc);

PoincareSim_nf_L2=[]        
PoincareSim_nf_L2 = DifferentialEquations.solve(ensemble_prob_nf_L2,
Vern9(),EnsembleThreads(),callback=cb_pos,
saveat=0.001, save_start=true,save_end=true,
abstol=1e-15,reltol=1e-15,trajectories=N_Inc);

## Save the result
@save "Results/L1PO_$(μ)_Ninc_$(N_Inc).jld2" par POL1 PoincareSim_nf_L1 EqP1 EqP2 EqP3 EngP1 EngP2 EngP3
@save "Results/L2PO_$(μ)_Ninc_$(N_Inc).jld2" par POL2 PoincareSim_nf_L2 EqP1 EqP2 EqP3 EngP1 EngP2 EngP3



