"""
# This file is used to generate the stable and unstable manifold PO of TBP problem.
# Coded By: KK
# Last Updated: 05/15/2022
"""
##-----------------------------
## First load necessary files and Pacakages
@time include("Functions/LoadPackages.jl")
@time include("Functions/TubeFunctionTBP.jl")

## Now generate the Monodromy matrix
@time include("TBP_GenerateMonoMatrix.jl")

# Define the time parameter to simulate the system
T=10.0;tspan=(0.0,T);dt=0.00001;

# Define the mass ratio of the system
μ=9.537e-4;par=[1-μ;μ];

# Calculate the L1 to L3 langrangian points
EqP1=FindL1(μ);EqP2=FindL2(μ);EqP3=FindL3(μ);

# Calculate the energy level at each Lagrangian points
EngP1=Heng([EqP1,0,0,0],par);EngP2=Heng([EqP2,0,0,0],par);EngP3=Heng([EqP3,0,0,0],par);

# Define the energy level to use
EngTarget=-1.5175

## Define the optimizer parameters
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

## Maximum iteration you would like to use
Niter=2000
# The tolerance you want
tol=1e-15
# Updating step
δ=1e-3
# Define initial delta value
Ax=1e-4

## Start optimizing...
# Solve for L1...
# Define initial guess 
x_bar=[EqP1+Ax,0,0,CalVoptEng(EngTarget,EqP1+Ax,par)]
# Calculate the initial condition
x_inc_L1=OptPOTBP(x_bar,par,EngTarget,tspan,cb,Niter,tol,δ)

# Solve for L2
# Define initial guess 
x_bar=[EqP2+Ax,0,0,CalVoptEng(EngTarget,EqP2+Ax,par)]
# Calculate the initial condition
x_inc_L2=OptPOTBP(x_bar,par,EngTarget,tspan,cb,Niter,tol,δ)

## Now calculate the P.O
# Calculate the PO of L1
TBProb=ODEProblem(TBPODE!,x_inc_L1,tspan,par)
sol=DifferentialEquations.solve(TBProb,Vern9(),
    saveat=0.001,tstops=T,abstol=1e-16,reltol=1e-16,callback=cb_pos)

PO1_Time=Array(sol.t)
PO1_Data=Array(sol)

# Calculate the PO of L2
TBProb=ODEProblem(TBPODE!,x_inc_L2,tspan,par)
sol=DifferentialEquations.solve(TBProb,Vern9(),
    saveat=0.001,tstops=T,abstol=1e-16,reltol=1e-16,callback=cb_pos)

PO2_Time=Array(sol.t)
PO2_Data=Array(sol)

## Finally save the PO for later use...
@save "Results/L1PO_μ_$(μ)_Eng_$(EngTarget)_Tube.jld2" μ par x_inc_L1 PO1_Time PO1_Data EqP1 EqP2 EqP3 EngP1 EngP2 EngP3
@save "Results/L2PO_μ_$(μ)_Eng_$(EngTarget)_Tube.jld2" μ par x_inc_L2 PO2_Time PO2_Data EqP1 EqP2 EqP3 EngP1 EngP2 EngP3




