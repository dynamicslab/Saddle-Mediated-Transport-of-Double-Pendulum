"""
# This file is used to generate the stable and unstable manifold of the PCR3BP problem.
# Coded By: KK
# Last Updated: 05/15/2022
"""
##-----------------------------
## First load necessary files and Pacakages
@time include("Functions/LoadPackages.jl")
@time include("Functions/TubeFunctionTBP.jl")
@time include("Functions/PlotSettings.jl")

## Now generate the Monodromy matrix
@time include("TBP_GenerateMonoMatrix.jl")

# Define the time parameter to simulate the system
T=10.0;tspan=(0.0,T);dt=0.00001;

# Load the PO data
μ=9.537e-4;EngTarget=-1.5175;

@load "Results/L1PO_μ_$(μ)_Eng_$(EngTarget)_Tube.jld2" μ par x_inc_L1 PO1_Time PO1_Data EqP1 EqP2 EqP3 EngP1 EngP2 EngP3
@load "Results/L2PO_μ_$(μ)_Eng_$(EngTarget)_Tube.jld2" μ par x_inc_L2 PO2_Time PO2_Data EqP1 EqP2 EqP3 EngP1 EngP2 EngP3


## Define the optimizer parameters
# When the particle hits the x axis, we will let the simulation stop
affect!(integrator) = terminate!(integrator)

#Define a continuous callback function
function condition(u,t,integrator) 
    # When u[2]=0, we will stop the integration
    return u[2]
end

# Now define the callback function
cb_pos = ContinuousCallback(condition,affect!,nothing,save_positions = (true,true))

## Define the new augmented initial state
nx=4

X0_PO1=vcat(x_inc_L1,Int.(reshape(Matrix(I,nx,nx),nx*nx,1)[:]))
X0_PO2=vcat(x_inc_L2,Int.(reshape(Matrix(I,nx,nx),nx*nx,1)[:]))

## Solve for the Monodromy matrix
MonoProbPO1=ODEProblem(MonoODE!,X0_PO1,tspan,par)
solPO1=DifferentialEquations.solve(MonoProbPO1,Vern9(),tstops=T,abstol=1e-16,reltol=1e-16,
    callback=cb_pos,saveat = 0.0001,save_start=true,save_end=true)

MonoProbPO2=ODEProblem(MonoODE!,X0_PO2,tspan,par)
solPO2=DifferentialEquations.solve(MonoProbPO2,Vern9(),tstops=T,abstol=1e-16,reltol=1e-16,
    callback=cb_pos,saveat = 0.0001,save_start=true,save_end=true)

TimePO1=Array(solPO1.t)
DataPO1=Array(solPO1)

TimePO2=Array(solPO2.t)
DataPO2=Array(solPO2)

## Gather the info of periodic data
Gap=100

PODataPO1=DataPO1[1:nx,1:Gap:end]
PODataPO2=DataPO2[1:nx,1:Gap:end]

if PODataPO1[:,end]!=DataPO1[1:nx,end]
    PODataPO1=hcat(PODataPO1,DataPO1[1:nx,end]);
end;

if PODataPO2[:,end]!=DataPO2[1:nx,end]
    PODataPO2=hcat(PODataPO2,DataPO2[1:nx,end]);
end;

## After we simulated for one period, let's compute the mono matrix at time T
MonoTPO1=reshape(DataPO1[nx+1:end,end],nx,nx)
MonoTPO2=reshape(DataPO2[nx+1:end,end],nx,nx)

# Then calculate the eigenvalue, we should have λ1=λ2=1, and λ3*λ4=1
println("The eigenvalue of the L1 PO is: ",eigvals(MonoTPO1))

# Then calculate the eigenvalue, we should have λ1=λ2=1, and λ3*λ4=1
println("The eigenvalue of the L2 PO is: ",eigvals(MonoTPO2))

# Then calculate the eigenvector
# Stable direction
println("The eigenvector of L1 PO is: ",real(eigvecs(MonoTPO1)))

# Then calculate the eigenvector
# Stable direction
println("The eigenvector of L2 PO is: ",real(eigvecs(MonoTPO2)))

## Finally, get the stable and unstable direction of L1 and L2 PO
# Stable direction
vecSPO1=real(eigvecs(MonoTPO1)[:,1])
vecSPO1=vecSPO1/norm(vecSPO1)
# Unstable direction
vecUPO1=real(eigvecs(MonoTPO1)[:,4])
vecUPO1=vecUPO1/norm(vecUPO1)

# Stable direction
vecSPO2=real(eigvecs(MonoTPO2)[:,1])
vecSPO2=vecSPO2/norm(vecSPO2)
# Unstable direction
vecUPO2=real(eigvecs(MonoTPO2)[:,4])
vecUPO2=vecUPO1/norm(vecUPO2)

## Now define the initial conditions we should use to simulate the manifold
# How far away we should start from the initial conditions
eps=1e-6

# Unstable manifold simulated forward & backward
XufPO1=PODataPO1.+eps*vecUPO1
XubPO1=PODataPO1.-eps*vecUPO1

# Stable manifold simulated forward & backward
XsfPO1=PODataPO1.+eps*vecSPO1
XsbPO1=PODataPO1.-eps*vecSPO1

# Unstable manifold simulated forward & backward
XufPO2=PODataPO2.+eps*vecUPO2
XubPO2=PODataPO2.-eps*vecUPO2

# Stable manifold simulated forward & backward
XsfPO2=PODataPO2.+eps*vecSPO2
XsbPO2=PODataPO2.-eps*vecSPO2

## Define callback funtion to stop the ODE solver
function condition1(u,t,integrator) 
    if u[1]<0 
        return u[2]
    else
        return false
    end
end

function condition2(u,t,integrator) 
    return u[1]-par[1]
end
    
function condition3(u,t,integrator) 
    if u[1]<0 
        return u[2]
    else
        return false
    end
end

function condition4(u,t,integrator) 
    if abs(u[2])<0.6
        return par[1]-u[1]
    else
        return false
    end
end

# Define a callback function
cb_1=ContinuousCallback(condition1,affect!,affect!,save_positions = (true,true))
cb_2=ContinuousCallback(condition2,affect!,affect!,save_positions = (true,true))
cb_3=ContinuousCallback(condition3,affect!,affect!,save_positions = (true,true))
cb_4=ContinuousCallback(condition4,affect!,affect!,save_positions = (true,true))

## Simulate the initial conditions 
## L1 orbit manifold
Select="L1"
backward=0
# This is the unstable manifold, with postitive ϵ, from mass1 (Sun) to mass2 (Jupiter)
Traj1PO1=SimulateManifoldTBP(XufPO1,par,backward,cb_2);
#Traj1PO1=SimulateManifoldTBP(XufPO1,par,backward,cb_1);
# This is the unstable manifold, with negative ϵ, from mass2 (Jupiter) to mass1 (Sun)
Traj2PO1=SimulateManifoldTBP(XubPO1,par,backward,cb_1);
#Traj2PO1=SimulateManifoldTBP(XubPO1,par,backward,cb_2);
backward=1
# This is the stable manifold, with postitive ϵ, from mass2 (Jupiter) to mass1 (Sun) 
Traj3PO1=SimulateManifoldTBP(XsfPO1,par,backward,cb_1);
#Traj3PO1=SimulateManifoldTBP(XsfPO1,par,backward,cb_2);
# This is the SimulateManifoldTBP manifold, with negative ϵ, from mass1 (Sun) to mass2 (Jupiter)
Traj4PO1=SimulateManifoldTBP(XsbPO1,par,backward,cb_2);
#Traj4PO1=SimulateManifoldTBP(XsbPO1,par,backward,cb_1);

## L2 orbit manifold
Select="L2"
backward=0
# This is the unstable manifold, from mass1 (Sun) to mass2 (Jupiter)
Traj1PO2=SimulateManifoldTBP(XufPO2,par,backward,cb_3);
#Traj1PO2=SimulateManifoldTBP(XufPO2,par,backward,cb_4);
# This is the unstable manifold, from mass2 (Jupiter) to mass1 (Sun)
Traj2PO2=SimulateManifoldTBP(XubPO2,par,backward,cb_4);
#Traj2PO2=SimulateManifoldTBP(XubPO2,par,backward,cb_3);
backward=1
# This is the stable manifold, from mass2 (Jupiter) to mass1 (Sun) 
Traj3PO2=SimulateManifoldTBP(XsfPO2,par,backward,cb_3);
#Traj3PO2=SimulateManifoldTBP(XsfPO2,par,backward,cb_4);
# This is the stable manifold, from mass1 (Sun) to mass2 (Jupiter)
Traj4PO2=SimulateManifoldTBP(XsbPO2,par,backward,cb_4);
#Traj4PO2=SimulateManifoldTBP(XsbPO2,par,backward,cb_3);

## Save all the Results
@save "Results/TBP_Manifold_$(μ)_Gap_$(Gap)_Eng_$(EngTarget).jld2" PODataPO1 PODataPO2 μ par Traj1PO1 Traj2PO1 Traj3PO1 Traj4PO1 Traj1PO2 Traj2PO2 Traj3PO2 Traj4PO2

## Plot the generated L1 PO manifold
ls=0.5 

HamPlotL1=Figure(resolution=(1300,1300))
ax1=Axis(HamPlotL1[1,1])
MakiePlotMonoTBP(ax1,Traj1PO1,col_p,ls)
MakiePlotMonoTBP(ax1,Traj2PO1,col_p,ls)
MakiePlotMonoTBP(ax1,Traj3PO1,col_o,ls)
MakiePlotMonoTBP(ax1,Traj4PO1,col_o,ls)
lines!(ax1,PODataPO1[1,:],PODataPO1[2,:],color=:black,linewidth=2*ls)
lines!(ax1,sin.(collect(0:0.01:2*pi)),cos.(collect(0:0.01:2*pi)),color=:black,linewidth=5*ls,linestyle=:dash)

## Plot the generated L2 PO manifold
HamPlotL2=Figure(resolution=(1300,1300))
ax1=Axis(HamPlotL2[1,1])
MakiePlotMonoTBP(ax1,Traj1PO2,col_p,ls)
MakiePlotMonoTBP(ax1,Traj2PO2,col_p,ls)
MakiePlotMonoTBP(ax1,Traj3PO2,col_o,ls)
MakiePlotMonoTBP(ax1,Traj4PO2,col_o,ls)
lines!(ax1,PODataPO2[1,:],PODataPO2[2,:],color=:black,linewidth=2*ls)
lines!(ax1,sin.(collect(0:0.01:2*pi)),cos.(collect(0:0.01:2*pi)),color=:black,linewidth=5*ls,linestyle=:dash)

## Plot all of them together
HamPlot=Figure(resolution=(1300,1300))
ax1=Axis(HamPlot[1,1])
MakiePlotMonoTBP(ax1,Traj1PO1,col_p,ls)
MakiePlotMonoTBP(ax1,Traj2PO1,col_p,ls)
MakiePlotMonoTBP(ax1,Traj3PO1,col_o,ls)
MakiePlotMonoTBP(ax1,Traj4PO1,col_o,ls)
lines!(ax1,PODataPO1[1,:],PODataPO1[2,:],color=:black,linewidth=2*ls)
MakiePlotMonoTBP(ax1,Traj1PO2,col_p,ls)
MakiePlotMonoTBP(ax1,Traj2PO2,col_p,ls)
MakiePlotMonoTBP(ax1,Traj3PO2,col_o,ls)
MakiePlotMonoTBP(ax1,Traj4PO2,col_o,ls)
lines!(ax1,PODataPO2[1,:],PODataPO2[2,:],color=:black,linewidth=2*ls)
lines!(ax1,sin.(collect(0:0.01:2*pi)),cos.(collect(0:0.01:2*pi)),color=:black,linewidth=5*ls,linestyle=:dash)
#
Makie.xlims!(ax1,-2.5,2.5)
Makie.ylims!(ax1,-2.5,2.5)

HamPlot



