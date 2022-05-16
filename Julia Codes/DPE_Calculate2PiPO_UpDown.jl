"""
# This file is used to generate the 2Pi PO of the double pendulum at certain energy level (Up-Down using experimental parameters)
# Coded By: KK
# Last Updated: 05/15/2022
"""
##-----------------------------
## First load necessary files and Pacakages
@time include("Functions/LoadPackages.jl")
@time include("Functions/TubeFunctionDP.jl")
@time include("Functions/PlotSettings.jl")

# Define the simulation parameters
T=2.0;tspan=(0.0,T);dt=0.001;

# Calculate the energy at different equlibrium points
EngEq1=HamiltonianEDP([0,0,0,0],p_nf_edp)
EngEq2=HamiltonianEDP([0,pi,0,0],p_nf_edp)
EngEq3=HamiltonianEDP([pi,0,0,0],p_nf_edp)
EngEq4=HamiltonianEDP([pi,pi,0,0],p_nf_edp)

println("The energy of down-down is: ",EngEq1)
println("The energy of down-up is: ",EngEq2)
println("The energy of up-down is: ",EngEq3)
println("The energy of up-up is: ",EngEq4)

## Define the energy level you needed
EngTarget=0.1845

## Define the callback function to be used
# When the particle hits the x axis, we will let the simulation stop
affect!(integrator) = terminate!(integrator)

#Define a continuous callback function
function conditionSaddle1(u,t,integrator) 
    return u[2]
end

function condition2(u,t,integrator) 
    return u[2]
end

# Now define the callback function
cb1 = ContinuousCallback(conditionSaddle1,affect!,nothing,save_positions = (true,true))
cb2 = ContinuousCallback(condition2,affect!,affect!,save_positions = (true,true))

## Define the optimization parameters
# Maximum iteration you would like to use
Niter=200
# The tolerance you want
tol=1e-28
# Define the optimizer
opt=[]#Flux.Optimise.ADAGrad(1e-2)

## Run and solve for the PO initial conditions
# Define which saddle you are tryying to use
index1=3
# For energy level Eng=0.18450092638625182
dtheta_opt=[2.208447190476409,2.8468988100236072]
# For energy level Eng=0.19
# dtheta_opt=[2.981532139256569,3.1020481378004927]
# Energy level Eng=0.2
# dtheta_opt=[3.9190069634750704,3.471236285823615]
## Solve for the first P.O
# (The learning rate needs to be adjusted manully, start with 1e-2 and then decrease the learningrate by 
# two orders of magnitude every 500 epochs)
δ=1e-14
x_inc_L2=DPEOptPOTwoPi(dtheta_opt,Niter,tol,EngTarget,p_nf_edp,cb1,dt,T,tspan,index1,δ,opt=opt)
dtheta_opt=x_inc_L2[3:4]

## Print the energy of the orbits
println("The identified trajectory energy is: ",HamiltonianEDP(x_inc_L2,p_nf_edp))

## Test run the solution
x_inc_L2=[pi,0,2.208447190476409,2.8468988100236072]
DP_Prob_Saddle2=ODEProblem(DPE_NoFric!,x_inc_L2,tspan,p_nf_edp)

##
sol_saddle2=DifferentialEquations.solve(DP_Prob_Saddle2,Vern9(),save_start=true,save_end=true,
    saveat=dt,abstol=1e-16,reltol=1e-16,callback=cb1)

##
PO2_Time=Array(sol_saddle2.t)
PO2_Data=Array(sol_saddle2)

## Plot the solution and double check its correct
PeriodPlot=Figure()
ax1=Axis(PeriodPlot[1,1],xlabel="Time",ylabel="theta")
Makie.lines!(ax1,PO2_Time,PO2_Data[1,:])
Makie.lines!(ax1,PO2_Time,PO2_Data[2,:])
display(PeriodPlot)

##
PeriodPlot=Figure()
ax1=Axis(PeriodPlot[1,1],xlabel="theta1",ylabel="theta2")
Makie.lines!(ax1,PO2_Data[1,:],PO2_Data[2,:])
display(PeriodPlot)

## Generate the animation of the double pendulum 
# Deinfe animation parameters
gap=5

# Define the colors of the pendulum
color1=col_p;color2=col_o;color3=:green;

# Define the size of the pendulum arm
mks1=38;mks2=38;ls1=20;ls2=2;Ratio=1;

# Define the frame rate of the animation
framerate=30

# Define the boundary of the pendulum animation
xlims1=-0.25;xlims2=0.25;ylims1=-0.25;ylims2=0.25

## Generate animation
DpPlot=Figure(resolution=(1920,1920))
ax1=Axis(DpPlot[1,1])
path="Figures/EDoublePendulum_PO2_UpDownPeriodic_2Pi_Animate.mp4";
DP_Animate(DpPlot,ax1,PO2_Data,p_nf_edp,gap,color1,color2,color3,cmap,mks1,mks2,ls1,Ratio,path,framerate,xlims1,xlims2,ylims1,ylims2)

## Get the location of the points at θ2=2*pi location
sol_saddle2=DifferentialEquations.solve(DP_Prob_Saddle2,Vern9(),save_start=true,save_end=true,
    saveat=dt,abstol=1e-16,reltol=1e-16,callback=cb2)

PO2_Time=Array(sol_saddle2.t)
PO2_Data=Array(sol_saddle2)

## Plot the solution and double check its correct
PeriodPlot=Figure()
ax1=Axis(PeriodPlot[1,1],xlabel="theta1",ylabel="theta2")
Makie.lines!(ax1,PO2_Time,PO2_Data[1,:])
Makie.lines!(ax1,PO2_Time,PO2_Data[2,:])
display(PeriodPlot)

##
PeriodPlot=Figure()
ax1=Axis(PeriodPlot[1,1],xlabel="theta1",ylabel="theta2")
Makie.lines!(ax1,PO2_Data[1,:],PO2_Data[2,:])
display(PeriodPlot)

## Finally, save the solution for later use
x_Two_Pi_L2=Array(sol_saddle2)[:,end]
x_Two_Pi_L2[1]=0;x_Two_Pi_L2[2]=0;
#@save "Results/DP_L2_TwoPi_PO_UpDown_Eng_$(EngTarget)_Tube.jld2" p_nf x_Two_Pi_L1 PO1_Time PO1_Data EngEq1 EngEq2 EngEq3 EngEq4


