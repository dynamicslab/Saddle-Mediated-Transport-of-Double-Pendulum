"""
# This file is used to generate the 2Pi PO of the double pendulum at certain energy level
# Coded By: KK
# Last Updated: 05/15/2022
"""
##-----------------------------
## First load necessary files and Pacakages
@time include("Functions/LoadPackages.jl")
@time include("Functions/TubeFunctionDP.jl")
@time include("Functions/PlotSettings.jl")

#set_theme!(theme_black())

# Define the simulation parameters
T=10.0;tspan=(0.0,T);dt=0.001;

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
EngTarget=-0.0698597125144575

## Define the callback function to be used
# When the particle hits the x axis, we will let the simulation stop
affect!(integrator) = terminate!(integrator)

#Define a continuous callback function
function conditionSaddle1(u,t,integrator) 
    return u[2]-3*pi
end

function condition2(u,t,integrator) 
    return u[2]-2*pi
end

# Now define the callback function
cb1 = ContinuousCallback(conditionSaddle1,affect!,nothing,save_positions = (true,true))
cb2 = ContinuousCallback(condition2,affect!,nothing,save_positions = (true,true))

## Define the optimization parameters
# Maximum iteration you would like to use
Niter=300
# The tolerance you want
tol=1e-28
# Define the optimizer
opt=Flux.Optimise.ADAM()

## Run and solve for the PO initial conditions
# Define which saddle you are tryying to use
index1=1
# For energy level: -0.0698597125144575
dtheta_opt=[6.252561555256689,10.74784938633]
# x_Two_Pi_L1=[0,0,-13.015939165310996,21.243558830077944]
# For energy level: -0.0966
#dtheta_opt=[5.893028292853365,9.044750203453109]
# For energy level: -0.16392354904190673
#dtheta_opt=[0.7395119798597876,3.252564148160083]
# For energy level: -0.17503434672230211
# dtheta_opt=[0.4900219680726047,0.517018965104989]
# For energy level: 0.2
# x_Two_Pi_L1=[0,0,-15.361143713745697,27.665633794061108]
# dtheta_opt=[9.65638354688606,20.5152473222608]
EngOpt=HamiltonianEDP([0,pi,dtheta_opt[1],dtheta_opt[2]],p_nf_edp)

## Solve for the first P.O
# (The learning rate needs to be adjusted manully, start with 1e-2 and then decrease the learningrate by 
# two orders of magnitude every 500 epochs)
δ=1e-2
opt=Flux.Optimise.ADAGrad(1e-14)
x_inc_L1=DPEOptPOTwoPi(dtheta_opt,Niter,tol,EngTarget,p_nf_edp,cb1,dt,T,tspan,index1,δ,opt=opt)
dtheta_opt=x_inc_L1[3:4]

## Print the energy of the orbits
println("The identified trajectory energy is: ",HamiltonianEDP(x_inc_L1,p_nf_edp))

## Test run the solution
DP_Prob_Saddle1=ODEProblem(DPE_NoFric!,x_inc_L1,tspan,p_nf_edp)

#
sol_saddle1=DifferentialEquations.solve(DP_Prob_Saddle1,Vern9(),save_start=true,save_end=true,
    saveat=dt,abstol=1e-16,reltol=1e-16,callback=cb1)

#
PO1_Time=Array(sol_saddle1.t)
PO1_Data=Array(sol_saddle1)

## Plot the solution and double check its correct
PeriodPlot=Figure()
ax1=Axis(PeriodPlot[1,1],xlabel="theta1",ylabel="theta2")
Makie.lines!(ax1,PO1_Time,PO1_Data[1,:])
Makie.lines!(ax1,PO1_Time,PO1_Data[2,:])
display(PeriodPlot)

##
PeriodPlot=Figure()
ax1=Axis(PeriodPlot[1,1],xlabel="theta1",ylabel="theta2")
Makie.lines!(ax1,PO1_Data[1,:],PO1_Data[2,:])
display(PeriodPlot)

## Generate the animation of the double pendulum 
# Deinfe animation parameters
gap=2

# Define the colors of the pendulum
color1=col_p;color2=col_o;color3=:green;

# Define the size of the pendulum arm
mks1=38;mks2=38;ls1=20;ls2=2;Ratio=1;

# Define the frame rate of the animation
framerate=30

# Define the boundary of the pendulum animation
xlims1=-0.25;xlims2=0.25;ylims1=-0.25;ylims2=0.25

## Generate the animation for it
DpPlot=Figure(resolution=(1920,1920))
ax1=Axis(DpPlot[1,1])
path="Figures/EDoublePendulum_PO2_DownUp_2Pi_Animate.mp4";
AnimationParameter=[p_nf_edp[1],p_nf_edp[2],0.17272,0.2286,p_nf_edp[end]]
DP_Animate(DpPlot,ax1,PO1_Data,p_nf_edp,gap,color1,color2,color3,cmap,mks1,mks2,ls1,Ratio,path,framerate,xlims1,xlims2,ylims1,ylims2)

## Get the location of the points at θ2=2*pi location
sol_saddle1=DifferentialEquations.solve(DP_Prob_Saddle1,Vern9(),save_start=true,save_end=true,
    saveat=dt,abstol=1e-16,reltol=1e-16,callback=cb2)

PO1_Time=Array(sol_saddle1.t)
PO1_Data=Array(sol_saddle1)

## Plot the solution and double check its correct
PeriodPlot=Figure()
ax1=Axis(PeriodPlot[1,1],xlabel="theta1",ylabel="theta2")
Makie.lines!(ax1,PO1_Time,PO1_Data[1,:])
Makie.lines!(ax1,PO1_Time,PO1_Data[2,:])
display(PeriodPlot)

##
PeriodPlot=Figure()
ax1=Axis(PeriodPlot[1,1],xlabel="theta1",ylabel="theta2")
Makie.lines!(ax1,PO1_Data[1,:],PO1_Data[2,:])
display(PeriodPlot)

## Finally, save the solution for later use
x_Two_Pi_L1=Array(sol_saddle1)[:,end]
x_Two_Pi_L1[1]=0;x_Two_Pi_L1[2]=0;
#@save "Results/EDP_L1_TwoPi_PO_Eng_$(EngTarget)_Tube.jld2" p_nf x_Two_Pi_L1 PO1_Time PO1_Data EngEq1 EngEq2 EngEq3 EngEq4


