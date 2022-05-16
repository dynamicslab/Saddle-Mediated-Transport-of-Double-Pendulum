"""
# This file is used to generate the PO of the double pendulum near saddle points (using experimental parameters)
# Coded By: KK
# Last Updated: 05/15/2022
"""
##-----------------------------
## First load necessary files and Pacakages
@time include("Functions/LoadPackages.jl")
@time include("Functions/TubeFunctionDP.jl")
@time include("Functions/PlotSettings.jl")

set_theme!(theme_black())

# Define the simulation parameters
T=1.8;tspan=(0.0,T);dt=0.001;

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
EngTarget=0.2

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


## Define the optimization parameters
# Maximum iteration you would like to use
Niter=1000
# The tolerance you want
tol=1e-24

## Run and solve for the PO initial conditions
# Define which saddle you are tryying to use
index1=1
index2=2

## Solve for the first P.O
dtheta_opt1=[14.75091676562244,8.84629055058124]
#[14.750916765499058,8.846290550507014]
#[14.750916762728504,8.846290548832497]
#[ 14.750916751349072,8.846290541994946]
#[14.750916702032384,8.846290512401495]
#[14.750917079128582,8.846290739176203]
#[14.750919359634997,8.846292109732534]
#[14.750916307575729,8.846290275035459]
HamiltonianEDP([0.0,pi,dtheta_opt1[1],dtheta_opt1[2]],p_nf_edp)

##
# Define the optimizer
opt=Flux.Optimise.ADAGrad(1e-10)
@time x_inc_L1=DPEOptPO(dtheta_opt1,Niter,tol,EngTarget,p_nf_edp,cb1,dt,T,tspan,index1,opt)
dtheta_opt1=x_inc_L1[3:4]

## Solve for the second P.O
dtheta_opt2=[1.2873939240043606,4.869279573811112]

HamiltonianEDP([pi,0,dtheta_opt2[1],dtheta_opt2[2]],p_nf_edp)

##
opt=Flux.Optimise.ADAM(1e-7)
@time x_inc_L2=DPEOptPO(dtheta_opt2,Niter,tol,EngTarget,p_nf_edp,cb2,dt,T,tspan,index2,opt)
dtheta_opt2=x_inc_L2[3:4]

## Print the energy of the orbits
println("The saddle L1 orbit energy is: ",HamiltonianEDP(x_inc_L1,p_nf_edp))

println("The saddle L2 orbit energy is: ",HamiltonianEDP(x_inc_L2,p_nf_edp))

## Test run the solution
DP_Prob_Saddle1=ODEProblem(DPE_NoFric!,x_inc_L1,tspan,p_nf_edp)
DP_Prob_Saddle2=ODEProblem(DPE_NoFric!,x_inc_L2,tspan,p_nf_edp)

##
sol_saddle1=DifferentialEquations.solve(DP_Prob_Saddle1,Vern9(),save_start=true,save_end=true,
    saveat=0.00001,abstol=1e-16,reltol=1e-16,callback=cb3)

##
sol_saddle2=DifferentialEquations.solve(DP_Prob_Saddle2,Vern9(),save_start=true,save_end=true,
    saveat=0.00001,abstol=1e-16,reltol=1e-16,callback=cb4)

##
PO1_Time=Array(sol_saddle1.t)
PO1_Data=Array(sol_saddle1)

##
EngVal1=HamiltonianEDP_Vectorized(PO1_Data,p_nf_edp)
EngVal2=zeros(size(PO1_Data,2))

DumTheta2=zeros(size(PO1_Data,2))
for i=1:size(PO1_Data,2)
    EngVal2[i]=HamiltonianEDP(PO1_Data[:,i],p_nf_edp)
    DumTheta2[i]=DPE_Cal_dTheta(EngTarget,PO1_Data[:,i],p_nf_edp,1;pos_or_neg=2)
end

##
PO2_Time=Array(sol_saddle2.t)
PO2_Data=Array(sol_saddle2)

##
EngVal3=HamiltonianEDP_Vectorized(PO2_Data,p_nf_edp)
EngVal4=zeros(size(PO2_Data,2))

for i=1:size(PO2_Data,2)
    EngVal4[i]=HamiltonianEDP(PO2_Data[:,i],p_nf_edp)
end

## Plot the solution and double check its correct
PeriodPlot=Figure()
ax1=Axis(PeriodPlot[1,1],xlabel="theta1",ylabel="dtheta1")
Makie.lines!(ax1,PO1_Data[1,:],PO1_Data[3,:])
display(PeriodPlot)

##
PeriodPlot=Figure()
ax1=Axis(PeriodPlot[1,1],xlabel="Time",ylabel="dtheta1")
Makie.lines!(ax1,PO1_Time,PO1_Data[1,:])
#Makie.lines!(ax1,PO1_Time,DumTheta2)
display(PeriodPlot)

## 
EngPlot=Figure()
ax1=Axis(EngPlot[1,1],xlabel="Time",ylabel="Energy")
Makie.lines!(ax1,PO1_Time,EngVal1)
Makie.lines!(ax1,PO1_Time,EngVal2)
display(EngPlot)

##
PeriodPlot=Figure()
ax1=Axis3(PeriodPlot[1,1],xlabel="theta1",ylabel="dtheta1",zlabel="dtheta2")
Makie.lines!(ax1,PO1_Data[1,:],PO1_Data[3,:],PO1_Data[4,:])
display(PeriodPlot)

## Generate the animation of the double pendulum 
# Deinfe animation parameters
gap=50

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
path="Figures/EDoublePendulum_PO1_DownUp_Eng_$(EngTarget)_Animate.mp4";
AnimationParameter=[p_nf_edp[1],p_nf_edp[2],0.17272,0.2286,p_nf_edp[end]]
DP_Animate(DpPlot,ax1,PO1_Data,p_nf_edp,gap,color1,color2,color3,cmap,mks1,mks2,ls1,Ratio,path,framerate,xlims1,xlims2,ylims1,ylims2)

##
PeriodPlot=Figure()
ax1=Axis(PeriodPlot[1,1],xlabel="theta2",ylabel="dtheta2")
Makie.lines!(ax1,PO2_Data[2,:],PO2_Data[4,:])
display(PeriodPlot)

##
PeriodPlot=Figure()
ax1=Axis3(PeriodPlot[1,1],xlabel="theta2",ylabel="dtheta1",zlabel="dtheta2")
Makie.lines!(ax1,PO2_Data[2,:],PO2_Data[3,:],PO2_Data[4,:])
display(PeriodPlot)

##
EngVal2=HamiltonianEDP_Vectorized(PO2_Data,p_nf_edp)

## 
EngPlot=Figure()
ax1=Axis(EngPlot[1,1],xlabel="Time",ylabel="Energy")
Makie.lines!(ax1,PO2_Time,EngVal3)
Makie.lines!(ax1,PO2_Time,EngVal4)
display(EngPlot)

## Generate the animation of the double pendulum 
# Deinfe animation parameters
gap=50

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
path="Figures/EDoublePendulum_PO2_DownUp_Eng_$(EngTarget)_Animate.mp4";
AnimationParameter=[p_nf_edp[1],p_nf_edp[2],0.17272,0.2286,p_nf_edp[end]]
DP_Animate(DpPlot,ax1,PO2_Data,p_nf_edp,gap,color1,color2,color3,cmap,mks1,mks2,ls1,Ratio,path,framerate,xlims1,xlims2,ylims1,ylims2)

## Finally, save the solution for later use
@save "Results/EDP_L1PO_Eng_$(EngTarget)_Tube.jld2" p_nf_edp x_inc_L1 PO1_Time PO1_Data EngEq1 EngEq2 EngEq3 EngEq4
@save "Results/EDP_L2PO_Eng_$(EngTarget)_Tube.jld2" p_nf_edp x_inc_L2 PO2_Time PO2_Data EngEq1 EngEq2 EngEq3 EngEq4


