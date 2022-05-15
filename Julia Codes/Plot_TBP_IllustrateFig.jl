"""
# This file will plot the skematic of TBP
# Coded By: KK
# Last Updated: 05/15/2022
"""
##-----------------------------
## First load necessary files and Pacakages
@time include("Functions/LoadPackages.jl")
@time include("Functions/TubeFunctionTBP.jl")
@time include("Functions/PlotSettings.jl")

## Define the simulation parameters
T=10.0;tspan=(0.0,T);μ=0.05;par=[1-μ;μ];dt=0.001;

# Define the initial condition and ODE problem
x0=[-1.7496856637955596,0,-0.00010505816385064136,1.0385268513422474]

TBProb=ODEProblem(TBPODE!,x0,tspan,par)
sol=DifferentialEquations.solve(TBProb,Vern9(),dt=dt,tstops=T,abstol=1e-16,reltol=1e-16)

Time=Array(sol.t)
Data=Array(sol)

## Plot the trajectory
pp1=Plots.scatter([-μ,1-μ],[0,0],color=:grey,label="",size=(800,400),tickfontsize=20,
    marker = (:circle, 10, 0.99, :crimson, stroke(2, 0.8, :black, :solid))) 
pp1=Plots.scatter!([EqP1,EqP2,EqP3,-μ+0.5,-μ+0.5],[0,0,0,sqrt(3)/2,-sqrt(3)/2],markersize=2,label="",size=(800,400),tickfontsize=20,
    marker = (:circle, 6, 0.99, :blue, stroke(2, 0.8, :black, :solid)))
pp1=Plots.plot!(Data[1,:],Data[2,:],linewidth=3,label="",size=(800,400),color=col_o,tickfontsize=20,framestyle=:none)

savefig(pp1,"Figures/3BP_Sktech.pdf")




