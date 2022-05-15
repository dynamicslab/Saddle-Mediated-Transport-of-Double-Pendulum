"""
# This file is used to plot the Heteroclinis orbits of TBP
# Coded By: KK
# Last Updated: 05/15/2022
"""
##-----------------------------
## First load necessary files and Pacakages
@time include("Functions/LoadPackages.jl")
@time include("Functions/TubeFunctionTBP.jl")
@time include("Functions/PlotSettings.jl")

## Load the Results
μ=9.537e-4;Gap=10;EngTarget=-1.5175;

@load "Results/TBP_Manifold_$(μ)_Gap_$(Gap)_Eng_$(EngTarget).jld2" PODataPO1 PODataPO2 μ par Traj1PO1 Traj2PO1 Traj3PO1 Traj4PO1 Traj1PO2 Traj2PO2 Traj3PO2 Traj4PO2

## Get the Poincare cut data
PDataL1_unstable_pos=GetPoincareData(Traj1PO1)#
PDataL1_unstable_neg=GetPoincareData(Traj2PO1)
PDataL1_stable_pos=GetPoincareData(Traj3PO1)
PDataL1_stable_neg=GetPoincareData(Traj4PO1)
#
PDataL2_unstable_pos=GetPoincareData(Traj1PO2)
PDataL2_unstable_neg=GetPoincareData(Traj2PO2)
PDataL2_stable_pos=GetPoincareData(Traj3PO2)
PDataL2_stable_neg=GetPoincareData(Traj4PO2)#

## Downsample the data if needed
# gap=10
# PDataL1_unstable_pos=PDataL1_unstable_pos[:,1:gap:end]
# PDataL2_stable_neg=PDataL2_stable_neg[:,1:gap:end]

## There's no  direct intersection of stable and unstable manifold at the first time. Thus, we need to simulate the final poincare cuts one more time.
# Now simulate all the initial conditions, and store the simulate values
# Define the ODE problem
tspan=(0.0,20.0)
tspanminus=(0.0,-10.0)
TBProb=ODEProblem(TBPODE!,zeros(4,1),tspan,par)
TBProbMinus=ODEProblem(TBPODE!,zeros(4,1),tspanminus,par)

# Prepare for ensemble problem for L1 point
function prob_func_L1(prob,i,repeat)
    remake(prob,u0=PDataL1_unstable_pos[:,i])
end

# Prepare for ensemble problem for L2 point
function prob_func_L2(prob,i,repeat)
    remake(prob,u0=PDataL2_stable_neg[:,i])
end

# Define ensemble problem
ensemble_prob_nf_L1 = EnsembleProblem(TBProb,prob_func=prob_func_L1)
ensemble_prob_nf_L2 = EnsembleProblem(TBProbMinus,prob_func=prob_func_L2)

## Define the event function
# When the particle hits the x axis, we will let the simulation stop
affect!(integrator) = terminate!(integrator)

#Define a continuous callback function
function condition2(u,t,integrator) 
    return u[1]-par[1]
end

function condition4(u,t,integrator) 
    if abs(u[2])<0.6
        return par[1]-u[1]
    else
        return false
    end
end

cb_2=ContinuousCallback(condition2,nothing,affect!,save_positions = (true,true))
cb_4=ContinuousCallback(condition2,affect!,nothing,save_positions = (true,true))

## Solve for ODE
L1Cut2 = DifferentialEquations.solve(ensemble_prob_nf_L1,
        Vern9(),EnsembleSerial(),callback=cb_2,
        saveat=0.001, save_start=true,save_end=true,
        abstol=1e-15,reltol=1e-15,trajectories=size(PDataL1_unstable_pos,2));

L2Cut2 = DifferentialEquations.solve(ensemble_prob_nf_L2,
        Vern9(),EnsembleSerial(),callback=cb_4,
        saveat=0.001, save_start=true,save_end=true,
        abstol=1e-15,reltol=1e-15,trajectories=size(PDataL2_stable_neg,2));

## Plot the manifold of the L1
ls=2
HamPlotL1=Scene(resolution=(1300,1300))
HamPlotL1=MakiePlotMonoTBP(HamPlotL1,L1Cut2,col_p,ls)

## Plot the Poincare cut at next iteration
L1Cut2Poincare=GetPoincareData(L1Cut2)
#
L1PoincareCut2=Scene(resolution=(1300,1300))
L1PoincareCut2=Makie.lines!(L1Cut2Poincare[2,:],L1Cut2Poincare[4,:])

## Plot the manifold of L2
ls=2
HamPlotL2=Scene(resolution=(1300,1300))
HamPlotL2=MakiePlotMonoTBP(HamPlotL2,L2Cut2,col_o,ls)

## Plot the Poincare cut at next iteration
L2Cut2Poincare=GetPoincareData(L2Cut2)
# 
L2PoincareCut2=Scene(resolution=(1300,1300))
L2PoincareCut2=Makie.lines!(L2Cut2Poincare[2,:],L2Cut2Poincare[4,:])

## Plot all poincare cut together
HeteroPoincare=Plots.plot(size=(400,400),tickfontsize=tickfontsize,framestyle=:box)
HeteroPoincare=Plots.plot!(PDataL1_unstable_pos[2,:],PDataL1_unstable_pos[4,:],
    color=col_p,label="",markersize=0.01)
HeteroPoincare=Plots.plot!(PDataL2_stable_neg[2,:],PDataL2_stable_neg[4,:],
    color=col_o,label="",markersize=0.01)
HeteroPoincare=Plots.plot!(L1Cut2Poincare[2,:],L1Cut2Poincare[4,:],
    color=col_p,label="",markersize=0.01)
HeteroPoincare=Plots.plot!(L2Cut2Poincare[2,:],L2Cut2Poincare[4,:],
    color=col_o,label="",markersize=0.01)
# HeteroPoincare=Plots.scatter!(PDataL1_unstable_neg[2,:],PDataL1_unstable_neg[4,:],
#     color=col_p,label="",markersize=0.01)
# HeteroPoincare=Plots.scatter!(PDataL2_stable_pos[2,:],PDataL2_stable_pos[4,:],
#     color=col_o,label="",markersize=0.01)
HeteroPoincare=Plots.hline!([0],color=:black,label="",ls=:dash)
HeteroPoincare=Plots.vline!([0],color=:black,label="",ls=:dash)
#
#savefig(HeteroPoincare,"Figures/TBP_HeteroPoincareCut2.pdf")
display(HeteroPoincare)

## Plot all poincare cut together
HeteroPoincare=Plots.plot(size=(400,400),tickfontsize=tickfontsize,framestyle=:box)
HeteroPoincare=Plots.plot!(L1Cut2Poincare[2,:],L1Cut2Poincare[4,:],
    color=col_p,label="",markersize=0.01)
HeteroPoincare=Plots.plot!(L2Cut2Poincare[2,:],L2Cut2Poincare[4,:],
    color=col_o,label="",markersize=0.01)
# HeteroPoincare=Plots.scatter!(PDataL1_unstable_neg[2,:],PDataL1_unstable_neg[4,:],
#     color=col_p,label="",markersize=0.01)
# HeteroPoincare=Plots.scatter!(PDataL2_stable_pos[2,:],PDataL2_stable_pos[4,:],
#     color=col_o,label="",markersize=0.01)
HeteroPoincare=Plots.hline!([0],color=:black,label="",ls=:dash)
HeteroPoincare=Plots.vline!([0],color=:black,label="",ls=:dash)
#
#savefig(HeteroPoincare,"Figures/TBP_HeteroPoincareCut2ZoomIn.pdf")
display(HeteroPoincare)

## Plot the Zoom In intersection
Plots.plot!(HeteroPoincare,xlim=(0.04,0.052),ylim=(-0.1,0.1),xticks = ([0.04,0.045,0.05], ["0.04", "0.045", "0.05"]))
savefig(HeteroPoincare,"Figures/TBP_HeteroPoincareCut2ZoomIn2.pdf")
display(HeteroPoincare)

## Now get the points in interest
L1Intersect=L1Cut2Poincare[:,(L1Cut2Poincare[2,:].>0.04337).*(L1Cut2Poincare[2,:].<0.04339).*(L1Cut2Poincare[4,:].<0.01).*(L1Cut2Poincare[4,:].>-0.01)]
L2Intersect=L2Cut2Poincare[:,(L2Cut2Poincare[2,:].>0.04337).*(L2Cut2Poincare[2,:].<0.04339).*(L2Cut2Poincare[4,:].<0.01).*(L2Cut2Poincare[4,:].>-0.01)]

DummyPlot=Scene()
DummyPlot=Makie.scatter!(L1Intersect[2,:],L1Intersect[4,:],color=col_p)
DummyPlot=Makie.scatter!(L2Intersect[2,:],L2Intersect[4,:],color=col_o)
display(DummyPlot)

## Now extract the tracjetories
IndexL1=findmax(L1Cut2Poincare[2,:].==L1Intersect[2,:])[2]
IndexL2=findmax(L2Cut2Poincare[2,:].==L2Intersect[2,:])[2]

TrajL1Cut1=Array(Traj1PO1[IndexL1])
TrajL1Cut2=Array(L1Cut2[IndexL1])
TrajL2Cut1=Array(Traj4PO2[IndexL2])
TrajL2Cut2=Array(L2Cut2[IndexL2])

## Finally, plot it
ls=0.5

HeteroPlot=Plots.plot(size=(400,400),tickfontsize=tickfontsize,framestyle=:box)
# Plot the L1 part
HeteroPlot=Plots.plot!(TrajL1Cut1[1,:],TrajL1Cut1[2,:],color=col_p,linewidth=3*ls,label="")
HeteroPlot=Plots.plot!(TrajL1Cut2[1,:],TrajL1Cut2[2,:],color=col_p,linewidth=3*ls,label="")
# Plot the L2 part
HeteroPlot=Plots.plot!(TrajL2Cut1[1,:],TrajL2Cut1[2,:],color=col_o,linewidth=3*ls,label="")
HeteroPlot=Plots.plot!(TrajL2Cut2[1,:],TrajL2Cut2[2,:],color=col_o,linewidth=3*ls,label="")
# Plot all L1 and L2 points, and the Jupiter
HeteroPlot=Plots.scatter!(HeteroPlot,[EqP1],[0],label="",
    marker = (:circle, mks1, 0.99, :blue, stroke(1.0, 0.9, :black, :solid)))
HeteroPlot=Plots.scatter!(HeteroPlot,[EqP2],[0],label="",
    marker = (:circle, mks1, 0.99, :blue, stroke(1.0, 0.9, :black, :solid)))
HeteroPlot=Plots.scatter!(HeteroPlot,[par[1]],[0],label="", 
    marker = (:circle, 2*mks1, 0.99, :red, stroke(1.0, 0.9, :black, :solid)))
HeteroPlot=Plots.plot!(sin.(collect(0:0.01:2*pi)),cos.(collect(0:0.01:2*pi)),
    color=:black,linewidth=2*ls,linestyle=:dash,label="")
HeteroPlot=Plots.plot!(HeteroPlot,xlim=(par[1]-0.1,par[1]+0.1),ylim=(-0.1,0.1))
savefig(HeteroPlot,"Figures/TBP_HomoclinicTraj.pdf")
display(HeteroPlot)

