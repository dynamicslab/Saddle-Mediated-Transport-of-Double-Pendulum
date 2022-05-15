"""
# This file is used to plot the L2 Homoclinic orbits
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
PDataL2_unstable_pos=GetPoincareData(Traj1PO2)
PDataL2_unstable_neg=GetPoincareData(Traj2PO2)
PDataL2_stable_pos=GetPoincareData(Traj3PO2)
PDataL2_stable_neg=GetPoincareData(Traj4PO2)

## Calculate the index of trajectory
L2IndexUnstable=findmax(Int64.(abs.(PDataL2_unstable_pos[3,:]).<=1e-5))
L2IndexStable=findmax(Int64.(abs.(PDataL2_stable_pos[3,:]).<=1e-5))

## Get the corresponding trajectory
HomoL2Data1=Array(Traj1PO2[L2IndexUnstable[2]])
HomoL2Data2=Array(Traj3PO2[L2IndexStable[2]])

# Now plot the symetric Homiclinic orbits
ls=0.5

HomoPlotL2=Plots.plot(size=(400,400),tickfontsize=tickfontsize,framestyle=:box)
HomoPlotL2=Plots.plot!(HomoL2Data1[1,:],HomoL2Data1[2,:],color=col_p,linewidth=3*ls,label="")
HomoPlotL2=Plots.plot!(HomoL2Data2[1,:],HomoL2Data2[2,:],color=col_o,linewidth=3*ls,label="")

HomoPlotL2=Plots.scatter!(HomoPlotL2,[par[1]],[0],label="",
    marker = (:circle, mks1, 0.99, :red, stroke(1.0, 0.9, :black, :solid)))
HomoPlotL2=Plots.plot!(sin.(collect(0:0.01:2*pi)),cos.(collect(0:0.01:2*pi)),
    color=:black,linewidth=2*ls,linestyle=:dash,label="")
savefig(HomoPlotL2,"Figures/TBP_HomoclinicL2.pdf")
display(HomoPlotL2)

## Note: In this lower energy case, the unsymmetric orbits are too close to the symetric one. 
# In order to see more significant changes in the orbits, you need to use a different energy level than -1.5175 
## To determine the unsymetric L2 Homoclinic orbits, we need to find the intersection of Poincare cut
# Select few data points near the intersection
# dummy1=PDataL2_unstable_pos[:,(PDataL2_unstable_pos[1,:].>-2.0605).*(PDataL2_unstable_pos[1,:].<-2.05905).*(PDataL2_unstable_pos[3,:].>0)]
# dummy2=PDataL2_stable_pos[:,(PDataL2_stable_pos[1,:].>-2.0605).*(PDataL2_stable_pos[1,:].<-2.05905).*(PDataL2_stable_pos[3,:].>0)]

# DummyPlot=Scene()
# DummyPlot=Makie.scatter!(dummy1[1,:],dummy1[3,:])
# DummyPlot=Makie.scatter!(dummy2[1,:],dummy2[3,:])
# display(DummyPlot)

# # Try to solve the first line
# line_model1 = Model(Ipopt.Optimizer)
# @variable(line_model1, slope1)
# @variable(line_model1, lk1)
# @objective(line_model1, Min, sum((dummy1[3,:]-slope1.*dummy1[1,:].-lk1).^2))
# optimize!(line_model1)
# slope1_val=JuMP.value(slope1)
# lk1_val=JuMP.value(lk1)

# # Then solve for the second line
# line_model2 = Model(Ipopt.Optimizer)
# @variable(line_model2, slope2)
# @variable(line_model2, lk2)
# @objective(line_model2, Min, sum((dummy2[3,:]-slope2.*dummy2[1,:].-lk2).^2))
# optimize!(line_model2)
# slope2_val=JuMP.value(slope2)
# lk2_val=JuMP.value(lk2)

# # Then solve for the intersections
# x_intersect=(lk1_val-lk2_val)/(slope2_val-slope1_val)
# dx_intersect=slope1_val*x_intersect+lk1_val

# # Now calculate the vy
# UnsymetricL2Homo=CalVy([x_intersect,0,dx_intersect,0],EngTarget,par)

# ## Now simulate the homoclinic point
# tspan=(0.0,20.0)
# tspanminus=(0.0,-120.0)
# cb_pos1 = ContinuousCallback(condition,affect!,nothing,save_positions = (true,true))
# cb_pos2 = ContinuousCallback(condition,nothing,nothing,save_positions = (true,true))
# TBProb=ODEProblem(TBPODE!,UnsymetricL2Homo,tspan,par)
# sol1=DifferentialEquations.solve(TBProb,Vern9(),
#     saveat=0.001,tstops=T,abstol=1e-16,reltol=1e-16,callback=cb_pos1)
# TBProb=ODEProblem(TBPODE!,UnsymetricL2Homo,tspanminus,par)
# sol2=DifferentialEquations.solve(TBProb,Vern9(),
#     saveat=0.001,tstops=T,abstol=1e-16,reltol=1e-16,callback=cb_pos2)

# sol_data1=Array(sol1)
# sol_data2=Array(sol2)

# DummyPlot=Scene(resolution=(800,800))
# DummyPlot=Makie.lines!(sol_data1[1,:],sol_data1[2,:])
# DummyPlot=Makie.lines!(sol_data2[1,:],sol_data2[2,:])
# display(DummyPlot)







