"""
# This file is used to plot illustration of TBP porblem's periodic orbits 
# Coded By: KK
# Last Updated: 05/15/2022
"""
##-----------------------------
## First load necessary files and Pacakages
@time include("Functions/LoadPackages.jl")
@time include("Functions/TubeFunctionTBP.jl")
@time include("Functions/PlotSettings.jl")
## Define parameter and load files
μ=9.537e-4; N_Inc=10;

@load "Results/L1PO_$(μ)_Ninc_$(N_Inc).jld2" par POL1 PoincareSim_nf_L1 EqP1 EqP2 EqP3 EngP1 EngP2 EngP3
@load "Results/L2PO_$(μ)_Ninc_$(N_Inc).jld2" par POL2 PoincareSim_nf_L2 EqP1 EqP2 EqP3 EngP1 EngP2 EngP3

## Now start plotting
# Define the size of marker
mks1=6
ls1=0.5

Time=Array(PoincareSim_nf_L1[5].t)
Data=Array(PoincareSim_nf_L1[5])

Data=[Data[:,1:10:end] Data[:,end]]
# Plot the L1 and L2 points
pp1=Plots.scatter([1-μ],[0],color=:grey,label="",size=(800,200),tickfontsize=20,
    marker = (:circle, 10, 0.99, :crimson, stroke(2, 0.8, :black, :solid))) 

pp1=Plots.scatter!([EqP1,EqP2],[0,0],markersize=mks1,label="",size=resol1,tickfontsize=tickfontsize,
    marker = (:circle, mks1, 0.99, :blue, stroke(1.0, 0.9, :black, :solid)))

# Extract the orbits
pp1=Plots.plot!(Data[1,:],Data[2,:],color=col_p,linewidth=3,label="",size=(800,200),tickfontsize=20,framestyle=:none)

savefig("Figures/TBP_POL1_Illustrate.pdf")

##
Time=Array(PoincareSim_nf_L2[8].t)
Data=Array(PoincareSim_nf_L2[8])

Data=[Data[:,1:10:end] Data[:,end]]
# Plot the L1 and L2 points
pp1=Plots.scatter([1-μ],[0],color=:grey,label="",size=(800,200),tickfontsize=20,
    marker = (:circle, 10, 0.99, :crimson, stroke(2, 0.8, :black, :solid))) 

pp1=Plots.scatter!([EqP1,EqP2],[0,0],markersize=mks1,label="",size=resol1,tickfontsize=tickfontsize,
    marker = (:circle, mks1, 0.99, :blue, stroke(1.0, 0.9, :black, :solid)))

# Extract the orbits
pp1=Plots.plot!(Data[1,:],Data[2,:],color=col_p,linewidth=3,label="",size=(800,200),tickfontsize=20,framestyle=:none)

savefig("Figures/TBP_POL2_Illustrate.pdf")



