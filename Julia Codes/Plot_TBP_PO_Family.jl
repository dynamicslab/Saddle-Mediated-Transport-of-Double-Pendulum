"""
# This file is used to plot varies periodic orbits of TBP porblem
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
mks1=2.5
ls1=0.5

# Plot the L1 and L2 points
pp1=Plots.scatter([EqP1,EqP2],[0,0],markersize=mks1,label="",size=resol1,tickfontsize=tickfontsize,
    marker = (:circle, mks1, 0.99, :blue, stroke(1.0, 0.9, :black, :solid)))

# Plot the remaining trajectory
@time for i=1:N_Inc
    # Extract the orbits
    Time=Array(PoincareSim_nf_L1[i].t)
    Data=Array(PoincareSim_nf_L1[i])
    
    Data=[Data[:,1:10:end] Data[:,end]]
    
    # Plot 
    pp1=Plots.plot!(Data[1,:],Data[2,:],linewidth=ls1,label="",color=col_p,
        size=resol1,tickfontsize=tickfontsize,framestyle=:box)
end

@time for i=1:N_Inc
    # Extract the orbits
    Time=Array(PoincareSim_nf_L2[i].t)
    Data=Array(PoincareSim_nf_L2[i])
    
    Data=[Data[:,1:10:end] Data[:,end]]
    
    # Plot 
    pp1=Plots.plot!(Data[1,:],Data[2,:],linewidth=ls1,label="",color=col_p,
        size=resol1,tickfontsize=tickfontsize,framestyle=:box)
end

#savefig(pp1,"Figures/TBP_PO_Family.pdf")

display(pp1)



