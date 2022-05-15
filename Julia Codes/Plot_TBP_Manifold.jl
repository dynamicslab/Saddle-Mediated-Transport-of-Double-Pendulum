"""
# This file is used to plot the L1 and L2 PO manifold.
# Coded By: KK
# Last Updated: 05/15/2022
"""
##-----------------------------
## First load necessary files and Pacakages
@time include("Functions/LoadPackages.jl")
@time include("Functions/TubeFunctionTBP.jl")
@time include("Functions/PlotSettings.jl")

## Load the Results
μ=9.537e-4;Gap=1000;EngTarget=-1.5175;

@load "Results/TBP_Manifold_$(μ)_Gap_$(Gap)_Eng_$(EngTarget).jld2" PODataPO1 PODataPO2 μ par Traj1PO1 Traj2PO1 Traj3PO1 Traj4PO1 Traj1PO2 Traj2PO2 Traj3PO2 Traj4PO2
## Now plot the L1 manifold tubes, that generates the homoclinic orbits and orbits inside the Hill's region
ls=0.5

HomoPlotL1=Plots.plot(size=(400,400),tickfontsize=tickfontsize,framestyle=:box)
HomoPlotL1=PlotsPlotMonoTBP(HomoPlotL1,Traj2PO1,col_p,ls)
HomoPlotL1=PlotsPlotMonoTBP(HomoPlotL1,Traj3PO1,col_o,ls)
# HomoPlotL1=PlotsPlotMonoTBP(HomoPlotL1,Traj1PO1,col_p,ls)
# HomoPlotL1=PlotsPlotMonoTBP(HomoPlotL1,Traj4PO1,col_o,ls)
HomoPlotL1=Plots.plot!(PODataPO1[1,:],PODataPO1[2,:],color=:black,linewidth=0.5*ls,label="")
HomoPlotL1=Plots.plot!(sin.(collect(0:0.01:2*pi)),cos.(collect(0:0.01:2*pi)),
    color=:black,linewidth=2*ls,linestyle=:dash,label="")
savefig(HomoPlotL1,"Figures/TBP_HomoTubeL1.pdf")
display(HomoPlotL1)

## Now plot the L2 manifold tubes, that generates the homoclinic orbits and orbits outside the Hill's region
HomoPlotL2=Plots.plot(size=(400,400),tickfontsize=tickfontsize,framestyle=:box)
HomoPlotL2=PlotsPlotMonoTBP(HomoPlotL2,Traj1PO2,col_p,ls)
HomoPlotL2=PlotsPlotMonoTBP(HomoPlotL2,Traj3PO2,col_o,ls)
# HomoPlotL2=PlotsPlotMonoTBP(HomoPlotL2,Traj2PO2,col_p,ls)
# HomoPlotL2=PlotsPlotMonoTBP(HomoPlotL2,Traj4PO2,col_o,ls)
HomoPlotL2=Plots.plot!(PODataPO2[1,:],PODataPO2[2,:],color=:black,linewidth=0.5*ls,label="")
HomoPlotL2=Plots.plot!(sin.(collect(0:0.01:2*pi)),cos.(collect(0:0.01:2*pi)),
    color=:black,linewidth=2*ls,linestyle=:dash,label="")
savefig(HomoPlotL2,"Figures/TBP_HomoTubeL2.pdf")
display(HomoPlotL2)

## Plot the zoom in details of the L1 manifold near the saddle
mks1=2.5
HomoPlotL1ZoomIn=HomoPlotL1
HomoPlotL1ZoomIn=Plots.plot(HomoPlotL1ZoomIn,
    xlim=(EqP1-0.1,EqP1+0.1),ylim=(-0.1,0.1))
HomoPlotL1ZoomIn=Plots.scatter!(HomoPlotL1ZoomIn,[EqP1],[0],label="",
    marker = (:circle, mks1, 0.99, :blue, stroke(1.0, 0.9, :black, :solid)))
savefig(HomoPlotL1ZoomIn,"Figures/TBP_HomoTubeL1ZoomIn.pdf")
display(HomoPlotL1ZoomIn)

## Plot the zoom in details of the L2 manifold near the saddle
HomoPlotL2ZoomIn=HomoPlotL2
HomoPlotL2ZoomIn=Plots.plot(HomoPlotL2ZoomIn,
    xlim=(EqP2-0.1,EqP2+0.1),ylim=(-0.1,0.1))
HomoPlotL2ZoomIn=Plots.scatter!(HomoPlotL2ZoomIn,[EqP2],[0],label="",
    marker = (:circle, mks1, 0.99, :blue, stroke(1.0, 0.9, :black, :solid)))
savefig(HomoPlotL2ZoomIn,"Figures/TBP_HomoTubeL2ZoomIn.pdf")
display(HomoPlotL2ZoomIn)

## Plot the zoom in plot of L1 and L2 tubes intersection
HeteroPlot=Plots.plot(size=(400,400),tickfontsize=tickfontsize,framestyle=:box,
    xlim=(EqP1-0.05,EqP2+0.05),ylim=(-0.05,0.05))
HeteroPlot=PlotsPlotMonoTBP(HeteroPlot,Traj1PO1,col_p,ls)
HeteroPlot=PlotsPlotMonoTBP(HeteroPlot,Traj4PO1,col_o,ls)
HeteroPlot=PlotsPlotMonoTBP(HeteroPlot,Traj2PO2,col_p,ls)
HeteroPlot=PlotsPlotMonoTBP(HeteroPlot,Traj4PO2,col_o,ls)
# HeteroPlot=PlotsPlotMonoTBP(HeteroPlot,Traj2PO1,col_p,ls)
# HeteroPlot=PlotsPlotMonoTBP(HeteroPlot,Traj3PO1,col_o,ls)
# HeteroPlot=PlotsPlotMonoTBP(HeteroPlot,Traj1PO2,col_p,ls)
# HeteroPlot=PlotsPlotMonoTBP(HeteroPlot,Traj3PO2,col_o,ls)
HeteroPlot=Plots.plot!(PODataPO1[1,:],PODataPO1[2,:],color=:black,linewidth=0.5*ls,label="")
HeteroPlot=Plots.plot!(PODataPO2[1,:],PODataPO2[2,:],color=:black,linewidth=0.5*ls,label="")
HeteroPlot=Plots.plot!(sin.(collect(0:0.01:2*pi)),cos.(collect(0:0.01:2*pi)),
    color=:black,linewidth=2*ls,linestyle=:dash,label="")
HeteroPlot=Plots.scatter!(HeteroPlot,[EqP1],[0],label="",
    marker = (:circle, mks1, 0.99, :blue, stroke(1.0, 0.9, :black, :solid)))
HeteroPlot=Plots.scatter!(HeteroPlot,[EqP2],[0],label="",
    marker = (:circle, mks1, 0.99, :blue, stroke(1.0, 0.9, :black, :solid)))

savefig(HeteroPlot,"Figures/TBP_HeteroTubeZoomIn.pdf")
display(HeteroPlot)

## Now plot all the manifold together
TubePlot=Plots.plot(size=(400,400),tickfontsize=tickfontsize,framestyle=:box)
TubePlot=PlotsPlotMonoTBP(TubePlot,Traj1PO1,col_p,ls)
TubePlot=PlotsPlotMonoTBP(TubePlot,Traj2PO1,col_p,ls)
TubePlot=PlotsPlotMonoTBP(TubePlot,Traj3PO1,col_o,ls)
TubePlot=PlotsPlotMonoTBP(TubePlot,Traj4PO1,col_o,ls)
TubePlot=Plots.plot!(PODataPO1[1,:],PODataPO1[2,:],color=:black,linewidth=0.5*ls,label="")
#
TubePlot=PlotsPlotMonoTBP(TubePlot,Traj1PO2,col_p,ls)
TubePlot=PlotsPlotMonoTBP(TubePlot,Traj2PO2,col_p,ls)
TubePlot=PlotsPlotMonoTBP(TubePlot,Traj3PO2,col_o,ls)
TubePlot=PlotsPlotMonoTBP(TubePlot,Traj4PO2,col_o,ls)
TubePlot=Plots.plot!(PODataPO2[1,:],PODataPO2[2,:],
    color=:black,linewidth=0.5*ls,label="")
#
TubePlot=Plots.plot!(sin.(collect(0:0.01:2*pi)),cos.(collect(0:0.01:2*pi)),
    color=:black,linewidth=2*ls,linestyle=:dash,label="")
#
savefig(TubePlot,"Figures/TBP_TubePlot.pdf")
display(TubePlot)

## Now plot zoom in plot of all the tubes
TubePlotZoomIn=Plots.plot(TubePlot,xlim=(EqP1-0.3,EqP2+0.3),ylim=(-0.2,0.2))
TubePlotZoomIn=Plots.scatter!(TubePlotZoomIn,[EqP1],[0],label="",
    marker = (:circle, mks1, 0.99, :blue, stroke(1.0, 0.9, :black, :solid)))
TubePlotZoomIn=Plots.scatter!(TubePlotZoomIn,[EqP2],[0],label="",
    marker = (:circle, mks1, 0.99, :blue, stroke(1.0, 0.9, :black, :solid)))
savefig(TubePlotZoomIn,"Figures/TBP_TubePlotZoomIn.pdf")
display(TubePlotZoomIn)


