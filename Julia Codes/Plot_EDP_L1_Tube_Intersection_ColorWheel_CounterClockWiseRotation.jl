"""
# This file is used to plot the tubes of double pendulum's L1 saddle points (Using experimental parameters), and generate a color wheel
# Coded By: KK
# Last Updated: 05/15/2022
"""
##-----------------------------
## First load necessary files and Pacakages
@time include("Functions/LoadPackages.jl")
@time include("Functions/TubeFunctionDP.jl")
@time include("Functions/PlotSettings.jl")

## Define the loading parameters
EngTarget=0.2;Tsim=3;SizeInc=9448;

## Load the POData
@load "Results/EDP_L1PO_Eng_$(EngTarget)_Tube.jld2" p_nf_edp x_inc_L1 PO1_Time PO1_Data EngEq1 EngEq2 EngEq3 EngEq4
@load "Results/EDP_L2PO_Eng_$(EngTarget)_Tube.jld2" p_nf_edp x_inc_L2 PO2_Time PO2_Data EngEq1 EngEq2 EngEq3 EngEq4

## load the calculation result
@time @load "Results/EDP_Eng_$(EngTarget)_Tsim_$(Tsim)_POIncNum_$(SizeInc)_F_P_UN.jld2" TrajVarInc1_pos_unstable_NoFric Inc1_1
@time @load "Results/EDP_Eng_$(EngTarget)_Tsim_$(Tsim)_POIncNum_$(SizeInc)_F_N_UN.jld2" TrajVarInc1_neg_unstable_NoFric Inc1_2
@time @load "Results/EDP_Eng_$(EngTarget)_Tsim_$(Tsim)_POIncNum_$(SizeInc)_B_P_S.jld2" TrajVarInc1_pos_stable_NoFric Inc1_3
@time @load "Results/EDP_Eng_$(EngTarget)_Tsim_$(Tsim)_POIncNum_$(SizeInc)_B_N_S.jld2" TrajVarInc1_neg_stable_NoFric Inc1_4

## Load saved cut data
@load "Results/EDP_Swipe_L1_Tube_Eng.jld2" StoredCutDataL1 EngListL1

## Use this if you want to try high energy level
Inc_neg_unstable=ExtractIncFromEndPoint(TrajVarInc1_neg_unstable_NoFric) 
Inc_neg_stable=ExtractIncFromEndPoint(TrajVarInc1_neg_stable_NoFric) 

[HamiltonianEDP(Inc_neg_stable[:,i],p_nf_edp) for i=1:size(Inc_neg_stable,2)]


## Use this if you want to try low energy level
#println("The energy level we are using now is: ",HamiltonianEDP(StoredCutDataL1[2],p_nf_edp))
#Inc_neg_unstable=-StoredCutDataL1[2]
#Inc_neg_stable=zeros(size(Inc_neg_unstable))
#Inc_neg_stable[1,:]=-Inc_neg_unstable[1,:]
#Inc_neg_stable[2,:]=-Inc_neg_unstable[2,:]
#Inc_neg_stable[3,:]=Inc_neg_unstable[3,:]
#Inc_neg_stable[4,:]=Inc_neg_unstable[4,:]

#[HamiltonianEDP(Inc_neg_stable[:,i],p_nf_edp) for i=1:size(Inc_neg_stable,2)]

## Now plot the end point we extracted from the previous calculation
mks=3
TubePlot=Figure(resolution=(1600,1600))
ax1=Axis(TubePlot[1,1],xlabel="θ1",ylabel="dθ1")
Makie.scatter!(ax1,Inc_neg_unstable[1,:],Inc_neg_unstable[3,:],color=col_p,strokewidth=0,markersize=mks)
Makie.scatter!(ax1,Inc_neg_stable[1,:],Inc_neg_stable[3,:],color=col_o,strokewidth=0,markersize=mks)
ax1.aspect=AxisAspect(1/1)
display(TubePlot)

## Once we made sure the tube Poincare cut is correct, let's start create a Polygon 
AxisNum1=1;AxisNum2=3
# Generate polygon
PolyPoinsX=vcat(Inc_neg_unstable[AxisNum1,:],Inc_neg_unstable[AxisNum1,1])
PolyPoinsY=vcat(Inc_neg_unstable[AxisNum2,:],Inc_neg_unstable[AxisNum2,1])
polygon_pos_unstable=SVector.(PolyPoinsX,PolyPoinsY)

# Generate stable polygon
PolyPoinsXs=vcat(Inc_neg_stable[AxisNum1,:],Inc_neg_stable[AxisNum1,1])
PolyPoinsYs=vcat(Inc_neg_stable[AxisNum2,:],Inc_neg_stable[AxisNum2,1])
polygon_pos_stable=SVector.(PolyPoinsXs,PolyPoinsYs)

## Sample point from the region
# Use this one for Eng=0.2
#xa = -2:0.03:2
#ya = -20.5:0.1:-2
# Use this one for Eng=-0.06985971252329815
xa = -1.8:0.01:1.8
ya = 3:0.05:20.5

Points = vec(SVector.(xa',ya))

insideIndexPolyCurve1 = [PolygonOps.inpolygon(p, polygon_pos_unstable; in=true, on=false, out=false) for p in Points]
insideIndexPolyCurve2 = [PolygonOps.inpolygon(p, polygon_pos_stable; in=true, on=false, out=false) for p in Points]

PointsPolyCurve1=Points[insideIndexPolyCurve1]
PointsPolyCurve2=Points[insideIndexPolyCurve2]

# Get the intersecting point
intersectIndexPolyCurve1 = [PolygonOps.inpolygon(p, polygon_pos_unstable; in=true, on=false, out=false) for p in PointsPolyCurve2]
Not_intersectIndexPolyCurve1=.!intersectIndexPolyCurve1
NotIntersecIncCurve1=PointsPolyCurve2[Not_intersectIndexPolyCurve1]

# The final points that lives in the intersection
IntersectPoints=PointsPolyCurve2[intersectIndexPolyCurve1]
IntersectPointsDummy=PointsPolyCurve2[intersectIndexPolyCurve1]

## Now design a color wheel
σx=1;σy=1

IntersectPointsNomalized=zeros(size(IntersectPointsDummy,1),2)

for i=1:size(IntersectPoints)[1]
    IntersectPointsNomalized[i,1]=IntersectPointsDummy[i][1]
    IntersectPointsNomalized[i,2]=IntersectPointsDummy[i][2]
end

# Normalize the data
IntersectPointsNomalized[:,1]=IntersectPointsNomalized[:,1]./maximum(abs.(IntersectPointsNomalized[:,1]))
IntersectPointsNomalized[:,2]=IntersectPointsNomalized[:,2]./maximum(abs.(IntersectPointsNomalized[:,2]))

IntersectPointsNomalized=vec(SVector.(IntersectPointsNomalized[:,1],IntersectPointsNomalized[:,2]))
##
#ColorWheel=GetColorWheel(IntersectPoints,[0,-15.361143713745697],σx,σy,HalfAngle=false);
ColorWheel=GetColorWheel(IntersectPointsNomalized,[0,0],σx,σy,HalfAngle=true);

## Plot the points 
mks=14
TickSize=125
SpineWidth=6
PolyPlot=Figure(resolution=(2600,1600))
#ax1=Axis(PolyPlot[1,1],xlabel="θ1",ylabel="dθ1",xlabelfont="DejaVu Sans")
ax1=Axis(PolyPlot[1,1],xticklabelfont="times",yticklabelfont="times",xticklabelsize=TickSize,yticklabelsize=TickSize,spinewidth=6)
Makie.scatter!(ax1,Inc_neg_unstable[1,:],Inc_neg_unstable[3,:],color=col_p,strokewidth=0,markersize=mks)
Makie.scatter!(ax1,Inc_neg_stable[1,:],Inc_neg_stable[3,:],color=col_o,strokewidth=0,markersize=mks)
Makie.scatter!(ax1,IntersectPoints,markersize=14,color=ColorWheel)
Makie.scatter!(ax1,[0],[15.361143713745697],markersize=14,color="black")
save("Figures/EDP_L1_ColorWheel_Sim0_CCW.png", PolyPlot,px_per_unit=100,pt_per_unit=100)
display(PolyPlot)

## Now find out the intial condition using energy relation
IncSim=zeros(4,size(IntersectPoints,1))
angle_index=2
for i=1:size(IntersectPoints,1)
    IncSim[1,i]=IntersectPoints[i][1]
    IncSim[3,i]=IntersectPoints[i][2]
    IncSim[4,i]=DPE_Cal_dTheta(EngTarget,IncSim[:,i],p_nf_edp,angle_index;pos_or_neg=0)
    #IncSim[4,i]=DPE_Cal_dTheta(EngTarget,IncSim[:,i],p_nf_edp,angle_index;pos_or_neg=1)
end

display(IncSim)

##
PolyPlot=Scene()
Makie.scatter!(PolyPlot,Inc_neg_unstable[1,:],Inc_neg_unstable[3,:],Inc_neg_unstable[4,:],color=col_p,strokewidth=0,markersize=1*mks)
Makie.scatter!(PolyPlot,Inc_neg_stable[1,:],Inc_neg_stable[3,:],Inc_neg_stable[4,:],color=col_o,strokewidth=0,markersize=1*mks)
Makie.scatter!(PolyPlot,IncSim[[1,3,4],:],markersize=1*mks,color=ColorWheel)
Makie.xlabel!("theta_1")
Makie.zlabel!("dtheta_1")
Makie.zlabel!("dtheta_2")
Makie.scale!(PolyPlot,1,0.1,0.1)
display(PolyPlot)

##
[HamiltonianEDP(IncSim[:,i],p_nf_edp) for i=1:size(IncSim,2)]

## Now simulate all those points forward in time, and stop the simulation when it reaches θ2=2π
# When the particle hits the x axis, we will let the simulation stop
function affect_cut!(integrator) 
    terminate!(integrator)
end

# Define the callback function
function cut_function1!(u,t,integrator) 
    # Define when to stop when theta1 hits the target
    return u[2]+2*pi
end

cb_cut1=ContinuousCallback(cut_function1!,affect_cut!,affect_cut!,save_positions = (true,false));

## Now simulate all those points 
Tsim=2.0;everystepsave=false;backward=0
Gap=1;tspan=(0.0,Tsim);
FinalSimInc=IncSim[:,1:Gap:end]

## Do not use EnsembleThreads() approach, which will crash the VS code
TrajInTube_End=SimulateManifold(FinalSimInc,p_nf_edp,backward,Tsim,everystepsave,cb_cut1;saveend=true,saveat=false,experimental=true)

## Get the final position of each point
TrajInTube_Data=ExtractIncFromEndPoint(TrajInTube_End)

## Determine how many points lives inside the intersection
Sim1Points = vec(SVector.(TrajInTube_Data[1,TrajInTube_Data[3,:].<0],TrajInTube_Data[3,TrajInTube_Data[3,:].<0]))
InsideIndex1 = [PolygonOps.inpolygon(p, polygon_pos_stable; in=true, on=false, out=false) for p in Sim1Points]
Sim1PointsInside=Sim1Points[InsideIndex1]

## Plot the points 
mks=14
TickSize=125
SpineWidth=6
PolyPlot=Figure(resolution=(2600,1600))
#ax1=Axis(PolyPlot[1,1],xlabel="θ1",ylabel="dθ1",xlabelfont="DejaVu Sans")
ax1=Axis(PolyPlot[1,1],xticklabelfont="times",yticklabelfont="times",xticklabelsize=TickSize,yticklabelsize=TickSize,spinewidth=6)
Makie.scatter!(ax1,Inc_neg_unstable[1,:],Inc_neg_unstable[3,:],color=col_p,strokewidth=0,markersize=mks)
Makie.scatter!(ax1,Inc_neg_stable[1,:],Inc_neg_stable[3,:],color=col_o,strokewidth=0,markersize=mks)
Makie.scatter!(ax1,TrajInTube_Data[[1,3],:],markersize=14,color=ColorWheel[1:Gap:end])
#Makie.scatter!(ax1,Sim1PointsInside,markersize=14,color="red")
Makie.scatter!(ax1,[0],[15.361143713745697],markersize=14,color="black")
save("Figures/EDP_L1_ColorWheel_Sim1_CCW.png",PolyPlot)
display(PolyPlot)

## Now simulate this again
TrajInTube_Data_NewInc=TrajInTube_Data
TrajInTube_Data_NewInc[2,:]=TrajInTube_Data_NewInc[2,:].+2*pi
TrajInTube_End_2=SimulateManifold(TrajInTube_Data_NewInc,p_nf_edp,backward,Tsim,everystepsave,cb_cut1;saveend=true,saveat=false,experimental=true)

## Get the final position of each point
# 38084 points for H=0.2, they all returned in the unstable manifold
TrajInTube_Data_2=ExtractIncFromEndPoint(TrajInTube_End_2)


## Extract the points that returned to the 2*pi position
TrajInTube_Data_2_Returned=TrajInTube_Data_2[:,(TrajInTube_Data_2[2,:].<-2*pi+0.01).*(TrajInTube_Data_2[2,:].>-2*pi-0.01).*(TrajInTube_Data_2[1,:].<pi+0.01).*(TrajInTube_Data_2[1,:].>-pi-0.01)]
ColorWheel_Returned=ColorWheel[(TrajInTube_Data_2[2,:].<-2*pi+0.01).*(TrajInTube_Data_2[2,:].>-2*pi-0.01).*(TrajInTube_Data_2[1,:].<pi+0.01).*(TrajInTube_Data_2[1,:].>-pi-0.01)]

## Determine how many points lives inside the intersection
# For H=0.2, 10768 points lies inside the intersection regions
Sim2Points = vec(SVector.(TrajInTube_Data_2_Returned[1,TrajInTube_Data_2_Returned[3,:].>0],TrajInTube_Data_2_Returned[3,TrajInTube_Data_2_Returned[3,:].>0]))
InsideIndex2 = [PolygonOps.inpolygon(p, polygon_pos_stable; in=true, on=false, out=false) for p in Sim2Points]
Sim2PointsInside=Sim2Points[InsideIndex2]

## Plot the points 
mks=14
TickSize=125
SpineWidth=6
PolyPlot=Figure(resolution=(2600,1600))
#ax1=Axis(PolyPlot[1,1],xlabel="θ1",ylabel="dθ1",xlabelfont="DejaVu Sans")
ax1=Axis(PolyPlot[1,1],xticklabelfont="times",yticklabelfont="times",xticklabelsize=TickSize,yticklabelsize=TickSize,spinewidth=6)
Makie.scatter!(ax1,Inc_neg_unstable[1,:],Inc_neg_unstable[3,:],color=col_p,strokewidth=0,markersize=mks)
Makie.scatter!(ax1,Inc_neg_stable[1,:],Inc_neg_stable[3,:],color=col_o,strokewidth=0,markersize=mks)
Makie.scatter!(ax1,TrajInTube_Data_2_Returned[[1,3],:],markersize=14,color=ColorWheel_Returned)
#Makie.scatter!(ax1,Sim2PointsInside,markersize=14,color="red")
Makie.scatter!(ax1,[0],[15.361143713745697],markersize=14,color="black")
save("Figures/EDP_L1_ColorWheel_Sim2_CCW.png",PolyPlot)
display(PolyPlot)

## Now simulate again
TrajInTube_Data_3_NewInc=TrajInTube_Data_2_Returned
TrajInTube_Data_3_NewInc[2,:]=TrajInTube_Data_3_NewInc[2,:].+2*pi
TrajInTube_End_3=SimulateManifold(TrajInTube_Data_3_NewInc,p_nf_edp,backward,Tsim,everystepsave,cb_cut1;saveend=true,saveat=false,experimental=true)

## Get the final position of each point
TrajInTube_Data_3=ExtractIncFromEndPoint(TrajInTube_End_3)

## Extract the points taht returned to the 2*pi position
# For H=0.2, 15297 points came back to the 2π intersections
TrajInTube_Data_3_Returned=TrajInTube_Data_3[:,(TrajInTube_Data_3[2,:].<-2*pi+0.01).*(TrajInTube_Data_3[2,:].>-2*pi-0.01).*(TrajInTube_Data_3[1,:].<pi+0.01).*(TrajInTube_Data_3[1,:].>-pi-0.01)]
ColorWheel_Returned_3=ColorWheel_Returned[(TrajInTube_Data_3[2,:].<-2*pi+0.01).*(TrajInTube_Data_3[2,:].>-2*pi-0.01).*(TrajInTube_Data_3[1,:].<pi+0.01).*(TrajInTube_Data_3[1,:].>-pi-0.01)]

## Determine how many points lives inside the intersection
# For H=0.2, there are 5823 points taht lives inside the intersection region
Sim3Points = vec(SVector.(TrajInTube_Data_3_Returned[1,TrajInTube_Data_3_Returned[3,:].>0],TrajInTube_Data_3_Returned[3,TrajInTube_Data_3_Returned[3,:].>0]))
InsideIndex3 = [PolygonOps.inpolygon(p, polygon_pos_stable; in=true, on=false, out=false) for p in Sim3Points]
Sim3PointsInside=Sim3Points[InsideIndex3]

## Plot the points 
mks=14
TickSize=125
SpineWidth=6
PolyPlot=Figure(resolution=(2600,1600))
#ax1=Axis(PolyPlot[1,1],xlabel="θ1",ylabel="dθ1",xlabelfont="DejaVu Sans")
ax1=Axis(PolyPlot[1,1],xticklabelfont="times",yticklabelfont="times",xticklabelsize=TickSize,yticklabelsize=TickSize,spinewidth=6)
Makie.scatter!(ax1,Inc_neg_unstable[1,:],Inc_neg_unstable[3,:],color=col_p,strokewidth=0,markersize=mks)
Makie.scatter!(ax1,Inc_neg_stable[1,:],Inc_neg_stable[3,:],color=col_o,strokewidth=0,markersize=mks)
Makie.scatter!(ax1,TrajInTube_Data_3_Returned[[1,3],:],markersize=14,color=ColorWheel_Returned_3)
Makie.scatter!(ax1,Sim3PointsInside,markersize=14,color="red")
Makie.scatter!(ax1,[0],[-15.361143713745697],markersize=14,color="black")
save("Figures/EDP_L1_ColorWheel_Sim3_CCW.png",PolyPlot)
display(PolyPlot)

## Run the periodic orbits forward in time
#Define a continuous callback function
function affect_cut!(integrator) 
    terminate!(integrator)
end
#
function conditionSaddle1(u,t,integrator) 
    return u[2]+pi
end
# Now define the callback function
cb1 = ContinuousCallback(conditionSaddle1,affect_cut!,affect_cut!,save_positions = (true,true))
# For H=−0.06985
#x_inc_L1=[0,pi,6.252561555256689,10.74784938633]
# For H=0.2
x_inc_L1=[0,pi,-9.65638354688606,-20.5152473222608]
DP_Prob_Saddle1=ODEProblem(DPE_NoFric!,x_inc_L1,(0,5.0),p_nf_edp)
#
sol_saddle1=DifferentialEquations.solve(DP_Prob_Saddle1,Vern9(),save_start=true,save_end=true,
    saveat=0.001,abstol=1e-16,reltol=1e-16,callback=cb1)
#
PO1_Time=Array(sol_saddle1.t)
PO1_Data=Array(sol_saddle1)

## Now plot this 2π periodic orbits, and save it
# Plot the positiveeps plot
PoincarPlot=Plots.plot(size=(650,400),tickfontsize=20)
PoincarPlot=Plots.plot!(PoincarPlot,PO1_Data[1,:],PO1_Data[3,:],color=col_p,
    label="",linewidth=2.5,framestyle=:box)
savefig(PoincarPlot,"Figures/EDP_2pi_PO1_1_vs_3_CCW.pdf") 

##
PoincarPlot=Plots.plot(size=(650,400),tickfontsize=20)
PoincarPlot=Plots.plot!(PoincarPlot,PO1_Data[2,:],PO1_Data[4,:],color=col_p,
    label="",linewidth=2.5,framestyle=:box)
savefig(PoincarPlot,"Figures/EDP_2pi_PO1_2_vs_4_CCW.pdf") 

##
PoincarPlot=Plots.plot(size=(650,400),tickfontsize=20,yticks = ([-π,0,π], ["-\\pi", "0", "\\pi"]),ylim=(-π-0.1,π+0.1))
PoincarPlot=Plots.plot!(PoincarPlot,PO1_Data[1,:],PO1_Data[2,:],color=col_p,
    label="",linewidth=2.5,framestyle=:box)
savefig(PoincarPlot,"Figures/EDP_2pi_PO1_1_vs_2_CCW.pdf") 
