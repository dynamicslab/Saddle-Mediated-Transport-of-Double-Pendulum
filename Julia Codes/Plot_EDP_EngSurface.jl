"""
# This file is used to plot the energy surface of double pendulum, with larger range on x and y axis
# Coded By: KK
# Last Updated: 05/15/2022
"""
##-----------------------------
## First load necessary files and Pacakages
@time include("Functions/LoadPackages.jl")
@time include("Functions/TubeFunctionDP.jl")
@time include("Functions/PlotSettings.jl")

# Define the simulation parameters
T=15.0;tspan=(0.0,T);dt=0.001;

# Calculate the energy at different equlibrium points
EngEq1=HamiltonianEDP([0,0,0,0],p_nf_edp)
EngEq2=HamiltonianEDP([0,pi,0,0],p_nf_edp)
EngEq3=HamiltonianEDP([pi,0,0,0],p_nf_edp)
EngEq4=HamiltonianEDP([pi,pi,0,0],p_nf_edp)

println("The energy of down-down is: ",EngEq1)
println("The energy of down-up is: ",EngEq2)
println("The energy of up-down is: ",EngEq3)
println("The energy of up-up is: ",EngEq4)

## Plot the energy surface of the double pendulum
Eng=EngEq3+0.1
dA=1
xs = LinRange(-pi-dA,pi+dA, 1000)
ys = LinRange(-pi-dA,pi+dA, 1000)
zs = [HamiltonianEDP([x,y,0,0],p_nf_edp) for x in xs, y in ys]

p1 = Point3(-pi-dA, -pi-dA, Eng)
p2 = Point3(pi+dA, pi+dA, Eng)
p3 = Point3(-pi-dA, pi+dA, Eng)
p4 = Point3(pi+dA, -pi-dA, Eng)

mk1=0.08
cl1=:blue
ticksize=34

##
#EngPlot=Scene()
# Plot the enegy surface
EngPlot=Makie.surface(xs,ys, zs,transparency = false, shading=true, show_axis=true,colormap=cmap,
     colorrange=(HamiltonianEDP([0,0,0,0],p_nf_edp),HamiltonianEDP([pi,pi,0,0],p_nf_edp)))
# Plot the equalibrium point
Makie.meshscatter!([0],[0],[HamiltonianEDP([0,0,0,0],p_nf_edp)],markersize = mk1, color =cl1,shading=false)
Makie.meshscatter!([0],[pi],[HamiltonianEDP([0,pi,0,0],p_nf_edp)],markersize = mk1, color =cl1,shading=false)
Makie.meshscatter!([pi],[0],[HamiltonianEDP([pi,0,0,0],p_nf_edp)],markersize = mk1, color =cl1,shading=false)
Makie.meshscatter!([pi],[pi],[HamiltonianEDP([pi,pi,0,0],p_nf_edp)],markersize = mk1, color =cl1,shading=false)
Makie.meshscatter!([0],[-pi],[HamiltonianEDP([0,pi,0,0],p_nf_edp)],markersize = mk1, color =cl1,shading=false)
Makie.meshscatter!([-pi],[0],[HamiltonianEDP([pi,0,0,0],p_nf_edp)],markersize = mk1, color =cl1,shading=false)
Makie.meshscatter!([-pi],[-pi],[HamiltonianEDP([pi,pi,0,0],p_nf_edp)],markersize = mk1, color =cl1,shading=false)
Makie.meshscatter!([-pi],[pi],[HamiltonianEDP([-pi,pi,0,0],p_nf_edp)],markersize = mk1, color =cl1,shading=false)
Makie.meshscatter!([pi],[-pi],[HamiltonianEDP([pi,-pi,0,0],p_nf_edp)],markersize = mk1, color =cl1,shading=false)
# # Plot the cut plane
Makie.mesh!([p1, p2, p3], color = (:grey, 0.8), transparency = true, shading=false)
Makie.mesh!([p1, p2, p4], color = (:grey, 0.8), transparency = true, shading=false)
#Makie.xlabel!("theta_1")
# Save plot
# save("Figures/DP_EngSurface.png", EngPlot,resolution = (2560,1960)) 
# Show
display(EngPlot)











