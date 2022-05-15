"""
# This file is used to plot the energy surface of TBP problem
# Coded By: KK
# Last Updated: 05/15/2022
"""
##-----------------------------
## First load necessary files and Pacakages
@time include("Functions/LoadPackages.jl")
@time include("Functions/TubeFunctionTBP.jl")
@time include("Functions/PlotSettings.jl")

## Define the mass ratio to use
μ=0.05;par=[1-μ;μ];

## Calculate the energy
xs = LinRange(-2,2, 500)
ys = LinRange(-2,2, 500)

zs = [U([x,y,0],par) for x in xs , y in ys]
zs[zs.<-2.5].=NaN

EqP1=FindL1(μ)
EqP2=FindL2(μ)
EqP3=FindL3(μ)

## Load the texture
Sun=load("Image/Sun.png")
Jupiter=load("Image/Jupiter.png")

## Plot
TBP_Plot=Scene()
mk1=0.05
cl1=:blue
TBP_Plot=Makie.mesh!(Sphere(Point3f0([-μ,0,-1.5]),0.3f0),color =Sun,shading=false,show_axis = false)
TBP_Plot=Makie.mesh!(Sphere(Point3f0([1-μ,0,-1.5]),0.1f0),color =Jupiter,shading=false,show_axis = false)
TBP_Plot=Makie.surface!(xs,ys,zs,shading=false,show_axis = false,colormap=cmap)
TBP_Plot=Makie.meshscatter!([EqP1],[0],[U([EqP1,0,0,0],par)],markersize = mk1, color =cl1)
TBP_Plot=Makie.meshscatter!([EqP2],[0],[U([EqP2,0,0,0],par)],markersize = mk1, color =cl1)
TBP_Plot=Makie.meshscatter!([EqP3],[0],[U([EqP3,0,0,0],par)],markersize = mk1, color =cl1)
TBP_Plot=Makie.meshscatter!([-μ+0.5],[sqrt(3)/2],[U([-μ+0.5,sqrt(3)/2,0,0],par)],markersize = mk1, color =cl1)
TBP_Plot=Makie.meshscatter!([-μ+0.5],[-sqrt(3)/2],[U([-μ+0.5,-sqrt(3)/2,0,0],par)],markersize = mk1, color =cl1)
update_cam!(TBP_Plot,Vec3f0(5, -6, 3), Vec3f0(0))
TBP_Plot.center = false # prevent to recenter on display
save("Figures/TBP_Engsurface.png", TBP_Plot,resolution = (2560,1960)) 
display(TBP_Plot)

## Plot the energy cut illustration
EngVal=-1.8
p1 = Point3(-2, -2, EngVal)
p2 = Point3(-2, 2, EngVal)
p3 = Point3(2, -2, EngVal)
p4 = Point3(2, 2, EngVal)

TBP_Plot=Scene()
TBP_Plot=Makie.mesh!(Sphere(Point3f0([-μ,0,-1.5]),0.3f0),color =Sun,shading=false,show_axis = false)
TBP_Plot=Makie.mesh!(Sphere(Point3f0([1-μ,0,-1.5]),0.1f0),color =Jupiter,shading=false,show_axis = false)
TBP_Plot=Makie.surface!(xs,ys,zs,shading=false,show_axis = false,colormap=cmap)
TBP_Plot=Makie.meshscatter!([EqP1],[0],[U([EqP1,0,0,0],par)],markersize = mk1, color =cl1)
TBP_Plot=Makie.meshscatter!([EqP2],[0],[U([EqP2,0,0,0],par)],markersize = mk1, color =cl1)
TBP_Plot=Makie.meshscatter!([EqP3],[0],[U([EqP3,0,0,0],par)],markersize = mk1, color =cl1)
TBP_Plot=Makie.meshscatter!([-μ+0.5],[sqrt(3)/2],[U([-μ+0.5,sqrt(3)/2,0,0],par)],markersize = mk1, color =cl1)
TBP_Plot=Makie.meshscatter!([-μ+0.5],[-sqrt(3)/2],[U([-μ+0.5,-sqrt(3)/2,0,0],par)],markersize = mk1, color =cl1)
TBP_Plot=Makie.mesh!([p1, p2, p3], color = (:grey,0.8), shading = false, transparency = true)
TBP_Plot=Makie.mesh!([p4, p2, p3], color = (:grey,0.8), shading = false, transparency = true)
update_cam!(TBP_Plot, Vec3f0(3, -5, 3), Vec3f0(0))
#TBP_Plot.center = false # prevent to recenter on display
save("Figures/TBP_EngsurfaceCut.png", TBP_Plot,resolution = (2560,1960)) 
display(TBP_Plot)



