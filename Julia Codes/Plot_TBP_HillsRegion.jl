"""
# This fille will be use to plot the Hill's region at different level
# Coded By: KK
# Last Updated: 05/15/2022
"""
##-----------------------------
## First load necessary files and Pacakages
@time include("Functions/LoadPackages.jl")
@time include("Functions/TubeFunctionTBP.jl")
@time include("Functions/PlotSettings.jl")

## Calcultae the enegy level
μ=0.05;par=[1-μ;μ];

xs = LinRange(-2,2, 500)
ys = LinRange(-2,2, 500)

zs = [U([x,y,0],par) for x in xs , y in ys]
zs[zs.<-2.5].=NaN

Sun=load("Image/Sun.png")
Jupiter=load("Image/Jupiter.png")

# EngVal to use -1.8,-1.7,-1.65,-1.54,-1.5
for EngVal in [-1.8,-1.7,-1.65,-1.54,-1.5]
    # Get the plane ploints
    p1 = Point3(-2, -2, EngVal)
    p2 = Point3(-2, 2, EngVal)
    p3 = Point3(2, -2, EngVal)
    p4 = Point3(2, 2, EngVal)

    TBP_Plot=Scene()
    #TBP_Plot=Makie.mesh!(Sphere(Point3f0([-μ,0,-1.5]),0.3f0),color =Sun,shading=true,show_axis = false)
    #TBP_Plot=Makie.mesh!(Sphere(Point3f0([1-μ,0,-1.5]),0.1f0),color =Jupiter,shading=true,show_axis = false)
    TBP_Plot=Makie.surface!(xs,ys,zs,shading=false,show_axis = false,colormap=cmap)
    TBP_Plot=Makie.meshscatter!([EqP1],[0],[0],markersize = mk1, color =cl1)
    TBP_Plot=Makie.meshscatter!([EqP2],[0],[0],markersize = mk1, color =cl1)
    TBP_Plot=Makie.meshscatter!([EqP3],[0],[0],markersize = mk1, color =cl1)
    TBP_Plot=Makie.meshscatter!([-μ+0.5],[sqrt(3)/2],[0],markersize = mk1, color =cl1) 
    TBP_Plot=Makie.meshscatter!([-μ+0.5],[-sqrt(3)/2],[0],markersize = mk1, color =cl1)
    TBP_Plot=Makie.mesh!([p1, p2, p3], color = (:white,1), shading = false, transparency = true)
    TBP_Plot=Makie.mesh!([p4, p2, p3], color = (:white,1), shading = false, transparency = true)
    cam2d!(TBP_Plot)
    
    #TBP_Plot.center = false # prevent to recenter on display
    #save("Figures/3BP_Hill_Eng_$(EngVal).png", TBP_Plot,resolution = (800,800)) 
    display(TBP_Plot)
end