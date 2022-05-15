"""
# This fille will be use to generate the anamation of the energy cut
# Coded By: KK
# Last Updated: 05/15/2022
"""
##-----------------------------
## First load necessary files and Pacakages
@time include("Functions/LoadPackages.jl")
@time include("Functions/TubeFunctionTBP.jl")

## Calcultae the enegy level
μ=0.05;par=[1-μ;μ];

xs = LinRange(-2,2, 500)
ys = LinRange(-2,2, 500)

zs = [U([x,y,0],par) for x in xs , y in ys]
zs[zs.<-2.5].=NaN

Sun=load("Image/Sun.png")
Jupiter=load("Image/Jupiter.png")

# EngVal to use -1.8,-1.7,-1.65,-1.54,-1.5
TBP_Plot=Scene(resolution=(2560,1960))
streams=VideoStream(TBP_Plot, framerate = 20)

for EngVal=-2.5:0.01:-1.5
    p1 = Point3(-2, -2, EngVal)
    p2 = Point3(-2, 2, EngVal)
    p3 = Point3(2, -2, EngVal)
    p4 = Point3(2, 2, EngVal)
    
    TBP_Plot=Makie.mesh!(Sphere(Point3f0([-μ,0,-1.5]),0.3f0),color =earth,shading=false,show_axis = false)
    TBP_Plot=Makie.mesh!(Sphere(Point3f0([1-μ,0,-1.5]),0.1f0),color =moon,shading=false,show_axis = false)
    TBP_Plot=Makie.surface!(xs,ys,zs,shading=false,show_axis = false,colormap=cmap)
    TBP_Plot=Makie.meshscatter!([EqP1],[0],[U([EqP1,0,0,0],par)],markersize = mk1, color =cl1)
    TBP_Plot=Makie.meshscatter!([EqP2],[0],[U([EqP2,0,0,0],par)],markersize = mk1, color =cl1)
    TBP_Plot=Makie.meshscatter!([EqP3],[0],[U([EqP3,0,0,0],par)],markersize = mk1, color =cl1)
    TBP_Plot=Makie.meshscatter!([-μ+0.5],[sqrt(3)/2],[U([-μ+0.5,sqrt(3)/2,0,0],par)],markersize = mk1, color =cl1)
    TBP_Plot=Makie.meshscatter!([-μ+0.5],[-sqrt(3)/2],[U([-μ+0.5,-sqrt(3)/2,0,0],par)],markersize = mk1, color =cl1)
    TBP_Plot=Makie.mesh!([p1, p2, p3], color = (:grey,0.8), shading = false, transparency = true)
    TBP_Plot=Makie.mesh!([p4, p2, p3], color = (:grey,0.8), shading = false, transparency = true)
    update_cam!(TBP_Plot, Vec3f0(6, -6, 5), Vec3f0(0))
    
    # Recored the frame
    recordframe!(streams)
    
    # Clean the plot without creating a new one
    while length(TBP_Plot) > 0
      delete!(TBP_Plot, TBP_Plot[end])
    end
end

save("Figures/TBP_EnergyCut_Animation.mp4",streams)

