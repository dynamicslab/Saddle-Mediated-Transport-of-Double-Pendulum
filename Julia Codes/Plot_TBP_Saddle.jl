"""
# This file is used to plot the vector filed of the saddles of PCR3BP.
# Coded By: KK
# Last Updated: 05/15/2022
"""
##-----------------------------
## First load necessary files and Pacakages
@time include("Functions/LoadPackages.jl")
@time include("Functions/TubeFunctionTBP.jl")
@time include("Functions/PlotSettings.jl")

## Define plotting parameters
μ=0.05;par=[1-μ;μ];dt=0.001;

## Plot the first saddle point
dA=0.05
xs = LinRange(EqP1-dA,EqP1+dA, 100)
ys = LinRange(-dA,dA, 100)
zs = [U([x,y,0],par) for x in xs , y in ys];

XYrange=0.189
xa = LinRange(EqP1-dA,EqP1+dA, 10)
ya = LinRange(-dA,dA, 10)
za = [U([x,y,0],par) for x in xa, y in ya];

u,v=ImageFiltering.imgradients(za, KernelFactors.ando3)
XA=repeat(xa,size(ya,1)) |> vec
YA=repeat(ya,size(xa,1)) |> vec

SaddlePlot=Figure()
ax1=Axis(SaddlePlot[1,1])
Makie.contour!(ax1,xs, ys, zs,levels=20,linewidth=8,show_axis=false,colormap=cmap,colorrange=(minimum(zs),maximum(zs)))
Makie.arrows!(ax1,xa, ya,u, v,arrowsize = 0.005, arrowcolor = :black,lengthscale =0.005,linewidth=12,normalize=true)
Makie.scatter!(ax1,[EqP1], [0],color=:blue,markersize=20,strokewidth=8)
Makie.xlims!(ax1,EqP1-dA,EqP1+dA)
Makie.ylims!(ax1,-dA,dA)
# Save plot
#save("Figures/TBP_Saddle1.png",SaddlePlot,resolution = (800,800)) 
# Show
display(SaddlePlot)

## Plot the second saddle
dA=0.05
xs = LinRange(EqP2-dA,EqP2+dA, 10)
ys = LinRange(-dA,dA, 10)
zs = [U([x,y,0],par) for x in xs , y in ys];

XYrange=0.189
xa = LinRange(EqP2-dA,EqP2+dA, 10)
ya = LinRange(-dA,dA, 10)
za = [U([x,y,0],par) for x in xa, y in ya];

u,v=ImageFiltering.imgradients(za, KernelFactors.ando3)
XA=repeat(xa,size(ya,1)) |> vec
YA=repeat(ya,size(xa,1)) |> vec

SaddlePlot=Scene()
SaddlePlot=Makie.contour!(xs, ys, zs,levels=10,linewidth=8,show_axis=false,colormap=cmap,
    colorrange=(minimum(zs),maximum(zs)))
SaddlePlot=Makie.arrows!(xa, ya,
    u, v,arrowsize = 0.005, arrowcolor = :black,lengthscale =0.005,linewidth=12,normalize=true)
SaddlePlot=Makie.scatter!([EqP2], [0],color=:blue,markersize=20,strokewidth=8)
Makie.xlims!(EqP2-dA,EqP2+dA)
Makie.ylims!(-dA,dA)
# Save plot
#save("Figures/TBP_Saddle2.png",SaddlePlot,resolution = (800,800)) 
# Show
display(SaddlePlot)





