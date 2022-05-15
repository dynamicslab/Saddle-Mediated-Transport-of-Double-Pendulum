"""
# This file is used to define the plooting setting, including the colors, resolutions, ect.
# Coded By: KK
# Last Updated: 05/15/2022
"""
## Define the different parameters here
# Orange color: Used to represent the unstable manifold
col_o=colorant"rgba(245,136,18,1)";
# Purple color: Used to represent the stable manifold
col_p=colorant"rgba(132,22,96,1)";

# The color map used to plot the enegy surface & Hill's region
cmap=ColorScheme(range(col_p,col_o));

# Define the ticker size
tickfontsize=20;

# Define the resolutions of plot
resol1=(650,400);
resol2=(2600,1600);

ENV["FREETYPE_ABSTRACTION_FONT_PATH"] ="C:\\WINDOWS\\Fonts"

# If you need to save to *.png, use following setting
#mks=14 # This set up the size of scatter 
#TickSize=125 # This set up the size of x and y ticks
#SpineWidth=6 # This set up the size of bounding box

MyTheme = GLMakie.Theme(fontsize = 125, resolution=resol2, linewidth=4, strokewidth=0, markersize=8,
    font="times")

GLMakie.set_theme!(MyTheme)
