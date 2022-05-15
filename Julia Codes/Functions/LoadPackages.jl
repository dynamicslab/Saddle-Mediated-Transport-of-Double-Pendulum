"""
# This file is used to load all the necessary Pacakages to run the program
# Coded By: KK
# Last Updated: 05/15/2022
"""

# Pacakages used for solving differential equation
using DifferentialEquations
using DiffEqFlux
using ModelingToolkit

# Pacakages used for Neural Network and Optimization
using JuMP
using Ipopt
using Flux
using DiffEqFlux
using Zygote
using ForwardDiff

# Other Pacakages
using LinearAlgebra
using JLD2
using FileIO
using ImageFiltering
using NBInclude
using GeometryBasics
using ControlSystems
using Optim
using DiffEqSensitivity

# Pacakages used for plotting
#using Makie
using GLMakie
#set_theme!(theme_black())
GLMakie.activate!()
#using AbstractPlotting
using Plots
pgfplotsx()
using Colors
using ColorSchemes
using Latexify
using LaTeXStrings
using PolygonOps
using StaticArrays

## The following code will using PyCall to import the CasaDi optimization Pacakages
#using PyCall
#using Conda

# Use the following code to set-up the Python environment on the machine
# PyCall.ENV["PYTHON"] = "C:\\Users\\kahdi\\Anaconda3\\envs\\tensorflow\\python.exe"
# Use the following code to use the Python come with Julia
# PyCall.ENV["PYTHON"] = ""
# Use this line to enable the pip_install
# Conda.pip_interop(true)
# Add the MPC Pacakageswe need
# Conda.pip("install", "casadi")
# Conda.pip("install", "do-mpc")

