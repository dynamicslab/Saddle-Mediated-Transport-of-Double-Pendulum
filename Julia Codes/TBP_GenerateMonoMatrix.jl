"""
# This file is used to generate the monodromy matrix
# Coded By: KK
# Last Updated: 05/15/2022
"""
##-----------------------------
## First load necessary files and Pacakages
@time include("Functions/LoadPackages.jl")
@time include("Functions/TubeFunctionTBP.jl")

# First create a symbolic variable and calculate the Jacobian matrix
nx=4
@ModelingToolkit.variables x_var[1:nx]
@ModelingToolkit.variables mu_var[1:2]
# Next, get the expression of the function RHS
t_var=0
RHS=TBPODE!(similar(x_var),x_var,mu_var,t_var)

# Then calculate the jacobian
DRHS=ModelingToolkit.jacobian(RHS,x_var)

# Transfer the calculated result into function
write("Functions/MonoMatrix.jl",string(build_function(DRHS,x_var,mu_var)[1]))
    
# Include those function
MonoMatrix=include("Functions/MonoMatrix.jl")
