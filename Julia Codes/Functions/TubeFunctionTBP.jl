"""
# This file is used to define the all the functions needed for the PCR3BP calculation.
# Coded By: KK
# Last Updated: 05/15/2022
"""

## Load Pacakages
using Ipopt
using JuMP
using Zygote
using DiffEqFlux
using DifferentialEquations

# Include those function
MonoMatrix=include("MonoMatrix.jl")

## This function will be used to calculate the potential energy of
# the CRTBP (Circular restricted three body problem). Coded using 
# Eq.2.3.10 in book of Shane Ross, r1 and r2 calculate using Eq in page 32
function U(q,par)
    r1square=(q[1]+par[2])^2+q[2]^2
    r2square=(q[1]-par[1])^2+q[2]^2
    Unum=-0.5*(par[1]*r1square+par[2]*r2square)-(par[1]/sqrt(r1square))-(par[2]/sqrt(r2square))
    return Unum
end

## This function will be used to calculate the Hamiltonian energy of three body problem
function Heng(u,par)
    x=u[1]
    y=u[2]
    vx=u[3]
    vy=u[4]

    r1=sqrt(((x-par[1]).^2)+(y.^2))
    r2=sqrt(((x+par[2]).^2)+(y.^2))
    
    EngVal=0.5*(vx^2+vy^2)+U(u,par)
    
    return EngVal
end

## This function will be used to simulate the EOM of CRTBP
function TBPODE!(du,u,p,t)
    r1=sqrt(u[2]^2+(u[1]-p[1])^2)
    r2=sqrt(u[2]^2+(u[1]+p[2])^2)
    
    du[1]=u[3]
    du[2]=u[4]
    du[3]=2*u[4]+u[1]-(u[1]-p[1])*p[2]/(r1^3)-p[1]*(u[1]+p[2])/(r2^3)
    du[4]=-2*u[3]+u[2]-(u[2]*p[2])/(r1)^(3)-p[1]*u[2]/(r2)^(3)
    
    return du
end

## The following functions will calculate the three langrangian points
# given the mass ratio of CRTBP
function FindL1(μ)
    rh=(μ/3)^(1/3)
    xg=rh*(1-(1/3)*rh-(1/9)*rh^2-(1/27)*rh^3)
    optModel= Model(Ipopt.Optimizer)
    JuMP.set_optimizer_attribute(optModel, "print_level", 0)
    @variable(optModel, xopt, start = xg)
    @NLobjective(optModel,Min,(xopt^5-(3-μ)*xopt^4+(3-2*μ)*xopt^3-μ*xopt^2+2*μ*xopt-μ)^2)
    JuMP.optimize!(optModel);
    
    return (1-μ)-value(xopt)
end

function FindL2(μ)
    rh=(μ/3)^(1/3)
    xg=rh*(1+(1/3)*rh+(1/9)*rh^2+(1/27)*rh^3)
    optModel= Model(Ipopt.Optimizer)
    JuMP.set_optimizer_attribute(optModel, "print_level", 0)
    @variable(optModel, xopt, start = xg)
    @NLobjective(optModel,Min,(xopt^5+(3-μ)*xopt^4+(3-2*μ)*xopt^3-μ*xopt^2-2*μ*xopt-μ)^2)
    JuMP.optimize!(optModel);
    
    return (1-μ)+value(xopt)
end

function FindL3(μ)
    optModel= Model(Ipopt.Optimizer)
    JuMP.set_optimizer_attribute(optModel, "print_level", 0)
    @variable(optModel, xopt, start =-1)
    @NLobjective(optModel,Min,(μ*((xopt+μ)-1)/abs((xopt+μ)-1)^3-(μ-1)*(xopt+μ)/abs(xopt+μ)^3-xopt)^2)
    JuMP.optimize!(optModel);
    
    return value(xopt)
end

## The following functions will be used to calculate the 
# periodict orbits of the CRTBP

# This function will calculate the velocity of the body given energy level
# and postion
function CalVoptEng(EngLevel,x,par)
    x_pos=x[1]

    μ1=par[1]
    μ2=par[2]
    
    pA=2*μ2/abs(x_pos-μ1)
    pB=2*μ1/abs(x_pos+μ2)
    
    Vopt=sqrt(2*EngLevel+x_pos^2+pA+pB+μ1*μ2)
    
    return Vopt
end

# This function will calculate the loss function of given 
# initial guess of the starting position of the three body problem
# Define the loss function we are trying to minimize
# par=[1-μ,μ]
function TBPlossVal(x,par,EngTarget,tspan,cb)
    x_pos=x[1]
    # Calculate the velocity based on the position and energy
    vy_opt=CalVoptEng(EngTarget,x_pos,par)
    # Define the initial energy
    x0g=[x_pos,0,0,vy_opt]
    TBProbPO=ODEProblem(TBPODE!,x0g,tspan,par)

    sol=Array(DifferentialEquations.solve(TBProbPO,Vern9(),
        dt=dt,tstops=T,abstol=1e-16,reltol=1e-16,
        callback=cb,save_everystep = false,save_start=false,save_end=true))

    loss=abs(sol[2])+abs(sol[3])
    return loss
end

# This function will be used to calculate the periodic orbits of CRTBP
function OptPOTBP(x,par,EngTarget,tspan,cb,N,tol,δ)
    x_pos=x[1]
    pin=1
    # Create a dummy function
    lossVal(x_pos)=TBPlossVal(x_pos,par,EngTarget,tspan,cb)

    while lossVal(x_pos)>tol && pin<N
        dx_pos= Zygote.forwarddiff(lossVal,x_pos)[1]
        #dvy= Zygote.gradient(lossVal,vy_opt)[1]
        val=lossVal(x_pos)
        x_pos=x_pos-δ*dx_pos

        println("The loss value is :$(val) at iteration $(pin)...")
        pin=pin+1
    end

    # Calculate the fianl intial conditions
    vy_opt=CalVoptEng(EngTarget,x_pos,par)
    x_po=[x_pos,0,0,vy_opt]

    if lossVal(x_pos)>tol
        println("Need more iterations to converge...")
    else
        println("The solution converged! Periodic orbits been found!")
        println("The identified initial contitions of periodic orbits is: ",x_po)
    end

    return x_po
end

# Define the augmented new ode
function MonoODE!(du,u,p,t)
    # The first four states are the normal ode terms
    x_n=u[1:nx]
    du[1:nx]=TBPODE!(similar(x_n),x_n,p,t)
    
    # The remaining terms comes from the augmented states
    x_m=reshape(u[nx+1:end],nx,nx)
    dMono=MonoMatrix(x_n,p)*x_m
    
    # Now assign the new derivative
    du[nx+1:end]=reshape(dMono,nx*nx,1)
    
    return du
end

# The following function will simulate multiple initial condition X0,
# and then calculate the manifold
function SimulateManifoldTBP(X0,par,backward,cb)
    # Now generate an EnsembleProblem to simulate all the initial conditions 
    # Define ODE problems
    Tsim=25.0
    dt=0.001
    tspan=(0.0,Tsim)
    if backward==1
        dt=-0.001
        tspan=(0.0,-Tsim)
    end
    TBProb=ODEProblem(TBPODE!,X0[:,1],tspan,par)
    
    # Prepare for ensemble problem
    function prob_func(prob,i,repeat)
      remake(prob,u0=X0[:,i])
    end

    # Define ensemble problem
    ensemble_prob_nf = EnsembleProblem(TBProb,prob_func=prob_func)

    # Now calcultae and plot the poincare map
    println("\t Calculating the manifold...")
    NC=size(X0,2)
    
    PoincareSim_nf=[]
    
    @time PoincareSim_nf = DifferentialEquations.solve(ensemble_prob_nf,
        Vern9(),EnsembleThreads(),callback=cb,
        saveat=0.01, save_start=true,save_end=true,
        abstol=1e-15,reltol=1e-15,trajectories=NC);

    return PoincareSim_nf
end

# The following function will plot the Manifold using Makie
function MakiePlotMonoTBP(MonoPlot,Traj,color,ls;gap=10)
    for i=1:size(Traj,3)
        TrajData=Array(Traj[i])
        lines!(MonoPlot,TrajData[1,1:gap:end],TrajData[2,1:gap:end],color=color,linewidth=ls)
    end
    
    return nothing
end

# The following function will plot the Manifold using Plots
function PlotsPlotMonoTBP(MonoPlot,Traj,color,ls)
    for i=1:size(Traj,3)
        TrajData=Array(Traj[i])
        MonoPlot=Plots.plot!(TrajData[1,:],TrajData[2,:],
            color=color,linewidth=ls,label="")
    end
    
    return MonoPlot
end

# The following function will get the Poincare Cut data of emsembled trajectory
function GetPoincareData(Traj)
    Nnum=size(Traj,3)
    DataCut=zeros(4,Nnum)

    for i=1:Nnum
        Data=Array(Traj[i])
        DataCut[:,i]=Data[:,end]
    end

    return DataCut
end

# Calculate the velocity in y direction given other three variables and energy
function CalVy(xbar,EngTarget,par)
    x=xbar[1]
    y=xbar[2]
    dx=xbar[3]

    μ1=par[1]
    μ2=par[2]

    part1=x^2-dx^2+y^2+2*EngTarget+μ1*μ2
    part2=2*μ2/sqrt(y^2+(x-μ1)^2)
    part3=2*μ1/sqrt(y^2+(x+μ2)^2)

    dy=sqrt(part1+part2+part3)

    x_new=[x,y,dx,dy]

    return x_new
end

# Calculate the velocity in x direction given other three variables and energy
function CalVx(xbar,EngTarget,par)
    x=xbar[1]
    y=xbar[2]
    dy=xbar[4]

    μ1=par[1]
    μ2=par[2]

    part1=x^2-dy^2+y^2+2*EngTarget+μ1*μ2
    part2=2*μ2/sqrt(y^2+(x-μ1)^2)
    part3=2*μ1/sqrt(y^2+(x+μ2)^2)

    dx=sqrt(part1+part2+part3)

    x_new=[x,y,dx,dy]

    return x_new
end




