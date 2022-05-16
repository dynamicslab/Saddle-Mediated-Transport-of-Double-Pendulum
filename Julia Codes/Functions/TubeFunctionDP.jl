"""
# This file includes all the functions needed to calculate the P.O of double pendulum and it's manifold
# Coded By: KK
# Last Updated: 05/15/2022
"""

## Define the parameters of the experimental double pendulum
# m1,m2,a1,a2,L1,I1,I2,k1,k2,g
p_f_edp=[0.0938439747983407,0.137595970499412,0.108565214610396,0.116779018445656,0.172719203817519,0.000437529429756910,0.00126882939257741,0.000237142782693581,1.00000019330558e-05,9.80858023362150]
# m1,m2,a1,a2,L1,I1,I2,g
p_nf_edp=[0.0938439747983407,0.137595970499412,0.108565214610396,0.116779018445656,0.172719203817519,0.000437529429756910,0.00126882939257741,9.80858023362150]

## The following function is the ODE of experimental double pendulum, unlike the simplyfied double pendulum, this double pendulum considers the moment of inertial
function DPE_NoFric!(du,u,p,t)
    # Define the parameters of double pendulum without any friction
    m1,m2,a1,a2,L1,I1,I2,g=p
    
    # The state is : θ1, Θ2, dθ1, dθ2
    du[1]=u[3]
    
    du[2]=u[4]
    
    du[3]=-(L1*a2^3*u[4]^2*m2^2*sin(u[1] - u[2]) + (L1*a2^2*g*m2^2*sin(u[1]))/2 + I2*L1*g*m2*sin(u[1]) + (L1*a2^2*g*m2^2*sin(u[1] - 2*u[2]))/2 + I2*a1*g*m1*sin(u[1]) + (L1^2*a2^2*u[3]^2*m2^2*sin(2*u[1] - 2*u[2]))/2 + I2*L1*a2*u[4]^2*m2*sin(u[1] - u[2]) + a1*a2^2*g*m1*m2*sin(u[1]))/(I1*I2 + L1^2*a2^2*m2^2 + I2*L1^2*m2 + I2*a1^2*m1 + I1*a2^2*m2 - L1^2*a2^2*m2^2*cos(u[1] - u[2])^2 + a1^2*a2^2*m1*m2) 
          #-(L1*a2^3*u[4]^2*m2^2*sin(u[1] - u[2]) + (L1*a2^2*g*m2^2*sin(u[1]))/2 + I2*L1*g*m2*sin(u[1]) + (L1*a2^2*g*m2^2*sin(u[1] - 2*u[2]))/2 + I2*a1*g*m1*sin(u[1]) + (L1^2*a2^2*u[3]^2*m2^2*sin(2*u[1] - 2*u[2]))/2 + I2*L1*a2*u[4]^2*m2*sin(u[1] - u[2]) + a1*a2^2*g*m1*m2*sin(u[1]))/(- L1^2*a2^2*m2^2*cos(u[1] - u[2])^2 + L1^2*a2^2*m2^2 + I2*L1^2*m2 + m1*a1^2*a2^2*m2 + I2*m1*a1^2 + I1*a2^2*m2 + I1*I2)
 
    du[4]=(a2*m2*(2*L1^3*u[3]^2*m2*sin(u[1] - u[2]) - 2*I1*g*sin(u[2]) - 2*L1^2*g*m2*sin(u[2]) + 2*I1*L1*u[3]^2*sin(u[1] - u[2]) - 2*a1^2*g*m1*sin(u[2]) + L1^2*a2*u[4]^2*m2*sin(2*u[1] - 2*u[2]) + 2*L1*a1^2*u[3]^2*m1*sin(u[1] - u[2]) + 2*L1^2*g*m2*cos(u[1] - u[2])*sin(u[1]) + 2*L1*a1*g*m1*cos(u[1] - u[2])*sin(u[1])))/(2*(I1*I2 + L1^2*a2^2*m2^2 + I2*L1^2*m2 + I2*a1^2*m1 + I1*a2^2*m2 - L1^2*a2^2*m2^2*cos(u[1] - u[2])^2 + a1^2*a2^2*m1*m2)) 
          #(a2*m2*(2*L1^3*u[3]^2*m2*sin(u[1] - u[2]) - 2*I1*g*sin(u[2]) - 2*L1^2*g*m2*sin(u[2]) + 2*I1*L1*u[3]^2*sin(u[1] - u[2]) - 2*a1^2*g*m1*sin(u[2]) + L1^2*a2*u[4]^2*m2*sin(2*u[1] - 2*u[2]) + 2*L1*a1^2*u[3]^2*m1*sin(u[1] - u[2]) + 2*L1^2*g*m2*cos(u[1] - u[2])*sin(u[1]) + 2*L1*a1*g*m1*cos(u[1] - u[2])*sin(u[1])))/(2*(- L1^2*a2^2*m2^2*cos(u[1] - u[2])^2 + L1^2*a2^2*m2^2 + I2*L1^2*m2 + m1*a1^2*a2^2*m2 + I2*m1*a1^2 + I1*a2^2*m2 + I1*I2))
 
    return du
end

## The following function is the ODE of the experimental double pendulum with friction
function DPE_Fric!(du,u,p,t)
    # Define the parameters of double pendulum without any friction
    m1,m2,a1,a2,L1,I1,I2,k1,k2,g=p
    
    # The state is : θ1, Θ2, dθ1, dθ2
    du[1]=u[3]
    
    du[2]=u[4]
    
    du[3]=-(I2*u[3]*k1 + I2*u[3]*k2 - I2*u[4]*k2 + a2^2*u[3]*k1*m2 + a2^2*u[3]*k2*m2 - a2^2*u[4]*k2*m2 + L1*a2^3*u[4]^2*m2^2*sin(u[1] - u[2]) + (L1*a2^2*g*m2^2*sin(u[1]))/2 + I2*L1*g*m2*sin(u[1]) + (L1*a2^2*g*m2^2*sin(u[1] - 2*u[2]))/2 + I2*a1*g*m1*sin(u[1]) + (L1^2*a2^2*u[3]^2*m2^2*sin(2*u[1] - 2*u[2]))/2 + L1*a2*u[3]*k2*m2*cos(u[1] - u[2]) - L1*a2*u[4]*k2*m2*cos(u[1] - u[2]) + I2*L1*a2*u[4]^2*m2*sin(u[1] - u[2]) + a1*a2^2*g*m1*m2*sin(u[1]))/(I1*I2 + L1^2*a2^2*m2^2 + I2*L1^2*m2 + I2*a1^2*m1 + I1*a2^2*m2 - L1^2*a2^2*m2^2*cos(u[1] - u[2])^2 + a1^2*a2^2*m1*m2) 
          #-(I2*u[3]*k1 + I2*u[3]*k2 - I2*u[4]*k2 + a2^2*u[3]*k1*m2 + a2^2*u[3]*k2*m2 - a2^2*u[4]*k2*m2 + L1*a2^3*u[4]^2*m2^2*sin(u[1] - u[2]) + (L1*a2^2*g*m2^2*sin(u[1]))/2 + I2*L1*g*m2*sin(u[1]) + (L1*a2^2*g*m2^2*sin(u[1] - 2*u[2]))/2 + I2*a1*g*m1*sin(u[1]) + (L1^2*a2^2*u[3]^2*m2^2*sin(2*u[1] - 2*u[2]))/2 + L1*a2*u[3]*k2*m2*cos(u[1] - u[2]) - L1*a2*u[4]*k2*m2*cos(u[1] - u[2]) + I2*L1*a2*u[4]^2*m2*sin(u[1] - u[2]) + a1*a2^2*g*m1*m2*sin(u[1]))/(- L1^2*a2^2*m2^2*cos(u[1] - u[2])^2 + L1^2*a2^2*m2^2 + I2*L1^2*m2 + m1*a1^2*a2^2*m2 + I2*m1*a1^2 + I1*a2^2*m2 + I1*I2)
 
    du[4]=(I1*u[3]*k2 - I1*u[4]*k2 + L1^2*u[3]*k2*m2 - L1^2*u[4]*k2*m2 + a1^2*u[3]*k2*m1 - a1^2*u[4]*k2*m1 + L1^3*a2*u[3]^2*m2^2*sin(u[1] - u[2]) - L1^2*a2*g*m2^2*sin(u[2]) - I1*a2*g*m2*sin(u[2]) + (L1^2*a2^2*u[4]^2*m2^2*sin(2*u[1] - 2*u[2]))/2 + L1*a2*u[3]*k1*m2*cos(u[1] - u[2]) + L1*a2*u[3]*k2*m2*cos(u[1] - u[2]) - L1*a2*u[4]*k2*m2*cos(u[1] - u[2]) + L1^2*a2*g*m2^2*cos(u[1] - u[2])*sin(u[1]) + I1*L1*a2*u[3]^2*m2*sin(u[1] - u[2]) - a1^2*a2*g*m1*m2*sin(u[2]) + L1*a1^2*a2*u[3]^2*m1*m2*sin(u[1] - u[2]) + L1*a1*a2*g*m1*m2*cos(u[1] - u[2])*sin(u[1]))/(I1*I2 + L1^2*a2^2*m2^2 + I2*L1^2*m2 + I2*a1^2*m1 + I1*a2^2*m2 - L1^2*a2^2*m2^2*cos(u[1] - u[2])^2 + a1^2*a2^2*m1*m2) 
          #(I1*u[3]*k2 - I1*u[4]*k2 + L1^2*u[3]*k2*m2 - L1^2*u[4]*k2*m2 + a1^2*u[3]*k2*m1 - a1^2*u[4]*k2*m1 + L1^3*a2*u[3]^2*m2^2*sin(u[1] - u[2]) - L1^2*a2*g*m2^2*sin(u[2]) - I1*a2*g*m2*sin(u[2]) + (L1^2*a2^2*u[4]^2*m2^2*sin(2*u[1] - 2*u[2]))/2 + L1*a2*u[3]*k1*m2*cos(u[1] - u[2]) + L1*a2*u[3]*k2*m2*cos(u[1] - u[2]) - L1*a2*u[4]*k2*m2*cos(u[1] - u[2]) + L1^2*a2*g*m2^2*cos(u[1] - u[2])*sin(u[1]) + I1*L1*a2*u[3]^2*m2*sin(u[1] - u[2]) - a1^2*a2*g*m1*m2*sin(u[2]) + L1*a1^2*a2*u[3]^2*m1*m2*sin(u[1] - u[2]) + L1*a1*a2*g*m1*m2*cos(u[1] - u[2])*sin(u[1]))/(- L1^2*a2^2*m2^2*cos(u[1] - u[2])^2 + L1^2*a2^2*m2^2 + I2*L1^2*m2 + m1*a1^2*a2^2*m2 + I2*m1*a1^2 + I1*a2^2*m2 + I1*I2)
 
    return du
end

## The following function will calculate the energy of the double pendulum
function HamiltonianEDP(u,p)
    # Define parameters
    m1,m2,a1,a2,L1,I1,I2,g=p
    
    # Extract states
    θ1,θ2,ω1,ω2=u
    
    # Kinetic energy
    T=(I1*ω1^2)/2 + (I2*ω2^2)/2 + (L1^2*ω1^2*m2)/2 + (a1^2*ω1^2*m1)/2 + (a2^2*ω2^2*m2)/2 + L1*a2*ω1*ω2*m2*cos(θ1 - θ2) 

    # Potential energy
    V=-g*m2*(L1*cos(θ1) + a2*cos(θ2)) - a1*g*m1*cos(θ1)
 
    # Total energy
    H=T+V

    return H
end

## The following function will calculate the energy of the double pendulum, input is a matrix, and the output is a vector
function HamiltonianEDP_Vectorized(u,p)
    # Define parameters
    m1,m2,a1,a2,L1,I1,I2,g=p
    
    # Extract states
    θ1=u[1,:]
    θ2=u[2,:]
    ω1=u[3,:]
    ω2=u[4,:]

    # Kinetic energy
    T=(m2*L1^2*ω1.^2)/2 + m2*cos.(θ1.-θ2)*L1*a2.*ω1.*ω2 + (m1*a1^2*ω1.^2)/2 + (m2*a2^2*ω2.^2)/2 + (I1*ω1.^2)/2 + (I2*ω2.^2)/2
    
    # Potential energy
    V=-g*(m2*(L1*cos.(θ1)+a2*cos.(θ2))+a1*m1*cos.(θ1))
 
    # Total energy
    H=T+V
    
    return H
end

## The following function is used to mod the double pendulum's angle so that it stays in -π to π
function ModAngle(solArray)
    if !isempty(solArray)
        for j=1:size(solArray,2)
            if solArray[2,j]>pi
                while solArray[2,j]>pi
                    solArray[2,j]=solArray[2,j]-2*pi
                end
            elseif solArray[2,j]<-pi
                while solArray[2,j]<-pi
                    solArray[2,j]=solArray[2,j]+2*pi
                end
            end
            if solArray[1,j]>pi
                while solArray[1,j]>pi
                    solArray[1,j]=solArray[1,j]-2*pi
                end
            elseif solArray[1,j]<-pi
                while solArray[1,j]<-pi
                    solArray[1,j]=solArray[1,j]+2*pi
                end
            end
        end
    end
    return solArray
end

## This function will calculate the value of dθ2 given values of energy level, θ1, θ2, and dθ1   
function DPE_Cal_dTheta(EngTarget,u,p,angle_index;pos_or_neg=1)
    # Extract the parameters of double pendulum
    m1,m2,a1,a2,L1,I1,I2,g=p
    
    # Calculate the angular velocity
    if angle_index==1
        # Calculate the dtheta_1
        if pos_or_neg==1
            dθ=-((2*EngTarget*I1 - I1*I2*u[4]^2 + 2*EngTarget*L1^2*m2 + 2*EngTarget*a1^2*m1 + 2*L1^3*g*m2^2*cos(u[1]) + 2*a1^3*g*m1^2*cos(u[1]) - (L1^2*a2^2*u[4]^2*m2^2)/2 - I2*L1^2*u[4]^2*m2 - I2*a1^2*u[4]^2*m1 - I1*a2^2*u[4]^2*m2 + 2*L1^2*a2*g*m2^2*cos(u[2]) - a1^2*a2^2*u[4]^2*m1*m2 + 2*I1*L1*g*m2*cos(u[1]) + 2*I1*a1*g*m1*cos(u[1]) + 2*I1*a2*g*m2*cos(u[2]) + (L1^2*a2^2*u[4]^2*m2^2*cos(2*u[1] - 2*u[2]))/2 + 2*L1*a1^2*g*m1*m2*cos(u[1]) + 2*L1^2*a1*g*m1*m2*cos(u[1]) + 2*a1^2*a2*g*m1*m2*cos(u[2]))^(1/2) + L1*a2*u[4]*m2*cos(u[1] - u[2]))/(m2*L1^2 + m1*a1^2 + I1)
        else
            dθ=((2*EngTarget*I1 - I1*I2*u[4]^2 + 2*EngTarget*L1^2*m2 + 2*EngTarget*a1^2*m1 + 2*L1^3*g*m2^2*cos(u[1]) + 2*a1^3*g*m1^2*cos(u[1]) - (L1^2*a2^2*u[4]^2*m2^2)/2 - I2*L1^2*u[4]^2*m2 - I2*a1^2*u[4]^2*m1 - I1*a2^2*u[4]^2*m2 + 2*L1^2*a2*g*m2^2*cos(u[2]) - a1^2*a2^2*u[4]^2*m1*m2 + 2*I1*L1*g*m2*cos(u[1]) + 2*I1*a1*g*m1*cos(u[1]) + 2*I1*a2*g*m2*cos(u[2]) + (L1^2*a2^2*u[4]^2*m2^2*cos(2*u[1] - 2*u[2]))/2 + 2*L1*a1^2*g*m1*m2*cos(u[1]) + 2*L1^2*a1*g*m1*m2*cos(u[1]) + 2*a1^2*a2*g*m1*m2*cos(u[2]))^(1/2) - L1*a2*u[4]*m2*cos(u[1] - u[2]))/(m2*L1^2 + m1*a1^2 + I1)
        end
    else
        # Calculate the dtheta_2
        if pos_or_neg==1
            dθ=((2*EngTarget*I2 - I1*I2*u[3]^2 + 2*EngTarget*a2^2*m2 + 2*a2^3*g*m2^2*cos(u[2]) - (L1^2*a2^2*u[3]^2*m2^2)/2 - I2*L1^2*u[3]^2*m2 - I2*a1^2*u[3]^2*m1 - I1*a2^2*u[3]^2*m2 + 2*L1*a2^2*g*m2^2*cos(u[1]) - a1^2*a2^2*u[3]^2*m1*m2 + 2*I2*L1*g*m2*cos(u[1]) + 2*I2*a1*g*m1*cos(u[1]) + 2*I2*a2*g*m2*cos(u[2]) + (L1^2*a2^2*u[3]^2*m2^2*cos(2*u[1] - 2*u[2]))/2 + 2*a1*a2^2*g*m1*m2*cos(u[1]))^(1/2) - L1*a2*u[3]*m2*cos(u[1] - u[2]))/(m2*a2^2 + I2)
        else
            dθ=-((2*EngTarget*I2 - I1*I2*u[3]^2 + 2*EngTarget*a2^2*m2 + 2*a2^3*g*m2^2*cos(u[2]) - (L1^2*a2^2*u[3]^2*m2^2)/2 - I2*L1^2*u[3]^2*m2 - I2*a1^2*u[3]^2*m1 - I1*a2^2*u[3]^2*m2 + 2*L1*a2^2*g*m2^2*cos(u[1]) - a1^2*a2^2*u[3]^2*m1*m2 + 2*I2*L1*g*m2*cos(u[1]) + 2*I2*a1*g*m1*cos(u[1]) + 2*I2*a2*g*m2*cos(u[2]) + (L1^2*a2^2*u[3]^2*m2^2*cos(2*u[1] - 2*u[2]))/2 + 2*a1*a2^2*g*m1*m2*cos(u[1]))^(1/2) + L1*a2*u[3]*m2*cos(u[1] - u[2]))/(m2*a2^2 + I2)
        end
    end
            
    return dθ
end

## The following function will calculate the loss value of the P.O. 
function DPElossVal_WithEng(dtheta_opt,EngTarget,p,cb,dt,T,tspan,index;angle_index=1,pos_or_neg=1)
    x0g=[]
    dtheta2=0
    if index==1
        dtheta2=DPE_Cal_dTheta(EngTarget,[0,pi,dtheta_opt,0],p,angle_index;pos_or_neg=pos_or_neg)
        x0g=[0,pi,dtheta_opt,dtheta2]
    else
        dtheta2=DPE_Cal_dTheta(EngTarget,[pi,0,dtheta_opt,0],p,angle_index;pos_or_neg=pos_or_neg)
        x0g=[pi,0,dtheta_opt,dtheta2]
    end

    DPE_ProbOpt=ODEProblem(DPE_NoFric!,x0g,tspan,p)

    sol=Array(DifferentialEquations.solve(DPE_ProbOpt,Vern9(),dt=dt,tstops=T,
        abstol=1e-16,reltol=1e-16,callback=cb,save_everystep = false,save_start=false,save_end=true))

    if index==1
        loss=norm(sol[3]+dtheta_opt)+norm(pi-sol[2])
    else
        loss=norm(sol[3]+dtheta_opt)+norm(pi-sol[1])
    end

    return loss
end

## Define a optimization function
function DPEOptPO_WithEng(dtheta_opt,N,tol,EngTarget,p,cb,dt,T,tspan,index,opt;angle_index=1,pos_or_neg=1)
    pin=1
    lossVal(dtheta_opt)=DPElossVal_WithEng(dtheta_opt,EngTarget,p,cb,dt,T,tspan,index;angle_index=angle_index,pos_or_neg=pos_or_neg)
    
    while lossVal(dtheta_opt)>tol && pin<N
        dvy= ForwardDiff.derivative(lossVal,dtheta_opt)
        val=lossVal(dtheta_opt)
        dtheta_opt=dtheta_opt-1e-3*dvy
        #Flux.Optimise.update!(opt, dtheta_opt, dvy)
        println("The loss value is :$(val) at iteration $(pin)...")
        pin=pin+1
    end

    dtheta2=0
    if index==1
        dtheta2=DPE_Cal_dTheta(EngTarget,[0,pi,dtheta_opt,0],p,angle_index;pos_or_neg=pos_or_neg)
        x0g=[0,pi,dtheta_opt,dtheta2]
    else
        dtheta2=DPE_Cal_dTheta(EngTarget,[pi,0,dtheta_opt,0],p,angle_index;pos_or_neg=pos_or_neg)
        x0g=[pi,0,dtheta_opt,dtheta2]
    end
    
    return x0g
end

## The following function will calculate the loss value of the P.O. 
function DPElossVal(dtheta_opt,EngTarget,p,cb,dt,T,tspan,index)
    x0g=[]
    if index==1
        x0g=[0,pi,dtheta_opt[1],dtheta_opt[2]]
    else
        x0g=[pi,0,dtheta_opt[1],dtheta_opt[2]]
    end

    DPE_ProbOpt=ODEProblem(DPE_NoFric!,x0g,tspan,p)

    sol=Array(DifferentialEquations.solve(DPE_ProbOpt,Vern9(),dt=dt,tstops=T,abstol=1e-16,reltol=1e-16,callback=cb,save_everystep = false,save_start=false,save_end=true))

    if index==1
        loss=sum((sol[3:4]+dtheta_opt).^2)+(EngTarget-HamiltonianEDP(x0g,p))^2#+norm(pi-sol[2])#+(EngTarget-HamiltonianEDP(x0g,p))^2
    else
        loss=sum((sol[3:4]+dtheta_opt).^2)+10*(EngTarget-HamiltonianEDP(x0g,p))^2#+norm(pi-sol[1])#+(EngTarget-HamiltonianEDP(x0g,p))^2
    end

    return loss
end

## Define an optimization function for calculating the PO
function DPEOptPO(dtheta_opt,N,tol,EngTarget,p,cb,dt,T,tspan,index,opt)
    pin=1
    lossVal(dtheta_opt)=DPElossVal(dtheta_opt,EngTarget,p,cb,dt,T,tspan,index)
    
    while lossVal(dtheta_opt)>tol && pin<N
        dvy= ForwardDiff.gradient(lossVal,dtheta_opt)
        val=lossVal(dtheta_opt)
        Flux.Optimise.update!(opt, dtheta_opt, dvy)
        println("The loss value is :$(val) at iteration $(pin)...")
        pin=pin+1
    end

    if index==1
        return x0g=[0,pi,dtheta_opt[1],dtheta_opt[2]]
    else
        return x0g=[pi,0,dtheta_opt[1],dtheta_opt[2]]
    end
    
    return x0g
end

## The following function will calculate the loss value of the P.O. (for the 2π rotational cases, and experimental double pendulum)
function DPElossValTwoPi(dtheta_opt,EngTarget,p,cb,dt,T,tspan,index)
    x0g=[]
    if index==1
        x0g=[0,pi,dtheta_opt[1],dtheta_opt[2]]
    else
        x0g=[pi,0,dtheta_opt[1],dtheta_opt[2]]
    end

    DP_ProbOpt=ODEProblem(DPE_NoFric!,x0g,tspan,p)

    sol=Array(DifferentialEquations.solve(DP_ProbOpt,Vern9(),dt=dt,tstops=T,
        abstol=1e-16,reltol=1e-16,callback=cb,save_everystep = false,save_start=false,save_end=true))

    if index==1
        return loss=sum(abs.(sol[3:4]-dtheta_opt))+abs(sol[1])+100*abs(EngTarget-HamiltonianEDP(x0g,p))
    elseif index==2
        return loss=abs(sol[1]-pi)+sum(abs.(sol[3:4]-dtheta_opt))+abs(EngTarget-HamiltonianEDP(x0g,p))
    elseif index==3
        return loss=sum(abs.(sol[3:4]-dtheta_opt))+100*abs(EngTarget-HamiltonianEDP(x0g,p))
    end
end

## Define a optimization function (for the 2π rotational cases, and experimental double pendulum)
function DPEOptPOTwoPi(dtheta_opt,N,tol,EngTarget,p,cb,dt,T,tspan,index,δ;opt=[])
    # (The learning rate needs to be adjusted manully, start with 1e-2 and then decrease the learningrate by 
    # two orders of magnitude every 500 epochs)
    pin=1
    lossVal(dtheta_opt)=DPElossValTwoPi(dtheta_opt,EngTarget,p,cb,dt,T,tspan,index)
    
    while lossVal(dtheta_opt)>tol && pin<N
        dvy= ForwardDiff.gradient(lossVal,dtheta_opt)
        val=lossVal(dtheta_opt)
        if opt==[]
            dtheta_opt=dtheta_opt-δ*dvy
        else
            Flux.Optimise.update!(opt, dtheta_opt, dvy)
        end
        println("The loss value is :$(val) at iteration $(pin)...")
        pin=pin+1
    end

    if index==1
        return x0g=[0,pi,dtheta_opt[1],dtheta_opt[2]]
    else
        return x0g=[pi,0,dtheta_opt[1],dtheta_opt[2]]
    end

end

## The following function will calculate the loss value of the special orbit that connect the two saddle points
function DPElossValSpecial(dtheta_opt,EngTarget,p,cb,dt,T,tspan)
    x0g=[]

    x0g=[0,pi,dtheta_opt[1],dtheta_opt[2]]

    DP_ProbOpt=ODEProblem(DPE_NoFric!,x0g,tspan,p)

    sol=Array(DifferentialEquations.solve(DP_ProbOpt,Vern9(),dt=dt,tstops=T,
        abstol=1e-16,reltol=1e-16,callback=cb,save_everystep = false,save_start=false,save_end=true))

    return loss=abs(sol[1]-pi)+100*abs(EngTarget-HamiltonianEDP(x0g,p))

end

## Define a optimization function for calculate the special orbit that connect the two saddle point
function DPEOptSpecial(dtheta_opt,N,tol,EngTarget,p,cb,dt,T,tspan,index1,δ;opt=[])
    # (The learning rate needs to be adjusted manully, start with 1e-2 and then decrease the learningrate by 
    # two orders of magnitude every 500 epochs)
    pin=1
    lossVal(dtheta_opt)=DPElossValSpecial(dtheta_opt,EngTarget,p,cb,dt,T,tspan)
    
    while lossVal(dtheta_opt)>tol && pin<N
        dvy= ForwardDiff.gradient(lossVal,dtheta_opt)
        val=lossVal(dtheta_opt)
        if opt==[]
            dtheta_opt=dtheta_opt-δ*dvy
        else
            Flux.Optimise.update!(opt, dtheta_opt, dvy)
        end
        println("The loss value is :$(val) at iteration $(pin)...")
        pin=pin+1
    end

    return x0g=[0,pi,dtheta_opt[1],dtheta_opt[2]]
end

## The following function will calculate the loss value of the P.O.
function DPElossValHomo(dtheta_opt,EngTarget,p,cb,dt,T,tspan,index)
    x0g=[]

    x0g=[0,0,dtheta_opt[1],dtheta_opt[2]]

    DPE_ProbOpt=ODEProblem(DPE_NoFric!,x0g,tspan,p)

    sol=Array(DifferentialEquations.solve(DPE_ProbOpt,Vern9(),dt=dt,tstops=T,abstol=1e-16,reltol=1e-16,callback=cb,save_everystep = false,save_start=false,save_end=true))

    #loss=abs(sol[1])+sum(abs.(sol[3:4]-dtheta_opt))+abs(EngTarget-HamiltonianEDP(x0g,p))
    loss=(sol[1])^2+sum((sol[3:4]-dtheta_opt).^2)+(EngTarget-HamiltonianEDP(x0g,p))^2

    return loss
end

## Define a optimization function (for the 2π rotational cases)
function DPEOptHomo(dtheta_opt,N,tol,EngTarget,p,cb,dt,T,tspan,index,δ;opt=[])
    # (The learning rate needs to be adjusted manully, start with 1e-2 and then decrease the learningrate by 
    # two orders of magnitude every 500 epochs)
    pin=1
    lossVal(dtheta_opt)=DPElossValHomo(dtheta_opt,EngTarget,p,cb,dt,T,tspan,index)
    
    while lossVal(dtheta_opt)>tol && pin<N
        dvy= ForwardDiff.gradient(lossVal,dtheta_opt)
        val=lossVal(dtheta_opt)
        if opt==[]
            dtheta_opt=dtheta_opt-δ*dvy
        else
            Flux.Optimise.update!(opt, dtheta_opt, dvy)
        end
        println("The loss value is :$(val) at iteration $(pin)...")
        pin=pin+1
    end

    x0g=[0,0,dtheta_opt[1],dtheta_opt[2]]

    return x0g
end

## The following function will calculate the Monodromy matrix of double pendulum
function CalSymMonoE(nx,p) 
    # First create a symbolic variable and calculate the Jacobian matrix
    @ModelingToolkit.variables x_var[1:nx]
    @ModelingToolkit.variables p_var[1:size(p)[1]]
    # Next, get the expression of the function RHS
    t_var=0
    RHS=DPE_NoFric!(similar(x_var),x_var,p_var,t_var)

    # Then calculate the jacobian
    DRHS=ModelingToolkit.jacobian(RHS,x_var)

    # Transfer the calculated result into function
    write("MonoMatrixDPE.jl",string(build_function(DRHS,x_var,p_var)[1]))

    # Include those function
    MonoMatrixDPE=include("MonoMatrixDPE.jl")
    
    return MonoMatrixDPE
end

## The following function will define the monodromy matrix
function MonoODEDPE!(du,u,p,t;nx=4)
    # The first four states are the normal ode terms
    x_n=u[1:nx]
    du[1:nx]=DPE_NoFric!(similar(x_n),x_n,p,t)
    
    # The remaining terms comes from the augmented states
    x_m=reshape(u[nx+1:end],nx,nx)
    dMono=MonoMatrixDPE(x_n,p)*x_m
    
    # Now assign the new derivative
    du[nx+1:end]=reshape(dMono,nx*nx,1)
    
    return du
end

## The following function will simulate the PO using monodromy matrix
function GetValMonoE(Xsaddle,nx,cb,tspan,p_nf;dt=0.0001)
    # Get the initial condition
    X0=vcat(Xsaddle,Int.(reshape(Matrix(I,nx,nx),nx*nx,1)[:]))
    # Get the Mono ode problem
    MonoProb=ODEProblem(MonoODEDPE!,X0,tspan,p_nf)
    # Solve the ode
    sol_Period_saddle=DifferentialEquations.solve(MonoProb,Vern9(),dt=dt,tstops=T,abstol=1e-16,reltol=1e-16,
    callback=cb,saveat=dt,save_everystep = true,save_start=true,save_end=true)
    # Get the array
    Time_Period=Array(sol_Period_saddle.t)
    Data_Period=Array(sol_Period_saddle)
    # Get the period data
    POData=Data_Period[1:nx,1:1:end]
    # After we simulated for one period, let's compute the mono matrix at time T
    MonoT=reshape(Data_Period[nx+1:end,end],nx,nx)
    # Then calculate the eigenvalue, we should have λ1=λ2=1, and λ3*λ4=1
    eigMono=eigvals(MonoT)
    # Then calculate the eigenvector
    eigVecMono=real(eigvecs(MonoT))
    
    return Data_Period,POData,eigMono,eigVecMono,MonoT
end

## Calculate the Monodromy matrix of the frictional Model
function CalSymMonoE_Fric(nx,p) 
    # First create a symbolic variable and calculate the Jacobian matrix
    @ModelingToolkit.variables x_var[1:nx]
    @ModelingToolkit.variables p_var[1:size(p)[1]]
    # Next, get the expression of the function RHS
    t_var=0
    RHS=DPE_Fric!(similar(x_var),x_var,p_var,t_var)

    # Then calculate the jacobian
    DRHS=ModelingToolkit.jacobian(RHS,x_var)

    # Transfer the calculated result into function
    write("MonoMatrixDPE_Fric.jl",string(build_function(DRHS,x_var,p_var)[1]))

    # Include those function
    MonoMatrixDPE_Fric=include("MonoMatrixDPE_Fric.jl")
    
    return MonoMatrixDPE_Fric
end

function MonoODEDPE_Fric!(du,u,p,t;nx=4)
    # The first four states are the normal ode terms
    x_n=u[1:nx]
    du[1:nx]=DPE_Fric!(similar(x_n),x_n,p,t)
    
    # The remaining terms comes from the augmented states
    x_m=reshape(u[nx+1:end],nx,nx)
    dMono=MonoMatrixDPE_Fric(x_n,p)*x_m
    
    # Now assign the new derivative
    du[nx+1:end]=reshape(dMono,nx*nx,1)
    
    return du
end

function GetValMonoE_Fric(Xsaddle,nx,cb,tspan,p_nf;dt=0.0001)
    # Get the initial condition
    X0=vcat(Xsaddle,Int.(reshape(Matrix(I,nx,nx),nx*nx,1)[:]))
    # Get the Mono ode problem
    MonoProb=ODEProblem(MonoODEDPE_Fric!,X0,tspan,p_nf)
    # Solve the ode
    sol_Period_saddle=DifferentialEquations.solve(MonoProb,Vern9(),dt=dt,tstops=T,abstol=1e-16,reltol=1e-16,
    callback=cb,saveat=dt,save_everystep = true,save_start=true,save_end=true)
    # Get the array
    Time_Period=Array(sol_Period_saddle.t)
    Data_Period=Array(sol_Period_saddle)
    # Get the period data
    POData=Data_Period[1:nx,1:1:end]
    # After we simulated for one period, let's compute the mono matrix at time T
    MonoT=reshape(Data_Period[nx+1:end,end],nx,nx)
    # Then calculate the eigenvalue, we should have λ1=λ2=1, and λ3*λ4=1
    eigMono=eigvals(MonoT)
    # Then calculate the eigenvector
    eigVecMono=real(eigvecs(MonoT))
    
    return Data_Period,POData,eigMono,eigVecMono,MonoT
end

## The following function will calculate the normalized stable and unstable eigenvector
function GetSatbleVector(MonoT,index1,index2)
    # Stable direction
    vecS=real(eigvecs(MonoT)[:,index1])
    vecS=vecS/norm(vecS)
    # Unstable direction
    vecU=real(eigvecs(MonoT)[:,index2])
    vecU=vecU/norm(vecU)
    
    return vecS,vecU
end

## The following function will calculate the generate initial condition given stable and unstable vector
function GetPeriodMono(PeriodicTraj,δV,positive,Stable,vecU,vecS,gap)
    k1=1
    Inc=PeriodicTraj
    
    for i=1:size(PeriodicTraj,2)
        if Stable==1
            Inc[:,i]=Inc[:,i]+positive*δV*vecS
        else
            Inc[:,i]=Inc[:,i]+positive*δV*vecU
        end
    end
    
    return Inc[:,1:gap:end]
end

## This function will be used to simulate a bounch of initial condition and gather tube information
function SimulateManifold(X0,par,backward,Tsim,everystepsave,cb;saveend=true,saveat=true,experimental=false)
    # Now generate an EnsembleProblem to simulate all the initial conditions 
    # Define ODE problems
    dt=0.00001
    tspan=(0.0,Tsim)
    
    if backward==1 
        dt=-0.00001
        tspan=(0.0,-Tsim)
    end
    
    DP_Mono_Prob=[]

    if size(par,1)>5 && experimental==false
        if backward==1 
            par[6:7]=abs.(par[6:7])
        else
            par[6:7]=abs.(par[6:7])
        end
        if experimental==false
            DP_Mono_Prob=ODEProblem(DP_Fric!,X0[:,1],tspan,par)
        else
            DP_Mono_Prob=ODEProblem(DPE_Fric!,X0[:,1],tspan,par)
        end
    else
        if experimental==false
            DP_Mono_Prob=ODEProblem(DP_NoFric!,X0[:,1],tspan,par)
        else
            DP_Mono_Prob=ODEProblem(DPE_NoFric!,X0[:,1],tspan,par)
        end
    end

    # Prepare for ensemble problem
    function prob_func(prob,i,repeat)
      remake(prob,u0=X0[:,i])
    end

    # Define ensemble problem
    ensemble_prob_nf = EnsembleProblem(DP_Mono_Prob,prob_func=prob_func)

    # Now calcultae and plot the poincare map
    println("\t Calculating...")
    NC=size(X0,2)
    if saveat==true
        # After the package update, somehow the EnsembleThreads() will creash the VS Code
        @time PoincareSim_nf = DifferentialEquations.solve(ensemble_prob_nf,
            Vern9(),EnsembleSerial(),trajectories=NC,dt=dt,callback=cb,
            save_everystep = everystepsave, 
            saveat=0.001,
            save_start=true,save_end=saveend,
            abstol=1e-15,reltol=1e-15);

        return PoincareSim_nf
    else
        @time PoincareSim_nf = DifferentialEquations.solve(ensemble_prob_nf,
            Vern9(),EnsembleSerial(),trajectories=NC,dt=dt,callback=cb,
            save_everystep = everystepsave, 
            save_start=true,save_end=saveend,
            abstol=1e-15,reltol=1e-15);

        return PoincareSim_nf
    end
end

## This function will be used to generate the animation of double pendulum
function DP_Animate(DpPlot,ax1,Data,p,gap,color1,color2,color3,cmap,mks1,mks2,ls1,Ratio,path,framerate,xlims1,xlims2,ylims1,ylims2)
    # Get the parameter of double pendulum
    M1,M2,L1,L2,g=p
    
    # Get the size of matrix
    Data=hcat(Data[:,1:gap:end-1],Data[:,end])
    sz1=size(Data,2)

    # Define the matrix
    X=zeros(5,sz1)

    # Hide the ticks
    hidespines!(ax1)
    hidedecorations!(ax1)
    
    # Set the plot limits
    limits!(ax1,xlims1,xlims2,ylims1,ylims2)

    # Run for loop
    sz1=0
    record(DpPlot,path,framerate=framerate) do io
        for i=1:size(Data,2)
            tstamp=Ratio*i/size(Data,2)
            # Get the location of point mass 1
            x1pos=L1*sin(Data[1,i])
            y1pos=-L1*cos(Data[1,i])
            # Get the location of point mass 2
            x2pos=x1pos+L2*sin(Data[2,i])
            y2pos=y1pos-L2*cos(Data[2,i])
            # Push the data
            sz1=sz1+1
            X[:,sz1]=[tstamp,x1pos,y1pos,x2pos,y2pos]
            
            # Bs and Be
            Bs=-0.1*maximum([xlims1,xlims2])
            Be=0.1*maximum([ylims1,ylims2])
            N=12
            dB=2*Be/N

            Makie.lines!(ax1,[Bs,Be],[0,0],color=(:black,1),linewidth=10,show_axis=false) 

            for j=1:N
                Makie.lines!(ax1,[Bs+(j-1)*dB,Bs+(j)*dB],[0,dB],color=(:black,1),linewidth=10,show_axis=false)
            end
            
            # Plot the arm
            Makie.lines!(ax1,[0,x1pos],[0,y1pos],color=:black,linewidth=ls1+8,show_axis=false)
            Makie.lines!(ax1,[0,x1pos],[0,y1pos],color=color1,linewidth=ls1,show_axis=false)
            
            Makie.lines!(ax1,[x1pos,x2pos],[y1pos,y2pos],color=:black,linewidth=ls1+8) 
            Makie.lines!(ax1,[x1pos,x2pos],[y1pos,y2pos],color=color2,linewidth=ls1) 
            
            # Plot the base
            Makie.scatter!(ax1,[0],[0],color=:aquamarine2,markersize=mks1,strokewidth=4)
            Makie.scatter!(ax1,[0],[0],color=:black,markersize=28,strokewidth=2,marker="x")
            
            # Plot the mass
            Makie.scatter!(ax1,[x1pos],[y1pos],color=color3,markersize=mks2,strokewidth=4)
            Makie.scatter!(ax1,[x2pos],[y2pos],color=color3,markersize=mks2,strokewidth=4)
            
            #DpPlot.center = false
            
            # Record frame
            Makie.recordframe!(io)
            
            # Empty current plots
            empty!(ax1)

            current_figure()
        end
    end
    
    return nothing#DpPlot,streams
end

## This function will be used to generate the animation of double pendulum on the cart
function DP_OnCart_Animate(DpPlot,ax1,Data,p,gap,color1,color2,color3,cmap,mks1,mks2,ls1,Ratio,path,framerate,xlims1,xlims2,ylims1,ylims2)
    # Get the parameter of double pendulum
    M1,M2,L1,L2,g=p
    
    # Get the size of matrix
    Data=hcat(Data[:,1:gap:end-1],Data[:,end])
    sz1=size(Data,2)
    print(size(Data))
    # Define the matrix
    X=zeros(5,sz1)
    
    # Hide the ticks
    hidespines!(ax1)
    hidedecorations!(ax1)

    # Set the plot limits
    limits!(ax1,xlims1,xlims2,ylims1,ylims2)
    
    # Run for loop
    sz1=0
    record(DpPlot,path,framerate=framerate) do io
        for i=1:size(Data,2)
            tstamp=Ratio*i/size(Data,2)
            # Get the location of point mass 1
            xbpos=Data[3,i]

            # Get the location of point mass 1
            x1pos=L1*sin(Data[1,i])+xbpos
            y1pos=-L1*cos(Data[1,i])

            # Get the location of point mass 2
            x2pos=x1pos+L2*sin(Data[2,i])
            y2pos=y1pos-L2*cos(Data[2,i])

            # Push the data
            sz1=sz1+1
            X[:,sz1]=[tstamp,x1pos,y1pos,x2pos,y2pos]
            
            # Bs and Be
            Bs=xlims1
            Be=xlims2
            N=100
            dB=2*Be/N
            
            Makie.lines!(ax1,[Bs,Be],[0,0],linewidth=5,show_axis=false,color=:gray) 
            Makie.lines!(ax1,[0,0],[0,ylims1],linewidth=5,linestyle=:dash,color=:gray) 

            # for j=1:N
            #     DpPlot=Makie.lines!([Bs+(j-1)*dB,Bs+(j)*dB],[0,dB],color=(:black,1),linewidth=5,show_axis=false)
            # end
            
            # Plot the arm
            Makie.poly!(ax1,Point2f0[(xbpos-0.05, -0.02), (xbpos+0.05, -0.02), (xbpos+0.05, 0.02), (xbpos-0.05,0.02)], color = :aquamarine2,strokecolor = :black, strokewidth = 4)
            
            Makie.lines!(ax1,[xbpos,x1pos],[0,y1pos],color=:black,linewidth=ls1+8,show_axis=false)
            Makie.lines!(ax1,[xbpos,x1pos],[0,y1pos],color=color1,linewidth=ls1,show_axis=false)
            
            Makie.lines!(ax1,[x1pos,x2pos],[y1pos,y2pos],color=:black,linewidth=ls1+8) 
            Makie.lines!(ax1,[x1pos,x2pos],[y1pos,y2pos],color=color2,linewidth=ls1) 
            
            # Plot the base
            Makie.scatter!(ax1,[xbpos],[0],color=:aquamarine2,markersize=mks1,strokewidth=4)
            Makie.scatter!(ax1,[xbpos],[0],color=:black,markersize=28,strokewidth=2,marker="x")
            
            # Plot the mass
            Makie.scatter!(ax1,[x1pos],[y1pos],color=color3,markersize=mks2,strokewidth=4)
            Makie.scatter!(ax1,[x2pos],[y2pos],color=color3,markersize=mks2,strokewidth=4)
            
            # Set the plot limits
            limits!(ax1,xlims1,xlims2,ylims1,ylims2)
            
            #DpPlot.center = false
            
            # Record frame
            Makie.recordframe!(io)
            
            # Empty current plots
            empty!(ax1)

            current_figure()
        end
    end 
    
    return nothing#DpPlot,streams
end

## This plot will be used to plot the PO using Makie
function PlotPeriod(ax1,Data,p_nf,Transform,AxisNum1,AxisNum2,AxisNum3,ThreeD)
    if Transform==1
        Data=CoordinateTransform(Data,p_nf)
    end
    
    if ThreeD==1
        Makie.lines!(ax1,Data[AxisNum1,:],Data[AxisNum2,:],Data[AxisNum3,:])
    else
        Makie.scatter!(ax1,[Data[AxisNum1,1]],[Data[AxisNum2,1]],marker='*',color=:orange,markersize=14)
        Makie.lines!(ax1,Data[AxisNum1,:],Data[AxisNum2,:])
    end

    return ax1
end


## This function will be used to plot the tube structure
function MakiePlotMono(MonoPlot,PeriodicTraj,Traj,color,ThreeD,AxisNum1,AxisNum2,AxisNum3,Mod,Transform,mks,lineplot,lw;gap=1)
    TrajData=[]
    if lineplot==0
        for i=1:gap:size(Traj,3)
            if i==1
                TrajData=Array(Traj[i])
            else
                TrajData=hcat(TrajData,Array(Traj[i]))
            end
        end

        if Mod==1
            TrajData=ModAngle(TrajData)
        end

        if Transform==1
            TrajData=CoordinateTransform(TrajData,p_nf)
            PeriodicTraj=CoordinateTransform(PeriodicTraj,p_nf)
        end

        if ThreeD==1
            Makie.lines!(MonoPlot,PeriodicTraj[AxisNum1,:],PeriodicTraj[AxisNum2,:],PeriodicTraj[AxisNum3,:],
                linewidth=lw)
            Makie.scatter!(MonoPlot,TrajData[AxisNum1,:],TrajData[AxisNum2,:],TrajData[AxisNum3,:],
                color=color,markersize=mks,strokewidth=0)
            Makie.xlabel!("theta_1")
            Makie.ylabel!("theta_2")
            if AxisNum3==3
                Makie.zlabel!("dtheta_1")
            elseif AxisNum3==4
                Makie.zlabel!("dtheta_2")
            end
        else
            Makie.scatter!(MonoPlot,[PeriodicTraj[AxisNum1,1]],[PeriodicTraj[AxisNum2,1]],marker='*',color=:orange,markersize=mks)
            Makie.lines!(MonoPlot,PeriodicTraj[AxisNum1,:],PeriodicTraj[AxisNum2,:],linewidth=lw)
            Makie.scatter!(MonoPlot,TrajData[AxisNum1,:],TrajData[AxisNum2,:],color=color,markersize=mks)
            Makie.xlabel!("theta_1")
            Makie.ylabel!("theta_2")
        end
    else
        if Transform==1
            PeriodicTraj=CoordinateTransform(PeriodicTraj,p_nf)
        end
        
        if ThreeD==1
            Makie.lines!(MonoPlot,PeriodicTraj[AxisNum1,:],PeriodicTraj[AxisNum2,:],PeriodicTraj[AxisNum3,:],
                linewidth=lw)
        else
            Makie.scatter!(MonoPlot,[PeriodicTraj[AxisNum1,1]],[PeriodicTraj[AxisNum2,1]],marker='*',color=:orange,markersize=mks)
            Makie.lines!(MonoPlot,PeriodicTraj[AxisNum1,:],PeriodicTraj[AxisNum2,:])
        end
        
        for i=1:gap:size(Traj,3)
            TrajData=Array(Traj[i])
            
            if Mod==1
                TrajData=ModAngle(TrajData)
            end

            if Transform==1
                TrajData=CoordinateTransform(TrajData,p_nf)
                PeriodicTraj=CoordinateTransform(PeriodicTraj,p_nf)
            end
            
            if ThreeD==1
                Makie.lines!(MonoPlot,TrajData[AxisNum1,:],TrajData[AxisNum2,:],TrajData[AxisNum3,:],
                    color=color,linewidth=lw,strokewidth=0)
                Makie.xlabel!("theta_1")
                Makie.ylabel!("theta_2")
            else
                Makie.lines!(MonoPlot,TrajData[AxisNum1,:],TrajData[AxisNum2,:],color=color,linewidth=lw)
                Makie.xlabel!("theta_1")
                Makie.ylabel!("theta_2")
            end
        end

        if AxisNum3==3
            Makie.zlabel!(MonoPlot,"dtheta_1")
        elseif AxisNum3==4
            Makie.zlabel!(MonoPlot,"dtheta_2")
        end
    end
    Makie.AxisAspect(1)
    return MonoPlot
end

## This function achieves similar objective as the function above, but using Plots.jl
function PlotsPlotMono(MonoPlot,PeriodicTraj,Traj,color,ThreeD,AxisNum1,AxisNum2,AxisNum3,Mod,Transform,mks,lineplot,lw,lw2,camview,gap_point;gap_line=1)
    TrajData=[]
    if lineplot==0
        for i=1:gap_line:size(Traj,3)
            if i==1
                TrajData=Array(Traj[i])
            else
                TrajData=hcat(TrajData,Array(Traj[i]))
            end
        end

        if Mod==1
            TrajData=ModAngle(TrajData)
        end

        if Transform==1
            TrajData=CoordinateTransform(TrajData,p_nf)
            PeriodicTraj=CoordinateTransform(PeriodicTraj,p_nf)
        end

        if ThreeD==1
            MonoPlot=Plots.plots!(MonoPlot,PeriodicTraj[AxisNum1,:],PeriodicTraj[AxisNum2,:],
                PeriodicTraj[AxisNum3,:],linewidth=lw,color="black",label="")
            MonoPlot=Plots.scatter!(TrajData[AxisNum1,:],TrajData[AxisNum2,:],TrajData[AxisNum3,:],
                color=color,markersize=mks,camera=camview)
        else
            MonoPlot=Plots.scatter!([PeriodicTraj[AxisNum1,1]],[PeriodicTraj[AxisNum2,1]],
                color=:blue,markersize=mks)
            MonoPlot=Plots.plot!(PeriodicTraj[AxisNum1,:],PeriodicTraj[AxisNum2,:],color="black")
            MonoPlot=Plots.scatter!(TrajData[AxisNum1,:],TrajData[AxisNum2,:],color=color,markersize=mks)
        end
    else
        PeriodicTraj=hcat(PeriodicTraj[:,1:gap_point:end],PeriodicTraj[:,end])

        if Transform==1
            PeriodicTraj=CoordinateTransform(PeriodicTraj,p_nf)
        end
        
        if ThreeD==1
            MonoPlot=Plots.plot!(MonoPlot,PeriodicTraj[AxisNum1,:],PeriodicTraj[AxisNum2,:],PeriodicTraj[AxisNum3,:],
                linewidth=lw2,color="black",label="")
        else
            MonoPlot=Plots.scatter!([PeriodicTraj[AxisNum1,1]],[PeriodicTraj[AxisNum2,1]],label="",
                marker = (:circle, 6, 0.99, :blue, stroke(2, 0.8, :black, :solid)),framestyle=:box)
            MonoPlot=Plots.plot!(PeriodicTraj[AxisNum1,:],PeriodicTraj[AxisNum2,:],
                linewidth=lw2,color="black",label="")
        end        
        
        for i=1:gap_line:size(Traj,3)
            TrajData=Array(Traj[i])
            
            TrajData=hcat(TrajData[:,1:gap_point:end],TrajData[:,end])

            if Mod==1
                TrajData=ModAngle(TrajData)
            end

            if Transform==1
                TrajData=CoordinateTransform(TrajData,p_nf)
                PeriodicTraj=CoordinateTransform(PeriodicTraj,p_nf)
            end
            
            if ThreeD==1
                MonoPlot=Plots.plot!(TrajData[AxisNum1,:],TrajData[AxisNum2,:],TrajData[AxisNum3,:],
                    color=color,linewidth=lw,label="",camera=camview)
            else
                MonoPlot=Plots.plot!(TrajData[AxisNum1,:],TrajData[AxisNum2,:],
                    color=color,linewidth=lw,label="")
            end
            
        end
        
    end
    
    return MonoPlot
end

## This function will be used to plot the Poincare cut
function MakiePlotFinalPoint(PoincarPlot,Traj,color,Axis;PlotTheta2=0,LinePlot=1)
    n=size(Traj,3)
    
    Data=zeros(4,n)
    j=0
    for i=1:n
        Data[:,i]=Array(Traj[i])[:,end]
        j=j+1
    end
    
    if LinePlot==0
        if PlotTheta2==0
            if Axis==1
                Makie.scatter!(PoincarPlot,Data[1,:],Data[3,:],color=color,markersize=20,strokewidth=0)
            else
                Makie.scatter!(PoincarPlot,Data[2,:],Data[4,:],color=color,markersize=20,strokewidth=0)
            end
        else
            if Axis==1
                Makie.scatter!(PoincarPlot,Data[1,:],Data[4,:],color=color,markersize=20,strokewidth=0)
            else
                Makie.scatter!(PoincarPlot,Data[2,:],Data[3,:],color=color,markersize=20,strokewidth=0)
            end
        end
    else
        if PlotTheta2==0
            if Axis==1
                Makie.lines!(PoincarPlot,Data[1,:],Data[3,:],color=color,markersize=20,strokewidth=0)
            else
                Makie.lines!(PoincarPlot,Data[2,:],Data[4,:],color=color,markersize=20,strokewidth=0)
            end
        else
            if Axis==1
                Makie.lines!(PoincarPlot,Data[1,:],Data[4,:],color=color,markersize=20,strokewidth=0)
            else
                Makie.lines!(PoincarPlot,Data[2,:],Data[3,:],color=color,markersize=20,strokewidth=0)
            end
        end
    end

    return Data,PoincarPlot
end

## This function will be used to plot the Poincare cut using Plots.jl
function PlotFinalPointPlots(PoincarPlot,Traj,color,Axis;PlotTheta2=0,LinePlot=1)
    n=size(Traj,3)
    
    Data=zeros(4,n)
    j=0
    for i=1:n
        Data[:,i]=Array(Traj[i])[:,end]
        j=j+1
    end
    if LinePlot==1
        if PlotTheta2==0
            if Axis==1
                PoincarPlot=Plots.plot!(PoincarPlot,Data[1,:],Data[3,:],color=color,
                    label="",linewidth=2.5,markerstrokewidth=0,markerstrokealpha=0,framestyle=:box)
            else
                PoincarPlot=Plots.plot!(PoincarPlot,Data[2,:],Data[4,:],color=color,
                    label="",linewidth=2.5,markerstrokewidth=0,markerstrokealpha=0,framestyle=:box)
            end
        else
            if Axis==1
                PoincarPlot=Plots.plot!(PoincarPlot,Data[1,:],Data[4,:],color=color,
                    label="",linewidth=2.5,markerstrokewidth=0,markerstrokealpha=0,framestyle=:box)
            else
                PoincarPlot=Plots.plot!(PoincarPlot,Data[2,:],Data[3,:],color=color,
                    label="",linewidth=2.5,markerstrokewidth=0,markerstrokealpha=0,framestyle=:box)
            end
        end
    else
        if PlotTheta2==0
            if Axis==1
                PoincarPlot=Plots.scatter!(PoincarPlot,Data[1,:],Data[3,:],color=color,
                    label="",markersize=1.5,markerstrokewidth=0,markerstrokealpha=0,framestyle=:box)
            else
                PoincarPlot=Plots.scatter!(PoincarPlot,Data[2,:],Data[4,:],color=color,
                    label="",markersize=1.5,markerstrokewidth=0,markerstrokealpha=0,framestyle=:box)
            end
        else
            if Axis==1
                PoincarPlot=Plots.scatter!(PoincarPlot,Data[1,:],Data[4,:],color=color,
                    label="",markersize=1.5,markerstrokewidth=0,markerstrokealpha=0,framestyle=:box)
            else
                PoincarPlot=Plots.scatter!(PoincarPlot,Data[2,:],Data[3,:],color=color,
                    label="",markersize=1.5,markerstrokewidth=0,markerstrokealpha=0,framestyle=:box)
            end
        end
    end
 
    return Data,PoincarPlot
end

## Define a function for plotting the Poincare cut
function MakiePlotPoincare(PoincareSim_nf,color,mks,AxisNum1,AxisNum2,AxisNum3)
    # Get the plotting handle
    PoincarePlot=Scene()

    # Pre-define the data
    DataPoincare=[]
    try
        DataPoincare=Array(PoincareSim_nf[1])
    catch

    end

    # Store the data
    for i=2:size(PoincareSim_nf,3)
        try
            # Get the data
            dummy=Array(PoincareSim_nf[i])
            # Store it
            DataPoincare=hcat(DataPoincare,dummy)
        catch
            
        end
    end

    if AxisNum3==0
        PoincarePlot=Makie.scatter!(PoincarePlot,DataPoincare[AxisNum1,:],DataPoincare[AxisNum2,:],
            color=color,markersize=mks,strokewidth=0)
    else
        PoincarePlot=Makie.scatter!(PoincarePlot,DataPoincare[AxisNum1,:],DataPoincare[AxisNum2,:],DataPoincare[AxisNum3,:],
            color=color,markersize=mks,strokewidth=0)
    end
    
    if AxisNum1==1
        Makie.xlabel!("θ1")
    elseif AxisNum1==2
        Makie.xlabel!("θ2")
    elseif AxisNum1==3
        Makie.xlabel!("dθ1")
    elseif AxisNum1==4
        Makie.xlabel!("dθ2")
    end

    if AxisNum2==1
        Makie.ylabel!("θ1")
    elseif AxisNum2==2
        Makie.ylabel!("θ2")
    elseif AxisNum2==3
        Makie.ylabel!("dθ1")
    elseif AxisNum2==4
        Makie.ylabel!("dθ2")
    end

    if AxisNum3==1
        Makie.zlabel!("θ1")
    elseif AxisNum3==2
        Makie.zlabel!("θ2")
    elseif AxisNum3==3
        Makie.zlabel!("dθ1")
    elseif AxisNum3==4
        Makie.zlabel!("dθ2")
    end

    return DataPoincare,PoincarePlot
end


## This function will extract the poincare cut from the simulation data
function ExtractIncFromEndPoint(Traj)
    N=size(Traj,3)

    Inc=zeros(4,N)

    for i=1:N
        Data=Array(Traj[i])

        Inc[:,i]=Data[:,end]
    end
    return Inc
end

## This function will generate random initial condition near the saddle point
function GetIncRand(x0,θup,θdown,dθup,dθdown,Np,ShowPlot,saddle)
    Inc=[]

    # Get the increment
    δθ=(θup-θdown)/Np
    δdθ=(dθup-dθdown)/Np

    # Get the possible valus of the initial conditions
    θArray=Array(θdown:δθ:θup)
    dθArray=Array(dθdown:δdθ:dθup)

    # Arrange them
    Inc=zeros(4,size(θArray,1)*size(dθArray,1))
    IncPlot=Makie.Scene()
    if saddle==1
        Inc[2,:].=pi
        Inc[4,:].=x0[4,1]
        # Run for loop
        k1=1
        for i in θArray
            for j in dθArray
                Inc[1,k1]=i+x0[1,1]
                Inc[3,k1]=j+x0[3,1]
                k1=k1+1
            end
        end
        if ShowPlot==1
            IncPlot=Makie.scatter!(Inc[1,:],Inc[2,:],Inc[3,:],markersize=8,strokewidth=0)
            Makie.xlabel!("theta_1")
            Makie.ylabel!("theta_2")
            Makie.zlabel!("dtheta_1")
            display(IncPlot)
        end
    else
        Inc[1,:].=pi
        Inc[3,:].=x0[3,1]
        # Run for loop
        k1=1
        for i in θArray
            for j in dθArray
                Inc[2,k1]=i
                Inc[4,k1]=j+x0[4,1]
                k1=k1+1
            end
        end
        if ShowPlot==1
            IncPlot=Makie.scatter!(Inc[1,:],Inc[2,:],Inc[4,:],markersize=8,strokewidth=0)
            Makie.xlabel!("theta_1")
            Makie.ylabel!("theta_2")
            Makie.zlabel!("dtheta_2")
            display(IncPlot)
        end
    end
    
    return Inc
end 

## This function is sued to create a color wheel
function GetColorWheel(Points,P_m,σx,σy;HalfAngle=false)
    IterNum=0
    if size(Points,1)==4
        AngleArray=zeros(size(Points,2),1)
        RatioArray=zeros(size(Points,2),1)
        IterNum=size(Points,2)
    else
        AngleArray=zeros(size(Points,1),1)
        RatioArray=zeros(size(Points,1),1)
        IterNum=size(Points,1)
    end

    for i=1:IterNum
        if size(Points,1)==4
            P_now=Points[:,i]
            x=P_now[1]-P_m[1]
            y=P_now[3]-P_m[2]
        else
            P_now=Points[i]
            x=P_now[1]-P_m[1]
            y=P_now[2]-P_m[2]
        end

        r_p=sqrt(x^2+y^2)
        angle_now=atand(y,x)
        if HalfAngle==true
            angle_now=angle_now*2
        end
        r_now=(σx*σy)/sqrt(σx^2*sind(angle_now)^2+σy^2*cosd(angle_now)^2)

        radius_ratio=r_p/r_now

        AngleArray[i,1]=angle_now
        RatioArray[i,1]=radius_ratio
    end

    ColorWheel=[HSV(AngleArray[i],RatioArray[i],1) for i in 1:IterNum]

    return ColorWheel
end





















