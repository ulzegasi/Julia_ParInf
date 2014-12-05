##
##
## ============================================================================================
##
##

function V_fast(theta,q)                                # Fast dynamics (harmonic and data terms)
    
    # Variables:
    # ---------------------------
    
    beta    = theta[1]
    tau     = theta[2]
    
    K       = T*exp(beta)
    gamma   = exp(tau)

    # Harmonic term:
    # ---------------------------
    
    tmp=0.
    tmp=q[1]^2+q[N]^2-2.0*q[1]*q[2]
    out=0.
    
    # Instead of
    #
    # for s = 2:N-1
    #     tmp += 2.0*(q[s]^2-q[s]*q[s+1])
    # end
    #
    # we use the (faster!) nested for cycles:

    for s=1:n
        tmp1 = 0.
        for k = 2:j
            tmp1 += 2.0*(q[(s-1)*j+k]^2-q[(s-1)*j+k]*q[(s-1)*j+k+1])
        end
        if s == 1 
            tmp += tmp1
        else 
            tmp += 2.0*( q[(s-1)*j+1]^2 - q[(s-1)*j+1]*q[(s-1)*j+2]) + tmp1
        end        
    end

    out = K*tmp/(2.0*gamma*dt)

    # Data points:
    # ---------------------------

    for s=0:n
        out += ( y[s+1] - r[s*j+1]*exp(q[s*j+1]) )^2/(2*sigma^2)
    end

    return out
end

##
##
## ============================================================================================
##
##

# Calculate derivatives of V_fast w.r.t. theta:
# ---------------------------

V_fast_theta      = function(theta) V_fast(theta,q) end
V_fast_der_theta  = forwarddiff_gradient(V_fast_theta, Float64, fadtype=:dual, n=s)

# Equivalent to:
#
# function V_fast_theta(theta) 
#    V_fast(theta,u) 
# end
#
# V_fast_der_theta = forwarddiff_gradient(V_fast_theta, Float64, fadtype=:dual, n=s)

##
##
## ============================================================================================
##
##

function V_fast_der(theta,q)                            # Derivatives of V_fast w.r.t. u

    # Variables:
    # ---------------------------

    beta     = theta[1]
    tau      = theta[2]
    
    K        = T*exp(beta)
    gamma    = exp(tau)

    out = zeros(N)

    # Harmonic term
    # ---------------------------

    # First "end point"
    # ---------------------------

    out[1] = K*(q[1]-q[2])/(gamma*dt)
    
    # Intermediate "end-point" components:
    # ---------------------------
    
    for s = 2:(N-1)
        out[s] = K*(2*q[s]-q[s-1]-q[s+1])/(gamma*dt)
    end

    # Last end point:
    # ---------------------------
    
    out[N] = K*(q[N]-q[N-1])/(gamma*dt)

    # Data (there is an exponent q[s*j+1] = u[s*j+1]):
    # ---------------------------
    
    for i=0:n
        out[i*j+1] -= r[i*j+1] * exp(q[i*j+1]) * 
                    ( y[i+1] - r[i*j+1]*exp(q[i*j+1]) )/sigma^2
    end

    return [V_fast_der_theta(theta), out]
end

##
##
## ============================================================================================
##
##

function V_slow(theta,q) 

    # Variables:
    # ---------------------------

    beta  = theta[1]
    tau   = theta[2]
    
    K     = T*exp(beta)
    gamma = exp(tau)
    

    # Definition of rho:
    # ---------------------------
    
    rho = Array(Float64,N)
    rho = K * lnr_der + gamma + 1

    # Potential, boundary terms:
    # ---------------------------
    out = 0.0
    out = ( exp(- q[N]) +  q[N]*rho[N-1] - exp(- q[1]) - q[1]*rho[1] )/gamma 
    
    # Potential, bulk terms:
    # ---------------------------
    
    for i=2:(N-1)
        out +=  dt * ( rho[i] - exp( -q[i] ) )^2     / (2*gamma*K) -
                dt * q[i] * ( (rho[i]-rho[i-1])/dt ) / gamma
    end
    
    out -= T * (2+gamma) / (4*K) + (N/2) * (beta-tau)
    
    return out
end

##
##
## ============================================================================================
##
##

# Calculate derivatives of V_slow w.r.t. theta:
# ---------------------------

V_slow_theta      = function(theta) V_slow(theta,q) end 
V_slow_der_theta  = forwarddiff_gradient(V_slow_theta, Float64, fadtype=:dual, n=s)

# Equivalent to:
# function V_slow_theta(theta) 
#     V_slow(theta,u) 
# end
# 
# V_slow_der_theta = forwarddiff_gradient(V_slow_theta, Float64, fadtype=:dual, n=s)

##
##
## ============================================================================================
##
##

function V_slow_der(theta,q)

   # Variables:
   # ---------------------------

    beta    = theta[1]
    tau     = theta[2]
    
    K       = T*exp(beta)
    gamma   =   exp(tau)
    
    
    # Definition of rho:
    # ---------------------------
    
    rho = Array(Float64,N)
    rho = K * lnr_der + gamma + 1

    # Calculation of partial derivatives of V_slow w.r.t. q:   
    # --------------------------- 
    
    out = Array(Float64,N)
    
    # First component:
    # ---------------------------
    
    out[1] = ( exp(-q[1]) - rho[1] )/gamma


    # Intermediate components:
    # ---------------------------
    
    for i=2:(N-1)
        out[i] = (dt/(gamma*K)) * ( 
                                ( rho[i] - exp(-q[i])) * 
                                exp(-q[i]) -
                                K*(rho[i]-rho[i-1])/dt  
                                )
    end

    # Last component:
    # ---------------------------

    out[N] = -( exp(-q[N]) - rho[N-1] )/gamma 
     
   
    return [V_slow_der_theta(theta),out]

end

##
##
## ============================================================================================
##
##

function RESPA(theta, q, p, mp, dtau, nn)               # (Tuckerman et al., JCP 97(3), 1992 )
    
    # Define small time step:
    # ---------------------------

    deltatau = dtau/nn 

    # First long-range step (dtau):
    # ---------------------------

    force = -V_slow_der(theta,q)
    for i=1:(N+s) 
        p[i] = p[i]  + (dtau/2)* force[i]
    end

    # Short-range steps (nn*deltatau):
    # ---------------------------

    for counter=1:nn                                    # Verlet integrator 
        force_old = -V_fast_der(theta,q)                # (Tuckerman et al., JCP 97(3), 1992, eq. 2.17)
        for i=1:s
            theta[i] = theta[i] + deltatau * ( p[i]  + (deltatau/2) * force_old[i] )/mp[i]
        end
        for i=1:N
            q[i] = q[i]  + deltatau * ( p[i+s]  + (deltatau/2) * force_old[i+s] )/mp[i+s]
        end
        force_new = -V_fast_der(theta,q)
        for i=1:(N+s) 
            p[i] = p[i]  + (deltatau/2)*( force_old[i] + force_new[i] )
        end
    end

    # Long-range step:
    # ---------------------------

    force = -V_slow_der(theta,q)
    for i=1:(N+s) 
        p[i] = p[i]  + (dtau/2)* force[i]
    end

end
