
function V_fast(theta,u)                                # Fast dynamics (harmonic and data terms)
    
    # Variables:
    # ---------------------------
    
    beta    = theta[1]
    tau     = theta[2]
    
    K       = T*exp(beta)
    gamma   = exp(tau)

    # Harmonic term:
    # ---------------------------
    
    tmp = 0.
    for s=1:n
        tmp1 = 0.
        for k=2:j
            tmp1 += k*u[(s-1)*j+k]^2/( (k-1)*dt )
        end
        tmp += ( u[(s-1)*j+1] - u[s*j+1] )^2/( j*dt ) + tmp1
    end
    
    out = K * tmp / ( 2*gamma )
    
    # Data:
    # ---------------------------

    for s=0:n
        out += ( y[s+1] - r[s*j+1]*exp(u[s*j+1]) )^2/(2*sigma^2)
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
V_fast_theta = function(theta) V_fast(theta,u) end
V_fast_der_theta  = forwarddiff_gradient(V_fast_theta, Float64, fadtype=:dual, n=s)
##
##
## ============================================================================================
##
##

function V_fast_der(theta,u)                            # Derivatives of V_fast w.r.t. {u}

    # Variables:
    # ---------------------------

    beta     = theta[1]
    tau      = theta[2]
    
    K        = T*exp(beta)
    gamma    = exp(tau)

    out = zeros(N)

    # Harmonic term
    # ---------------------------

    # First "end-point" component
    # ---------------------------
    
    out[1] = K*( u[1] - u[j+1] )/( gamma*j*dt )

    # Intermediate "end-point" components:
    # ---------------------------
    
    for s=1:(n-1)
        out[s*j+1] = ( K/(gamma*j*dt) ) * 
                        ( 2*u[s*j+1] - u[(s+1)*j+1] - u[(s-1)*j+1] )

    end
    
    # Staging components:
    # ---------------------------
    
    for s=0:(n-1)
       for k=2:j
           out[s*j+k] = (K*k/(gamma*dt*(k-1))) * u[s*j+k]
       end
    end
    
    # Last "end-point" component:
    # ---------------------------
    
    out[N] = - (K/(gamma*j*dt)) * ( u[(n-1)*j+1] - u[N] )
    
    # Data (there is an exponent q[s*j+1] = u[s*j+1]):
    # ---------------------------
    
    for s=0:n
        out[s*j+1] -= r[s*j+1] * exp(u[s*j+1]) * 
                    ( y[s+1] - r[s*j+1]*exp(u[s*j+1]) )/sigma^2
    end

    ########################### Derivative w.r.t. theta ######################
    # V_fast_der_theta_2 = Array(Float64,s)
	#
    # tmp = 0.
    # for s=1:n
    #     tmp1 = 0.
    #     for k=2:j
    #         tmp1 += k*u[(s-1)*j+k]^2/( (k-1)*dt )
    #     end
    #     tmp += ( u[(s-1)*j+1] - u[s*j+1] )^2/( j*dt ) + tmp1
    # end
	# 
    # V_fast_der_theta_2[1] = K * tmp / ( 2*gamma )
    # V_fast_der_theta_2[2] = - K * tmp / ( 2*gamma )
    ########################################################################################    

    # return [V_fast_der_theta_2, out]
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
    
    out = ( exp(- q[N]) +  q[N]*rho[N-1] - exp(- q[1]) - q[1]*rho[1] )/gamma 
    
    # Potential, bulk terms:
    # ---------------------------
    
    for i=2:(N-1)
        out +=  dt * ( rho[i] - exp( -q[i] ) )^2     / (2.0*gamma*K) -
                dt * q[i] * ( (rho[i]-rho[i-1])/dt ) / gamma
    end
    
    out -= T * (2.0+gamma) / (4.0*K) + (N/2.0) * (beta-tau) -
            (0.5*beta^2+beta*(log(T)-5.0)) - (tau^2+tau)
    
    return out
end

##
##
## ============================================================================================
##
##
# Calculate derivatives of V_slow w.r.t. theta:
# --------------------------- 
V_slow_theta = function(theta) V_slow(theta,q) end
V_slow_der_theta  = forwarddiff_gradient(V_slow_theta, Float64, fadtype=:dual, n=s)
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
     
    ########################### Derivative w.r.t. theta ######################
    
    # V_slow_der_theta_2 = Array(Float64,s)
	# 
    # tmp1=0.
    # for i=2:(N-1)
    #     tmp1 += 1/(2*K*gamma*exp(2*q[i])) * (-dt * ((gamma+1)*exp(q[i])-1)^2 - 
    #         K^2*exp(2*q[i])*(-dt*lnr_der[i]^2 +2*q[i]*(lnr_der[i]-lnr_der[i-1]) ) ) 
    # end
	# 
    # tmp1 += -K/gamma*(q[1]*lnr_der[1] - q[N]*lnr_der[N-1]) + T/(4*K)*(gamma+2) - N/2
	# 
    # tmp2=0.
    # for i=2:(N-1)
    #     tmp2 += 1/(2*K*gamma*exp(2*q[i])) * ( dt * (gamma^2*exp(2*q[i])-(exp(q[i])-1)^2) + 
    #         K*exp(q[i])*(-2*dt*(exp(q[i])-1)*lnr_der[i] + 
    #             K*exp(q[i])*(-dt*lnr_der[i]^2+2*q[i]*(lnr_der[i]-lnr_der[i-1])) ) )
    # end
	# 
    # tmp2 += 1/gamma * ( exp(-q[1])-exp(-q[N]) + q[1] +
    #  K*lnr_der[1]*q[1]-q[N]*(1+K*lnr_der[N-1]) ) - (1/4)*(T*gamma/K-2*N) 
	# 
    # V_slow_der_theta_2[1] = tmp1
    # V_slow_der_theta_2[2] = tmp2

    ########################################################################################


    # Variable transformations (q -> u): 
    # (Tuckerman et al., JCP 99 (1993), eqs. 2.22-23)
    # ---------------------------
    
    out_u = zeros(N)
    
    # First component (s=0, k=1):
    # ---------------------------
   
   tmp = 0.
   for l=2:j
       tmp += (j-l+1)*out[l]/j
   end
   out_u[1] = out[1] + tmp
    
   # Other "end-point" components (k=1): 
   # ---------------------------
    
    for s=1:(n-1)
        tmp = 0.
        for l=2:j
            tmp += (j-l+1)*out[s*j+l]/j + 
                    (l-1)*out[(s-1)*j+l]/j
        end
        out_u[s*j+1] = out[s*j+1] + tmp
    end
    
    # Staging components:
    # ---------------------------
    
    for s=0:(n-1)
        for k=2:j
            out_u[s*j+k] = out[s*j+k] + (k-2)*out_u[s*j+k-1]/(k-1)
        end
    end
    
    # Last component:
    # ---------------------------
    
    tmp = 0.
    for l=2:j
        tmp +=  (l-1)*out[(n-1)*j+l]/j
    end
    out_u[N] = out[N] + tmp    

    return [V_slow_der_theta(theta),out_u]
    # return [V_slow_der_theta_2, out_u]

end

##
##
## ============================================================================================
##
##

function RESPA(theta, u, q, p, mp, dtau, nn)               # (Tuckerman et al., JCP 97(3), 1992 )
    
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
    
    for counter_fast_respa = 1:nn                       # Verlet integrator 
        
        force_old = -V_fast_der(theta,u)                # (Tuckerman et al., JCP 97(3), 1992, eq. 2.17)
        for i=1:s
            theta[i] = theta[i] + deltatau * ( p[i]  + (deltatau/2) * force_old[i] )/mp[i]
        end
        for i=1:N
            u[i] = u[i]  + deltatau * ( p[i+s]  + (deltatau/2) * force_old[i+s] )/mp[i+s]
        end
        force_new = -V_fast_der(theta,u)
        for i=1:(N+s) 
            p[i] = p[i]  + (deltatau/2)*( force_old[i] + force_new[i] )
        end
    end
    
    # Back transformations (u -> q) to update the {q} variables :
    # (Tuckerman et al., JCP 99 (1993), eq. 2.19)
    # ---------------------------
    ss=n-1
    while ss>=0
        k=j+1
        q[ss*j+k] = u[ss*j+k]
        while k>2
           k -= 1
           q[ss*j+k] = u[ss*j+k] + u[ss*j+1]/k + (k-1)*q[ss*j+k+1]/k
        end
        ss -= 1
    end
    q[1] = u[1]


    # Long-range step:
    # ---------------------------

    force = -V_slow_der(theta,q)
    for i=1:(N+s) 
        p[i] = p[i]  + (dtau/2)* force[i]
    end

end