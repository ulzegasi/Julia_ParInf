##########################################################################
## Fast potential as a function
##########################################################################

function V_fast_fun(theta,u)    # Fast dynamics (harmonic and data terms)
    
    # Variables:
    # ---------------------------
    
    out = 0.0

    beta    = theta[1]
    tau     = theta[2]
    
    K       = T*exp(beta)
    gamma   = exp(tau)

    out = ( y[1] - r[1]*exp(u[1]) )^2/(2*sigma^2)
    for i=1:n
        out += K * ( u[(i-1)*j+1] - u[i*j+1] )^2/( j*dt )   / ( 2*gamma ) +
                    ( y[i+1] - r[i*j+1]*exp(u[i*j+1]) )^2/(2*sigma^2)
        for k=2:j
            out += k*u[(i-1)*j+k]^2/( (k-1)*dt ) * K/(2*gamma)
        end
    end

    return out
end

##
##
# Calculate derivatives of V_fast w.r.t. the parameters theta:
# --------------------------- 
# V_fast_fun_theta      = function(theta) V_fast_fun(theta,u) end
# V_fast_fun_der_theta  = forwarddiff_gradient(V_fast_fun_theta, Float64, fadtype=:dual, n=s)

# V_fast_fun_u      = function(u) V_fast_fun(theta,u) end
# V_fast_fun_der_u  = forwarddiff_gradient(V_fast_fun_u, Float64, fadtype=:dual, n=N)
##
##

##########################################################################
## Fast potential as an expression
##########################################################################

V_fast = quote       # Fast dynamics (harmonic and data terms)

    # Variables:
    # ---------------------------
    
    bet      = theta[1]
    tau      = theta[2]
    
    K        = T*exp(bet)
    gamm     = exp(tau)

    # Harmonic term:
    # ---------------------------
    out_fast = 0.0
    out_temp = ( y[1] - r[1]*exp(u[1]) )^2/(2*sigma^2)
    for i=1:n
        out_temp += K * ( u[(i-1)*j+1] - u[i*j+1] )^2/( j*dt )   / ( 2*gamm ) +
                    ( y[i+1] - r[i*j+1]*exp(u[i*j+1]) )^2/(2*sigma^2)
    end
    for i=1:n
        for k=2:j
            out_temp += k*u[(i-1)*j+k]^2/( (k-1)*dt )*K/(2*gamm)
        end
    end
    
    out_fast = out_temp
    
end                         

##
##
# Calculate derivatives of V_fast w.r.t. the parameters theta and the coordinates {u}
# --------------------------- 
V_fast_der_theta = rdiff(V_fast, outsym=:out_fast, theta = ones(Float64, s) ) 
V_fast_der_u     = rdiff(V_fast, outsym=:out_fast, u = ones(Float64, N) ) 
@eval V_fast_der(theta,u) = [($V_fast_der_theta)[2][1],($V_fast_der_u)[2][1]]
##
##

##########################################################################
## Slow potential as a function
##########################################################################

function V_slow_fun(theta,u) 

    # Variables:
    # ---------------------------

    beta  = theta[1]
    tau   = theta[2]
    
    K     = T*exp(beta)
    gamma = exp(tau)
    

    # Definition of rho:
    # ---------------------------
    
    # rho = Array(Float64,N)
    # rho = K * lnr_der + gamma + 1

    # Potential, boundary terms:
    # ---------------------------
    
    out = ( exp(- u[N]) + u[N]*(K * lnr_der[N-1] + gamma + 1) -
            exp(- u[1]) - u[1]*(K * lnr_der[1]   + gamma + 1) ) / gamma

    
    # Potential, bulk terms:
    # ---------------------------

    for i = 1:(n-1)
        out +=  dt * ( (K * lnr_der[i*j+1] + gamma + 1) - exp( -u[i*j+1] ) )^2 / (2.0*gamma*K) -
                        dt * u[i*j+1] * ( K*(lnr_der[i*j+1]-lnr_der[i*j])/dt ) / gamma
    end

    for i = 1:n
        for k = 2:j
            tmp = (j-k+1)*u[(i-1)*j+1]/j
            for l = k : (j+1)
                tmp += (k-1)*u[(i-1)*j+l]/(l-1)
            end
            out += dt * ( (K * lnr_der[(i-1)*j+k] + gamma + 1) - exp( -tmp ) )^2 / (2.0*gamma*K) -
                            dt * tmp * ( K*(lnr_der[(i-1)*j+k]-lnr_der[(i-1)*j+k-1])/dt ) / gamma
        end
    end

    # for i=2:(N-1)
    #     out +=  dt * ( rho[i] - exp( -q[i] ) )^2     / (2.0*gamma*K) -
    #             dt * q[i] * ( (rho[i]-rho[i-1])/dt ) / gamma
    # end
    
    out -= T * (2.0+gamma) / (4.0*K) + (N/2.0) * (beta-tau) -
            (0.5*beta^2+beta*(log(T)-5.0)) - (tau^2+tau)
    
    return out
end

##
##
# Calculate derivatives of V_slow w.r.t. the parameters theta:
# --------------------------- 
# V_slow_fun_theta      = function(theta) V_slow_fun(theta,u) end
# V_slow_fun_der_theta  = forwarddiff_gradient(V_slow_fun_theta, Float64, fadtype=:dual, n=s)

# V_slow_fun_u      = function(u) V_slow_fun(theta,u) end
# V_slow_fun_der_u  = forwarddiff_gradient(V_slow_fun_u, Float64, fadtype=:dual, n=N)
##
##

##########################################################################
## Slow potential as an expression
##########################################################################

V_slow = quote

    # Variables:
    # ---------------------------

    bet   = theta[1]
    tau   = theta[2]
    
    K     = T*exp(bet)
    gamm  = exp(tau)


    # Potential, boundary terms:
    # ---------------------------
    
    out_slow = ( exp(- u[N]) + u[N]*(K * lnr_der[N-1] + gamm + 1) -
                 exp(- u[1]) - u[1]*(K * lnr_der[1]   + gamm + 1) )/gamm
    
    # Potential, bulk terms:
    # ---------------------------
    
    for i = 1:(n-1)
        out_slow +=  dt * ( (K * lnr_der[i*j+1] + gamm + 1) - exp( -u[i*j+1] ) )^2 / (2.0*gamm*K) -
                        dt * u[i*j+1] * ( K*(lnr_der[i*j+1]-lnr_der[i*j])/dt ) / gamm
    end

    for i = 1:n
        for k = 2:j
            tmp = (j-k+1)*u[(i-1)*j+1]/j
            for l = k : (j+1)
                tmp += (k-1)*u[(i-1)*j+l]/(l-1)
            end
            out_slow += dt * ( (K * lnr_der[(i-1)*j+k] + gamm + 1) - exp( -tmp ) )^2 / (2.0*gamm*K) -
                            dt * tmp * ( K*(lnr_der[(i-1)*j+k]-lnr_der[(i-1)*j+k-1])/dt ) / gamm
        end
    end

    # for i = 1:n
    #     for k = 2:j
    #         out_slow += dt * ( (K * lnr_der[(i-1)*j+k] + gamma + 1) - exp( -q[(i-1)*j+k] ) )^2 / (2.0*gamma*K) -
    #                         dt * q[(i-1)*j+k] * ( K*(lnr_der[(i-1)*j+k]-lnr_der[(i-1)*j+k-1])/dt ) / gamma
    #     end
    # end

    out_slow -= ( T * (2.0+gamm) / (4.0*K) + (N/2.0) * (bet-tau) -
            (0.5*bet^2+bet*(log(T)-5.0)) - (tau^2+tau) )

end

##
##
# Calculate derivatives of V_slow w.r.t. the parameters theta and the coordinates {q}
# --------------------------- 
V_slow_der_theta = rdiff(V_slow, outsym=:out_slow, theta = ones(Float64, s) ) 
V_slow_der_u     = rdiff(V_slow, outsym=:out_slow, u = ones(Float64, N) ) 
@eval V_slow_der(theta,u) = [($V_slow_der_theta)[2][1],($V_slow_der_u)[2][1]]
##
##

##########################################################################
## RESPA
##########################################################################

function RESPA(theta, u, p, mp, dtau, nn, deltatau)     # (Tuckerman et al., JCP 97(3), 1992 )
    

    # First long-range step (dtau):
    # ---------------------------
    
    force = -V_slow_der(theta,u)
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
    # ss=n-1
    # while ss>=0
    #     k=j+1
    #     q[ss*j+k] = u[ss*j+k]
    #     while k>2
    #        k -= 1
    #        q[ss*j+k] = u[ss*j+k] + u[ss*j+1]/k + (k-1)*q[ss*j+k+1]/k
    #     end
    #     ss -= 1
    # end
    # q[1] = u[1]


    # Long-range step:
    # ---------------------------
    force = -V_slow_der(theta,u)
    for i=1:(N+s) 
        p[i] = p[i]  + (dtau/2)* force[i]
    end

end