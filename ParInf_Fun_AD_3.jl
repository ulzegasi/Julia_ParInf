
function V_fast(theta,u)    # Fast dynamics (harmonic and data terms)
    
    # Variables:
    # ---------------------------
    
    beta    = theta[1]
    tau     = theta[2]
    
    K       = T*exp(beta)
    gamma   = exp(tau)

    ###############################################################################
    ###############################################################################
    ################### OLD SERIAL CODE ############################# 
    ###Harmonic term:
    ### ---------------------------
    ##
    ##tmp = 0.
    ##for s=1:n
    ##    tmp1 = 0.
    ##    for k=2:j
    ##        tmp1 += k*u[(s-1)*j+k]^2/( (k-1)*dt )
    ##    end
    ##    tmp += ( u[(s-1)*j+1] - u[s*j+1] )^2/(  j*dt ) + tmp1
    ##end
    ##
    ##out = K * tmp / ( 2*gamma )
    ##
    ### Data:
    ### ---------------------------
    ##
    ##for s=0:n
    ##    out += ( y[s+1] - r[s*j+1]*exp(u[s*j+1]) )^2/(2*sigma^2)
    ##end
    ##
    ###############################################################################
    ###############################################################################
    
    ## Parallel part:
    tmp = 0.0
    tmp = 
        sum(map(fetch, 
            {@spawnat wp ( uind = [localindexes(u)[1]]; loctmp = 0.0; 
                for jwp = 1:size(localpart(u),1)
                    if !in(uind[jwp],ty)
                        k = uind[jwp]-ty[findfirst(x->x>uind[jwp],ty)-1]+1;
                        loctmp += (K/(2*gamma))*k*(localpart(u)[jwp])^2 / ((k-1)*dt)
                    elseif in(uind[jwp],ty)
                        loctmp += ( y[((uind[jwp]-1)/j)+1]-r[uind[jwp]]*exp(localpart(u)[jwp]) )^2/(2*sigma^2)
                    end
                end;
                loctmp ) for wp in workers()}))

    ## On local process (1):
    tmp1 = 0.0
    for i = 1:n
        tmp1 += (K/(2*gamma))*( u[(i-1)*j+1] - u[i*j+1] )^2 / (  j*dt )
    end

    out = tmp + tmp1
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

function V_fast_der(theta,u,der_out)    # Derivatives of V_fast w.r.t. {u}

    # Variables:
    # ---------------------------

    beta     = theta[1]
    tau      = theta[2]
    
    K        = T*exp(beta)
    gamma    = exp(tau)

    for wp = workers()                          # parallel on worker processes 
        @spawnat wp ( uind = [localindexes(u)[1]];
            for jwp = 1:size(localpart(der_out),1) 
                if !in(uind[jwp],ty)
                    k = uind[jwp]-ty[findfirst(x->x>uind[jwp],ty)-1]+1;
                    localpart(der_out)[jwp] = (K*k/(gamma*dt*(k-1))) * localpart(u)[jwp]
                elseif in(uind[jwp],ty[2:end-1])
                    localpart(der_out)[jwp] = ( K/(gamma*j*dt) ) *
                    ( 2*localpart(u)[jwp]-u[ty[findin(ty,uind[jwp])-1][1]] - u[ty[findin(ty,uind[jwp])+1][1]] )-
                            r[uind[jwp]] * exp(localpart(u)[jwp])*( y[((uind[jwp]-1)/j)+1] -
                                r[uind[jwp]]*exp(localpart(u)[jwp]) )/sigma^2
                elseif uind[jwp] == 1
                    localpart(der_out)[jwp] = K*( localpart(u)[jwp] - u[ty[2]] )/( gamma*j*dt )-
                        r[1] * exp(localpart(u)[jwp]) * ( y[1] - r[1]*exp(localpart(u)[jwp]) )/sigma^2
                elseif uind[jwp] == N 
                    localpart(der_out)[jwp] = K*( localpart(u)[jwp] - u[ty[end-1]] )/( gamma*j*dt )-
                        r[N] * exp(localpart(u)[jwp]) * ( y[n+1] - r[N]*exp(localpart(u)[jwp]) )/sigma^2
                end
            end
        )
    end
    ###################################################################
    ############### OLD SERIAL ########################################
    ## # Harmonic term
    ## # ---------------------------
    ##
    ## # First "end-point" component
    ## # ---------------------------
    ##
    ## der_out[1] = K*( u[1] - u[j+1] )/( gamma*j*dt )
    ##
    ## # Intermediate "end-point" components:
    ## # ---------------------------
    ##
    ## for s=1:(n-1)
    ##     out[s*j+1] = ( K/(gamma*j*dt) ) * 
    ##                     ( 2*u[s*j+1] - u[(s+1)*j+1] - u[(s-1)*j+1] )
    ## end
    ##
    ## # Staging components:
    ## # ---------------------------
    ##
    ## for s=0:(n-1)
    ##    for k=2:j
    ##        out[s*j+k] = (K*k/(gamma*dt*(k-1))) * u[s*j+k]
    ##    end
    ## end
    ##
    ## # Last "end-point" component:
    ## # ---------------------------
    ##
    ## der_out[N] = - (K/(gamma*j*dt)) * ( u[(n-1)*j+1] - u[N] )
    ##
    ## # Data (there is an exponent q[s*j+1] = u[s*j+1]):
    ## # ---------------------------
    ##
    ## for s=0:n
    ##     out[s*j+1] -= r[s*j+1] * exp(u[s*j+1]) * 
    ##                 ( y[s+1] - r[s*j+1]*exp(u[s*j+1]) )/sigma^2
    ## end
    ###################################################################
    ###################################################################

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
    return (V_fast_der_theta(theta), der_out)

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

    tmp  = 0.0
    tmp1 = 0.0

    # Potential, bulk terms:
    # ---------------------------
    ######################################################################
    ## Possible parallel version
    ##
    ## tmp = 
    ##     sum(map(fetch, 
    ##         {@spawnat wp ( uind = [localindexes(q)[1]]; loctmp = 0.0; 
    ##             for jwp = 1:size(localpart(q),1)
    ##                 if uind[jwp] != 1 && uind[jwp] != N
    ##                     loctmp += dt*( rho[uind[jwp]]-exp(-localpart(q)[jwp]) )^2 / (2.0*gamma*K) -
    ##                             localpart(q)[jwp]* (rho[uind[jwp]]-rho[uind[jwp]-1])/ gamma
    ##                 end
    ##             end;
    ##             loctmp ) for wp in workers()}))
    ######################################################################

    ######################################################################
    ## OLD SERIAL
    for i=2:(N-1)
        tmp +=  dt * ( rho[i] - exp( -q[i] ) )^2     / (2.0*gamma*K) -
                dt * q[i] * ( (rho[i]-rho[i-1])/dt ) / gamma
    end
    ######################################################################

    # Potential, boundary terms:
    # ---------------------------
    
    tmp1 = ( exp(- q[N]) +  q[N]*rho[N-1] - exp(- q[1]) - q[1]*rho[1] )/gamma 
    
    out = tmp + tmp1
    
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

    ## dout_u = distribute(out_u)    

    return (V_slow_der_theta(theta),out_u)
    # return [V_slow_der_theta_2, out_u]

end

##
##
## ============================================================================================
##
##

function RESPA(theta, u, q, p, mp, dtau, nn, counter_respa, time_respa_f,time_respa_f1,time_respa_f2,time_respa_f3,time_respa_f4,der_out)               # (Tuckerman et al., JCP 97(3), 1992 )
    
    # Define small time step:
    # ---------------------------

    deltatau = dtau / nn

    # First long-range step (dtau):
    # ---------------------------

    force_theta, force = V_slow_der(theta,q)
    #### BE CAREFUL: the sign of the force is not correct!
    #### It is not possible to write -DArray.
    #### We will have to use a - sign with every occurrence of force or force_theta in the following

    ## if  (force.indexes) != (p.indexes) 
    ##     error("Error, the indexes of Force and p are not distributed correctly !")
    ## end
    for i = 1:s 
        p_theta[i] = p_theta[i] + (dtau/2) * (-force_theta[i])
    end
    for wp = workers()                          # parallel on worker processes 
        @spawnat wp ( pind = [localindexes(p)[1]];
            for jwp = 1:size(localpart(p),1) 
                localpart(p)[jwp] += (dtau/2) * (-force[pind[jwp]])
            end
        )
    end
    ##########################################################
    ## OLD SERIAL VERSION
    ## for i = 1:(N+s) 
    ##     p[i] = p[i] + (dtau/2) * force[i]
    ## end
    ##########################################################

    # Short-range steps (nn*deltatau):
    # ---------------------------
    time_fast_respa = time()
    for counter_fast_respa = 1:nn          # Verlet integrator (Tuckerman et al., JCP 97(3), 1992, eq. 2.17)
        
        time_fast_respa_0 = time()
        force_theta_old, force_old = V_fast_der(theta,u,der_out)
        #### BE CAREFUL: the sign of the force is not correct!
        #### It is not possible to write -DArray.
        #### We will have to use a - sign with every occurrence of the force in the following
        time_respa_f1[counter_respa,counter_fast_respa] = time()-time_fast_respa_0
        time_fast_respa_1 = time()
        for i=1:s
            theta[i] = theta[i] + deltatau*(p_theta[i]+(deltatau/2)*(-force_theta_old[i]))/mp_theta[i]
        end

        for wp = workers()                          # parallel on worker processes 
            @spawnat wp (
                for jwp = 1:size(localpart(u),1) 
                    localpart(u)[jwp] += deltatau*(localpart(p)[jwp] +
                        (deltatau/2)*(-localpart(force_old)[jwp]) )/localpart(mp)[jwp]
                end
            )
        end
        ##########################################################
        ## OLD SERIAL VERSION
        ## for i=1:N
        ##     u[i] = u[i]  + deltatau * ( p[i+s]  + (deltatau/2) * force_old[i+s] )/mp[i+s]
        ## end
        ##########################################################
        time_respa_f2[counter_respa,counter_fast_respa] = time()-time_fast_respa_1
        time_fast_respa_2 = time()
        force_theta_new, force_new = V_fast_der(theta,u,der_out)
        #### BE CAREFUL: the sign of the force is not correct!
        #### It is not possible to write -DArray.
        #### We will have to use a - sign with every occurrence of the force in the following
        time_respa_f3[counter_respa,counter_fast_respa] = time()-time_fast_respa_2
        time_fast_respa_3 = time()
        for i=1:s 
            p_theta[i] += (deltatau/2)*(-force_theta_old[i]-force_theta_new[i])
        end
        for wp = workers()                          # parallel on worker processes 
            @spawnat wp (
                for jwp = 1:size(localpart(p),1) 
                    localpart(p)[jwp] += (deltatau/2)*(-localpart(force_old)[jwp] -
                        localpart(force_new)[jwp])
                end
            )
        end
        ##########################################################
        ## OLD SERIAL VERSION
        ## for i=1:(N+s) 
        ##     p[i] = p[i]  + (deltatau/2)*( force_old[i] + force_new[i] )
        ## end
        ##########################################################
        time_respa_f4[counter_respa,counter_fast_respa] = time()-time_fast_respa_3
    end
    time_respa_f[counter_respa] = time()-time_fast_respa

    # Back transformations (u -> q) to update the {q} variables :
    # (Tuckerman et al., JCP 99 (1993), eq. 2.19)
    # ---------------------------
    
    #for wp = workers()                          # parallel on worker processes 
    #    @spawnat wp (uind = [localindexes(u)[1]];
    #        for jwp = 1:size(localpart(u),1) 
    #            if in(uind[jwp],ty)
    #                localpart(q)[jwp] = localpart(u)[jwp]
    #            end
    #        end
    #    )
    #end

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

    force_theta, force = V_slow_der(theta,q)
    #### BE CAREFUL: the sign of the force is not correct!
    #### It is not possible to write -DArray.
    #### We will have to use a - sign with every occurrence of the force in the following

    for i = 1:s 
        p_theta[i] = p_theta[i] + (dtau/2) * (-force_theta[i])
    end
    for wp = workers()                          # parallel on worker processes 
        @spawnat wp ( pind = [localindexes(p)[1]];
            for jwp = 1:size(localpart(p),1) 
                localpart(p)[jwp] += (dtau/2) * (-force[pind[jwp]])
            end
        )
    end

end