##########################################################################
## Potentials as functions
##########################################################################
function V_N_fun(theta,u)   
    
    # Variables:
    # ---------------------------

    bet = theta[1]
    gam = theta[2]
    
    K   = T*gam/(bet^2)

    out = 0.0
    for s = 1:n
        for k = 2:j
            out += (1/2) * (T/dt) * k/(k-1) * u[(s-1)*j+k]^2
        end
    end

    return out
end
##########################################################################
function V_n_fun(theta,u)   
    
    # Variables:
    # ---------------------------

    bet = theta[1]
    gam = theta[2]
    
    K   = T*gam/(bet^2)

    out = (bq[n+1]-bet*u[n*j+1])^2 / (2*sigma^2)
    for s = 1:n
        out += (1/2) * (T/dt) * (1/j) * ( u[(s-1)*j+1] - u[s*j+1] )^2 -
                (bq[s]-bet*u[(s-1)*j+1])^2 / (2*sigma^2)
    end

    return out
end
##########################################################################
function V_1_fun(theta,u)   
    
    # Variables:
    # ---------------------------

    bet = theta[1]
    gam = theta[2]
    
    K   = T*gam/(bet^2)

    rho = (2+gam)*bet/(2*gam)

    out = (1/gam)*exp(- bet*u[N]) + u[N]* ((T/bet) * lnr_der[N] + rho ) -
            (1/gam)*exp(- bet*u[1]) - u[1]* ((T/bet) * lnr_der[2] + rho ) #=+
            (dt/T) * (1/2) * ( (T/bet) * lnr_der[2] + rho -
            (bet/gam)*exp(-bet*u[1]) )^2 - (dt/T) * (1/2) * (bet^2/gam) * exp(-bet*u[1]) -
            (dt/T) * T * u[1] * (1/dt) * (T/bet) * (lnr_der[3] - lnr_der[2])=#

    for s = 1:(n-1)
        out +=  (dt/T) * (1/2) * ( ((T/bet) * lnr_der[s*j+1] + rho ) - (bet/gam)*exp(-bet*u[s*j+1]) )^2  -
                (dt/T) * (1/2) * (bet^2/gam) * exp(-bet*u[s*j+1]) -
                (dt/T) * T * u[s*j+1] * (1/dt) * (T/bet) * (lnr_der[s*j+1] - lnr_der[s*j]) 
    end

    for s = 1:n
        for k = 2:j
            tmp = (j-k+1)*u[(s-1)*j+1]/j
            for l = k:(j+1)
                tmp += (k-1)*u[(s-1)*j+l]/(l-1)
            end
            out += (dt/T) * (1/2) * ( ((T/bet) * lnr_der[(s-1)*j+k] + rho ) - (bet/gam)*exp(-bet*tmp) )^2  -
                    (dt/T) * (1/2) * (bet^2/gam) * exp(-bet*tmp) -
                    (dt/T) * T * tmp * (1/dt) * (T/bet) * (lnr_der[(s-1)*j+k] - lnr_der[(s-1)*j+k-1]) 
        end
    end

    return out
end
##########################################################################
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
## Potentials as expressions
##########################################################################
V_N = quote     
    bet = theta[1]
    gam = theta[2]
    K   = T*gam/(bet^2)
    out_N    = 0.0
    out_temp = 0.0
    for s = 1:n
        for k = 2:j
            out_temp += (1/2) * (T/dt) * k/(k-1) * u[(s-1)*j+k]^2
        end
    end
    out_N = out_temp
end                         
##########################################################################
V_n = quote
    bet = theta[1]
    gam = theta[2]
    K   = T*gam/(bet^2)
    out_n    = 0.0
    out_temp = (bq[n+1]-bet*u[n*j+1])^2 / (2*sigma^2)
    for s = 1:n
        out_temp += (1/2) * (T/dt) * (1/j) * ( u[(s-1)*j+1] - u[s*j+1] )^2 -
                (bq[s]-bet*u[(s-1)*j+1])^2 / (2*sigma^2)
    end
    out_n = out_temp
end
##########################################################################
V_1 = quote
    bet = theta[1]
    gam = theta[2]
    K   = T*gam/(bet^2)
    rho = (2+gam)*bet/(2*gam)
    out_1    = 0.0
    out_temp = (1/gam)*exp(- bet*u[N]) + u[N]* ((T/bet) * lnr_der[N] + rho ) -
            (1/gam)*exp(- bet*u[1]) - u[1]* ((T/bet) * lnr_der[2] + rho ) #=+
            (dt/T) * (1/2) * ( (T/bet) * lnr_der[2] + rho -
            (bet/gam)*exp(-bet*u[1]) )^2 - (dt/T) * (1/2) * (bet^2/gam) * exp(-bet*u[1]) -
            (dt/T) * T * u[1] * (1/dt) * (T/bet) * (lnr_der[3] - lnr_der[2])=#
    for s = 1:(n-1)
        out_temp += (dt/T) * (1/2) * ( ((T/bet) * lnr_der[s*j+1] + rho ) - (bet/gam)*exp(-bet*u[s*j+1]) )^2  -
                    (dt/T) * (1/2) * (bet^2/gam) * exp(-bet*u[s*j+1]) -
                    (dt/T) * T * u[s*j+1] * (1/dt) * (T/bet) * (lnr_der[s*j+1] - lnr_der[s*j]) 
    end
    for s = 1:n
        for k = 2:j
            tmp = (j-k+1)*u[(s-1)*j+1]/j
            for l = k:(j+1)
                tmp += (k-1)*u[(s-1)*j+l]/(l-1)
            end
            out_temp += (dt/T) * (1/2) * ( ((T/bet) * lnr_der[(s-1)*j+k] + rho ) - (bet/gam)*exp(-bet*tmp) )^2  -
                        (dt/T) * (1/2) * (bet^2/gam) * exp(-bet*tmp) -
                        (dt/T) * T * tmp * (1/dt) * (T/bet) * (lnr_der[(s-1)*j+k] - lnr_der[(s-1)*j+k-1]) 
        end
    end
    out_1 = out_temp
end
##
##
##########################################################################
## Derivatives
##########################################################################
# dV_N_theta = rdiff(V_N, outsym=:out_N, theta = ones(Float64, params)) 
# dV_N_u     = rdiff(V_N, outsym=:out_N, u     = ones(Float64, N)) 
# @eval dV_N(theta,u) = [($dV_N_theta)[2][1],($dV_N_u)[2][1]]

dV_n_theta = rdiff(V_n, outsym=:out_n, theta = ones(Float64, params)) 
dV_n_u     = rdiff(V_n, outsym=:out_n, u     = ones(Float64, N)) 

dV_1_theta = rdiff(V_1, outsym=:out_1, theta = ones(Float64, params)) 
dV_1_u     = rdiff(V_1, outsym=:out_1, u     = ones(Float64, N)) 

@eval dV_n(theta,u) = [($dV_n_theta)[2][1],($dV_n_u)[2][1]]
@eval dV_1(theta,u) = [($dV_1_theta)[2][1],($dV_1_u)[2][1]]
# or
@eval dV(theta,u) = [($dV_n_theta)[2][1]+($dV_1_theta)[2][1],($dV_n_u)[2][1]+($dV_1_u)[2][1]]

##
##
##########################################################################
## Napa
##########################################################################

function napa(theta, u, counter)     # (Tuckerman et al., JCP 97(3), 1992 )
    
    # Fast outer propagator (V_N), step dtau/2 
    ##########################################
    timef = time()
    for s = 1:n 
        for k = 2:j 
            u_old = u[(s-1)*j+k]
            (u[(s-1)*j+k] *= cos(w_stg(k)*dtau/2.0)) + p[params+(s-1)*j+k]*sin(w_stg(k)*dtau/2.0) / (m_stg*w_stg(k))
            (p[params+(s-1)*j+k] *= cos(w_stg(k)*dtau/2.0)) - m_stg*w_stg(k)*u_old*sin(w_stg(k)*dtau/2.0)
        end
    end
    time_respa_f[counter-nsample_burnin] += (time()-timef)

    # Slow inner propagator (V_n, V_1), step dtau
    #############################################
    time1 = time()
    force_old = -dV(theta,u)
    for s = 1:(n+1)
        u[(s-1)*j+1] += dtau * ( p[params+(s-1)*j+1]  + (dtau/2.0) * force_old[params+(s-1)*j+1] ) / m_bdy
    end 
    for i = 1:params
        theta[i] += dtau * ( p[i]  + (dtau/2.0) * force_old[i] ) / mp[i]
    end
    force_new = -dV(theta,u)
    for i = 1:(params+N)
        p[i] += (dtau/2)*( force_old[i] + force_new[i] )
    end  

    time_respa_s[counter-nsample_burnin] += (time()-time1)

    # Again, fast outer propagator (V_N), step dtau/2 
    #################################################
    timef = time()
    for s = 1:n 
        for k = 2:j 
            u_old = u[(s-1)*j+k]
            (u[(s-1)*j+k] *= cos(w_stg(k)*dtau/2.0)) + p[params+(s-1)*j+k]*sin(w_stg(k)*dtau/2.0) / (m_stg*w_stg(k))
            (p[params+(s-1)*j+k] *= cos(w_stg(k)*dtau/2.0)) - m_stg*w_stg(k)*u_old*sin(w_stg(k)*dtau/2.0)
        end
    end
    time_respa_f[counter-nsample_burnin] += (time()-timef)
    
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

end