##########################################################################
## Fast potential as a function
##########################################################################

function V_fast_fun(theta,u)    # Fast dynamics (harmonic and data terms)
    
    # Variables:
    # ---------------------------
    
    out = 0.0

    # beta    = theta[1]
    # tau     = theta[2]
    
    K  = T*exp(theta[1])
    g  = exp(theta[2])

    out = ( y[1] - r[1]*exp(u[1]) )^2/(2*sigma^2)
    for i=1:n
        out += K * ( u[(i-1)*j+1] - u[i*j+1] )^2/( j*dt )   / ( 2*g ) +
                    ( y[i+1] - r[i*j+1]*exp(u[i*j+1]) )^2/(2*sigma^2)
        for k=2:j
            out += k*u[(i-1)*j+k]^2/( (k-1)*dt ) * K/(2*g)
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
## Fast potential as ... many functions
##########################################################################

@everywhere function V_fast_dat(u, i::Int64)      # i = 0:n

    out = ( y[i+1] - r[i*j+1]*exp(u[i*j+1]) )^2/(2*sigma^2)
    return out

end


@everywhere function V_fast_bdy(theta, u, i::Int64)    # i = 1:n

    K = T*exp(theta[1])
    g = exp(theta[2])
    
    out = K * ( u[(i-1)*j+1] - u[i*j+1] )^2/( j*dt ) / ( 2*g )
    
    return out
end

@everywhere function V_fast_stg(theta, u, i::Int64, k::Int64)  # i = 1:n, k = 2:j
    
    K = T*exp(theta[1])
    g = exp(theta[2])

    out = k*u[(i-1)*j+k]^2/( (k-1)*dt ) * K/(2*g)
    
    return out
end

@everywhere function V_fast(theta, u) 
    V_fast_sum = V_fast_dat(u, y, r, 0)
    for i = 1:n
        V_fast_sum += V_fast_bdy(theta, u, i) + V_fast_dat(u, y, r, i)
        for k = 2:j 
            V_fast_sum += V_fast_stg(theta, u, i, k)
        end
    end
    return V_fast_sum
end

out = dzeros(2)
for p in workers()
    @spawnat p localpart(out)[1] = sum([fix(localpart(d),i) for i = 1:3])
end
out

@everywhere V_fast_der_dat = rdiff(V_fast_dat,(ones(N), 1            ))
@everywhere V_fast_der_bdy = rdiff(V_fast_bdy,(ones(2), ones(N), 1   ))
@everywhere V_fast_der_stg = rdiff(V_fast_stg,(ones(2), ones(N), 1, 1))


##########################################################################
## Fast potential as an expression (-> array of expressions)
##########################################################################

V_fast            = DArray(I->[:($i) for i in I[1]],(nworkers(),),workers(),[nworkers()])
V_fast_der_u_expr = DArray(I->[:($i) for i in I[1]],(nworkers(),),workers(),[nworkers()])
V_fast_der_t_expr = DArray(I->[:($i) for i in I[1]],(nworkers(),),workers(),[nworkers()])

@spawnat first_worker ( 

    # firstind = bdy_indexes[1][1];
    # lastind  = bdy_indexes[1][end];
    lastind  = length(bdy_indexes[1]);

    localpart(V_fast)[1] = quote
        K = T*exp(theta[1])
        g =   exp(theta[2])

        out_fast = ( localpart(dy)[1] - localpart(dr)[1]*exp(localpart(du)[1]) )^2/(2*sigma^2)
        out_temp = 0.0
        for ii = 1:($lastind)
            out_temp += K * ( localpart(du)[(ii-1)*j+1] - localpart(du)[ii*j+1] )^2/( j*dt ) / ( 2*g ) +
                                ( localpart(dy)[ii+1] - localpart(dr)[ii*j+1]*exp(localpart(du)[ii*j+1]) )^2/(2*sigma^2)
        end
        for ii = 1:($lastind)
            for k = 2:j
                out_temp += k*localpart(du)[(ii-1)*j+k]^2/( (k-1)*dt )*K/(2*g)
            end
        end
        out_fast += out_temp
    end
)

for i = 2:num_workers
    @spawnat workers()[i] (
        # firstind = bdy_indexes[i][1];
        # lastind  = bdy_indexes[i][end];
        lastind  = length(bdy_indexes[i]);
        # du_local = dulocal;
        # dr_local = drlocal;
        # dy_local = dylocal;

        localpart(V_fast)[1] = quote
            K = T*exp(theta[1])
            g =   exp(theta[2])
    
            out_fast = 0.0
            out_temp = 0.0
            for ii = 1:($lastind)
                out_temp += K * ( dulocal[(ii-1)*j+1] - dulocal[ii*j+1] )^2/( j*dt ) / ( 2*g ) +
                                    ( dylocal[ii+1] - drlocal[ii*j+1]*exp(dulocal[ii*j+1]) )^2/(2*sigma^2) 
            end
            for ii = 1:($lastind)
                for k = 2:j
                    out_temp += k*dulocal[(ii-1)*j+k]^2/( (k-1)*dt )*K/(2*g)
                end
            end
            out_fast += out_temp 
        end
    )
end
V_fast
##########################################################################
## Derivatives of V_fast 
##########################################################################
for i = 1:num_workers
    @spawnat workers()[i] (
        localpart(V_fast_der_t_expr)[1] = rdiff(localpart(V_fast)[1], outsym=:out_fast, theta = ones(Float64, s));
        localpart(V_fast_der_u_expr)[1] = rdiff(localpart(V_fast)[1], outsym=:out_fast, du = ones(Float64, du_indexes[i][end]));
    )
end

V_fast_der_theta(theta) = begin
sum(
    map(
        fetch, {@spawnat p (
            @eval ($(localpart(V_fast_der_t_expr)[1]))[2][1]) for p in workers()}
        )
    )
end

V_fast_der_u(u) = begin
sum(
    map(
        fetch, {@spawnat p (
            @eval ($(localpart(V_fast_der_u_expr)[1]))[2][1]) for p in workers()}
        )
    )
end

V_fast_der(theta,u) = [V_fast_der_theta(theta), V_fast_der_u(u)]

##
##
##########################################################################
## Slow potential as a function
##########################################################################

function V_slow_fun(theta,u) 

    # Variables:
    # ---------------------------

    # beta  = theta[1]
    # tau   = theta[2]
    
    K = T*exp(theta[1])
    g = exp(theta[2])
    

    # Definition of rho:
    # ---------------------------
    
    # rho = Array(Float64,N)
    # rho = K * lnr_der + g + 1

    # Potential, boundary terms:
    # ---------------------------
    
    out = ( exp(- u[N]) + u[N]*(K * lnr_der[N-1] + g + 1) -
            exp(- u[1]) - u[1]*(K * lnr_der[1]   + g + 1) ) / g

    
    # Potential, bulk terms:
    # ---------------------------

    for i = 1:(n-1)
        out +=  dt * ( (K * lnr_der[i*j+1] + g + 1) - exp( -u[i*j+1] ) )^2 / (2.0*g*K) -
                        dt * u[i*j+1] * ( K*(lnr_der[i*j+1]-lnr_der[i*j])/dt ) / g
    end

    for i = 1:n
        for k = 2:j
            tmp = (j-k+1)*u[(i-1)*j+1]/j
            for l = k : (j+1)
                tmp += (k-1)*u[(i-1)*j+l]/(l-1)
            end
            out += dt * ( (K * lnr_der[(i-1)*j+k] + g + 1) - exp( -tmp ) )^2 / (2.0*g*K) -
                            dt * tmp * ( K*(lnr_der[(i-1)*j+k]-lnr_der[(i-1)*j+k-1])/dt ) / g
        end
    end

    # for i=2:(N-1)
    #     out +=  dt * ( rho[i] - exp( -q[i] ) )^2     / (2.0*g*K) -
    #             dt * q[i] * ( (rho[i]-rho[i-1])/dt ) / g
    # end
    
    out -= T * (2.0+g) / (4.0*K) + (N/2.0) * (theta[1]-theta[2]) -
            (0.5*theta[1]^2+theta[1]*(log(T)-5.0)) - (theta[2]^2+theta[2])
    
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
## Slow potential as an expression (-> array of expressions)
##########################################################################
V_slow            = DArray(I->[:($i) for i in I[1]],(nworkers(),),workers(),[nworkers()])
V_slow_der_u_expr = DArray(I->[:($i) for i in I[1]],(nworkers(),),workers(),[nworkers()])
V_slow_der_t_expr = DArray(I->[:($i) for i in I[1]],(nworkers(),),workers(),[nworkers()])

@spawnat first_worker ( 

    firstind = bdy_indexes[1][1];
    lastind  = bdy_indexes[1][end];

    localpart(V_slow)[1] = quote
        
        # Variables:
    # ---------------------------

    # bet   = theta[1]
    # tau   = theta[2]
    
    K  = T*exp(theta[1])
    g  = exp(theta[2])

    # Potential, boundary terms:
    # ---------------------------
    
    out_slow = ( - exp(- u[1]) - u[1]*(K * lnr_der[1]   + g + 1) )/g
    
    # Potential, bulk terms:
    # ---------------------------
    
    for i = ($firstind):($lastind)
        out_slow +=  dt * ( (K * lnr_der[i*j+1] + g + 1) - exp( -u[i*j+1] ) )^2 / (2.0*g*K) -
                        dt * u[i*j+1] * ( K*(lnr_der[i*j+1]-lnr_der[i*j])/dt ) / g
    end

    out_temp = 0.0
    for i = ($firstind):($lastind)
        for k = 2:j
            tmp = (j-k+1)*u[(i-1)*j+1]/j
            for l = k : (j+1)
                tmp += (k-1)*u[(i-1)*j+l]/(l-1)
            end
            out_temp += dt * ( (K * lnr_der[(i-1)*j+k] + g + 1) - exp( -tmp ) )^2 / (2.0*g*K) -
                            dt * tmp * ( K*(lnr_der[(i-1)*j+k]-lnr_der[(i-1)*j+k-1])/dt ) / g
        end
    end
    out_slow += out_temp

    end
)

@spawnat last_worker (

    firstind = bdy_indexes[end][1];
    lastind  = bdy_indexes[end][end];

    localpart(V_slow)[1] = quote

        # Variables:
        # ---------------------------

        # bet   = theta[1]
        # tau   = theta[2]
        
        K  = T*exp(theta[1])
        g  = exp(theta[2])

        # Potential, boundary terms:
        # ---------------------------
        
        out_slow = ( exp(- u[N]) + u[N]*(K * lnr_der[N-1] + g + 1) )/g
        
        # Potential, bulk terms:
        # ---------------------------
        
        for i = ($firstind):($lastind-1)
            out_slow +=  dt * ( (K * lnr_der[i*j+1] + g + 1) - exp( -u[i*j+1] ) )^2 / (2.0*g*K) -
                            dt * u[i*j+1] * ( K*(lnr_der[i*j+1]-lnr_der[i*j])/dt ) / g
        end

        out_temp = 0.0
        for i = ($firstind):($lastind)
            for k = 2:j
                tmp = (j-k+1)*u[(i-1)*j+1]/j
                for l = k : (j+1)
                    tmp += (k-1)*u[(i-1)*j+l]/(l-1)
                end
                out_temp += dt * ( (K * lnr_der[(i-1)*j+k] + g + 1) - exp( -tmp ) )^2 / (2.0*g*K) -
                                dt * tmp * ( K*(lnr_der[(i-1)*j+k]-lnr_der[(i-1)*j+k-1])/dt ) / g
            end
        end
        out_slow += out_temp

        out_slow -= ( T * (2.0+g) / (4.0*K) + (N/2.0) * (theta[1]-theta[2]) -
                (0.5*theta[1]^2+theta[1]*(log(T)-5.0)) - (theta[2]^2+theta[2]) )

    end
)

for i = 2:(num_workers-1)
    @spawnat workers()[i] (
        firstind = bdy_indexes[i][1];
        lastind  = bdy_indexes[i][end];
        localpart(V_slow)[1] = quote
            
            # Variables:
            # ---------------------------

            # bet   = theta[1]
            # tau   = theta[2]
            
            K  = T*exp(theta[1])
            g  = exp(theta[2])

            
            # Potential, bulk terms:
            # ---------------------------

            out_slow = 0.0
            
            for i = ($firstind):($lastind)
                out_slow +=  dt * ( (K * lnr_der[i*j+1] + g + 1) - exp( -u[i*j+1] ) )^2 / (2.0*g*K) -
                                dt * u[i*j+1] * ( K*(lnr_der[i*j+1]-lnr_der[i*j])/dt ) / g
            end

            out_temp = 0.0
            for i = ($firstind):($lastind)
                for k = 2:j
                    tmp = (j-k+1)*u[(i-1)*j+1]/j
                    for l = k : (j+1)
                        tmp += (k-1)*u[(i-1)*j+l]/(l-1)
                    end
                    out_temp += dt * ( (K * lnr_der[(i-1)*j+k] + g + 1) - exp( -tmp ) )^2 / (2.0*g*K) -
                                    dt * tmp * ( K*(lnr_der[(i-1)*j+k]-lnr_der[(i-1)*j+k-1])/dt ) / g
                end
            end
            out_slow += out_temp
    
        end 
    )
end

##########################################################################
## Derivatives of V_slow
##########################################################################
for i = 1:num_workers
    @spawnat workers()[i] (
        localpart(V_slow_der_t_expr)[1] = rdiff(localpart(V_slow)[1], outsym=:out_slow, theta = ones(Float64, s) );
        localpart(V_slow_der_u_expr)[1] = rdiff(localpart(V_slow)[1], outsym=:out_slow, u     = ones(Float64, N) );
    )
end

V_slow_der_theta(theta) = begin
sum(
    map(
        fetch, {@spawnat p (
            @eval ($(localpart(V_slow_der_t_expr)[1]))[2][1]) for p in workers()}
        )
    )
end

V_slow_der_u(u) = begin
sum(
    map(
        fetch, {@spawnat p (
            @eval ($(localpart(V_slow_der_u_expr)[1]))[2][1]) for p in workers()}
        )
    )
end

V_slow_der(theta,u) = [V_slow_der_theta(theta), V_slow_der_u(u)]

##
##
##########################################################################
## RESPA
##########################################################################

function RESPA(theta, u, p, mp, dtau, nn, deltatau,counter)     # (Tuckerman et al., JCP 97(3), 1992 )
    

    # First long-range step (dtau):
    # ---------------------------
    time1 = time()
    force = -V_slow_der(theta,u)
    for i=1:(N+s) 
        p[i] = p[i]  + (dtau/2)* force[i]
    end
    time_respa_s[counter-nsample_burnin] += (time()-time1) 
    # Short-range steps (nn*deltatau):
    # ---------------------------

    for counter_fast_respa = 1:nn                       # Verlet integrator 
        timefd = time()
        force_old = -V_fast_der(theta,u)               # (Tuckerman et al., JCP 97(3), 1992, eq. 2.17)
        time_respa_f_d[counter-nsample_burnin] += (time()-timefd) 
        timefu = time()
        for i=1:s
            theta[i] = theta[i] + deltatau * ( p[i]  + (deltatau/2) * force_old[i] )/mp[i]
        end
        time0 = time()
        for i=1:N
            u[i] = u[i]  + deltatau * ( p[i+s]  + (deltatau/2) * force_old[i+s] )/mp[i+s]
        end
        time_respa_f_u[counter-nsample_burnin] += (time()-timefu)
        timeglobu = time()
        for pid in workers()
            remotecall(pid, x->(global u; u = x; nothing), u);
            remotecall(pid, x->(global theta; theta = x; nothing), theta)
        end
        time_respa_globu[counter-nsample_burnin] += (time()-timeglobu)
        timefd = time()
        force_new = -V_fast_der(theta,u)
        time_respa_f_d[counter-nsample_burnin] += (time()-timefd)
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
    time2 = time()
    force = -V_slow_der(theta,u)
    for i=1:(N+s) 
        p[i] = p[i]  + (dtau/2)* force[i]
    end
    time_respa_s[counter-nsample_burnin] += (time()-time2)
end