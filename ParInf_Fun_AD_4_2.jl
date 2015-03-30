##########################################################################
## Conversion DArray -> local Array 
##########################################################################
function convert_du(du::DArray)
    return float64(convert(Array,du))
end

##########################################################################
## Conversion DArray -> local Array 
##########################################################################
function reduce_array(du_array::Array)
    a::Array{Float64,1} = du_array[chk_indexes[1][1]]
    for p = 2:num_workers
        a = append!(a,du_array[(chk_indexes[p][1][2]):(chk_indexes[p][1][end])])
    end
    return a
end

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

function V_fast_fun_proc1(theta::Array,du_array::Array)
    K = T*exp(theta[1])
    g =   exp(theta[2])
    out_fast = ( y[1] - r[1]*exp(du_array[1]) )^2/(2*sigma^2)
    out_temp = 0.0
    for p = 1:num_workers
        for i = bdy_indexes[p] 
            out_temp += K * ( du_array[(i-1)*j+p] - du_array[i*j+p] )^2/( j*dt ) / ( 2*g ) +
                            ( y[i+1] - r[i*j+1]*exp(du_array[i*j+p]) )^2/(2*sigma^2)
        end
        for i = bdy_indexes[p] 
            for k = 2:j
                out_temp += k*du_array[(i-1)*j+k+p-1]^2/( (k-1)*dt )*K/(2*g)
            end
        end
    end
    out_fast += out_temp
    return out_fast
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
## Fast potential as an expression (-> array of expressions)
##########################################################################

V_fast            = DArray(I->[:($i) for i in I[1]],(nworkers(),),workers(),[nworkers()])
V_fast_temp       = DArray(I->[:($i) for i in I[1]],(nworkers(),),workers(),[nworkers()])
V_fast_der_u_expr = DArray(I->[:($i) for i in I[1]],(nworkers(),),workers(),[nworkers()])

V_fast_proc1 = quote
    K = T*exp(theta[1])
    g =   exp(theta[2])
    out_fast = ( y[1] - r[1]*exp(du_array[1]) )^2/(2*sigma^2)
    out_temp = 0.0
    for p = 1:num_workers
        for i = bdy_indexes[p] 
            out_temp += K * ( du_array[(i-1)*j+p] - du_array[i*j+p] )^2/( j*dt ) / ( 2*g ) +
                            ( y[i+1] - r[i*j+1]*exp(du_array[i*j+p]) )^2/(2*sigma^2)
        end
        for i = bdy_indexes[p] 
            for k = 2:j
                out_temp += k*du_array[(i-1)*j+k+p-1]^2/( (k-1)*dt )*K/(2*g)
            end
        end
    end
    out_fast += out_temp
end

@spawnat first_worker ( 

    lastind  = length(bdy_indexes[1]);

    localpart(V_fast)[1] = quote
        K = T*exp(theta[1])
        g =   exp(theta[2])

        out_fast = ( localpart(dy)[1] - localpart(dr)[1]*exp(localpart(du)[1]) )^2/(2*sigma^2)
        out_temp = 0.0
        for i = 1:($lastind)
            out_temp += K * ( localpart(du)[(i-1)*j+1] - localpart(du)[i*j+1] )^2/( j*dt ) / ( 2*g ) +
                                ( localpart(dy)[i+1] - localpart(dr)[i*j+1]*exp(localpart(du)[i*j+1]) )^2/(2*sigma^2)
        end
        for i = 1:($lastind)
            for k = 2:j
                out_temp += k*localpart(du)[(i-1)*j+k]^2/( (k-1)*dt )*K/(2*g)
            end
        end
        out_fast += out_temp
    end
)

for p = 2:num_workers
    @spawnat workers()[p] (
        
        lastind  = length(bdy_indexes[p]);

        localpart(V_fast)[1] = quote
            K = T*exp(theta[1])
            g =   exp(theta[2])
    
            out_fast = 0.0
            out_temp = 0.0
            for i = 1:($lastind)
                out_temp += K * ( localpart(du)[(i-1)*j+1] - localpart(du)[i*j+1] )^2/( j*dt ) / ( 2*g ) +
                                    ( localpart(dy)[i+1] - localpart(dr)[i*j+1]*exp(localpart(du)[i*j+1]) )^2/(2*sigma^2) 
            end
            for i = 1:($lastind)
                for k = 2:j
                    out_temp += k*localpart(du)[(i-1)*j+k]^2/( (k-1)*dt )*K/(2*g)
                end
            end
            out_fast += out_temp 
        end
    )
end

function V_fast_eval()
    return sum([@fetchfrom workers()[p] (eval(localpart(V_fast)[1])) for p = 1:num_workers])
end

# @time(begin  
#     for jj = 1:100
#         eval(V_fast_proc1)
#     end
# end)
# @time(begin
#     for jj = 1:100
#         V_fast_eval()
#     end
# end)
# @time(begin  # This is much faster!
#     for jj = 1:100
#         V_fast_fun_proc1(theta,du_array)
#     end
# end)


##########################################################################
## Derivatives of V_fast 
##########################################################################

## The potential is re-written in the variables u, y and r instead of the corresponding local parts du, dy and dr
## Local parts will be re-introduced after applying rdiff()
## ## Be careful! ##, the indeces are untouched and therefore already set to be used again with localpart()  

for p = 1:num_workers
    @spawnat workers()[p] (
            temp  = string(localpart(V_fast)[1]);
            temp1 = replace(temp,  "localpart(du)", "u");
            temp2 = replace(temp1, "localpart(dy)", "y");
            temp3 = replace(temp2, "localpart(dr)", "r");
            localpart(V_fast_temp)[1] = parse(temp3);
    )
end


for p = 1:num_workers
    @spawnat workers()[p] (
        begin
            localpart(V_fast_der_u_expr)[1] = rdiff(localpart(V_fast_temp)[1], outsym=:out_fast, u     = ones(Float64, N));
        end)
end

## The expressions of the derivatives are re-parsed to re-introduce localparts()
for p = 1:num_workers
    @spawnat workers()[p] (
            temp = string(localpart(V_fast_der_u_expr)[1]);
            temp1 = replace(temp, "u", "localpart(du)");
            temp2 = replace(temp1, "y", "localpart(dy)");
            temp3 = replace(temp2, "r[", "localpart(dr)[");
            localpart(V_fast_der_u_expr)[1] = parse(temp3);
    )
end
## Derivatives can be evaluated
# @fetchfrom 2 eval(localpart(V_fast_der_u_expr)[1])[2][1]
# @fetchfrom 3 eval(localpart(V_fast_der_u_expr)[1])[2][1]
# @fetchfrom 4 eval(localpart(V_fast_der_u_expr)[1])[2][1]

## Derivatives in {u}
for p = 1:num_workers
    @spawnat workers()[p] @eval V_fast_der_u() = $(localpart(V_fast_der_u_expr)[1])[2][1]
end

## Derivatives in theta
res_fast = rdiff(V_fast_proc1, outsym=:out_fast, theta = ones(Float64, s))
@eval V_fast_der_theta() = $(res_fast)[2][1] 

## Derivatives can be evaluated
## V_fast_der_theta() = sum([@fetchfrom p eval(localpart(V_fast_der_t_expr)[1])[2][1] for p in workers()])

# @time(begin
#     for jj = 1:100
#         eval(res_fast)[2][1]
#     end    
# end)
# @time(begin   # This is the fastest!
#     for jj = 1:100
#         V_fast_der_theta()
#     end    
# end)



@everywhere function fast_force_borders()
    return [(@fetchfrom workers()[1] -V_fast_der_u()[1]);
    [(@fetchfrom workers()[p] -V_fast_der_u()[end]) + (@fetchfrom workers()[p+1] -V_fast_der_u()[1]) for p = 1:(num_workers-1)];
    (@fetchfrom workers()[end] -V_fast_der_u()[end])]
end

@everywhere function fast_force_eval(darray_force::DArray)
    force_borders_temp = fast_force_borders()
    for p = 1:num_workers
        @spawnat workers()[p] (begin
            localpart(darray_force)[1]         = force_borders_temp[p]
            localpart(darray_force)[2:(end-1)] = -V_fast_der_u()[2:(end-1)]
            localpart(darray_force)[end]       = force_borders_temp[p+1]
        end)
    end
end

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

function V_slow_fun_proc1(theta::Array,du_array::Array) 
    K  = T*exp(theta[1])
    g  =   exp(theta[2])

    out_slow = ( - exp(- du_array[1]) - du_array[1]*(K * lnr_der[1]   + g + 1) )/g +
                ( exp(- du_array[du_size]) + du_array[du_size]*(K * lnr_der[N-1] + g + 1) )/g
    
    for p = 1:(num_workers-1)
        for i = bdy_indexes[p]
            out_slow +=  dt * ( (K * lnr_der[i*j+1] + g + 1) - exp( -du_array[i*j+p] ) )^2 / (2.0*g*K) -
                         dt * du_array[i*j+p] * ( K*(lnr_der[i*j+1]-lnr_der[i*j])/dt ) / g
        end
    end
    last_p = num_workers
    for i = bdy_indexes[last_p][1]:(bdy_indexes[last_p][2]-1)
        out_slow +=  dt * ( (K * lnr_der[i*j+1] + g + 1) - exp( -du_array[i*j+last_p] ) )^2 / (2.0*g*K) -
                         dt * du_array[i*j+last_p] * ( K*(lnr_der[i*j+1]-lnr_der[i*j])/dt ) / g
    end
    for p = 1:num_workers
        out_temp = 0.0
        for i = bdy_indexes[p]
            for k = 2:j
                tmp = (j-k+1)*du_array[(i-1)*j+p]/j
                for l = k : (j+1)
                    tmp += (k-1)*du_array[(i-1)*j+l+p-1]/(l-1)
                end
                out_temp += dt * ( (K * lnr_der[(i-1)*j+k] + g + 1) - exp( -tmp ) )^2 / (2.0*g*K) -
                                dt * tmp * ( K*(lnr_der[(i-1)*j+k]-lnr_der[(i-1)*j+k-1])/dt ) / g
            end
        end
        out_slow += out_temp
    end
    
    out_slow -= ( T * (2.0+g) / (4.0*K) + (N/2.0) * (theta[1]-theta[2]) -
                (0.5*theta[1]^2+theta[1]*(log(T)-5.0)) - (theta[2]^2+theta[2]) )
    return out_slow
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
V_slow_temp       = DArray(I->[:($i) for i in I[1]],(nworkers(),),workers(),[nworkers()])
V_slow_der_u_expr = DArray(I->[:($i) for i in I[1]],(nworkers(),),workers(),[nworkers()])

V_slow_proc1 = quote
    K  = T*exp(theta[1])
    g  =   exp(theta[2])

    out_slow = ( - exp(- du_array[1]) - du_array[1]*(K * lnr_der[1]   + g + 1) )/g +
                ( exp(- du_array[$du_size]) + du_array[$du_size]*(K * lnr_der[N-1] + g + 1) )/g
    
    for p = 1:(num_workers-1)
        for i = bdy_indexes[p]
            out_slow +=  dt * ( (K * lnr_der[i*j+1] + g + 1) - exp( -du_array[i*j+p] ) )^2 / (2.0*g*K) -
                         dt * du_array[i*j+p] * ( K*(lnr_der[i*j+1]-lnr_der[i*j])/dt ) / g
        end
    end
    last_p = num_workers
    for i = bdy_indexes[last_p][1]:(bdy_indexes[last_p][2]-1)
        out_slow +=  dt * ( (K * lnr_der[i*j+1] + g + 1) - exp( -du_array[i*j+last_p] ) )^2 / (2.0*g*K) -
                         dt * du_array[i*j+last_p] * ( K*(lnr_der[i*j+1]-lnr_der[i*j])/dt ) / g
    end
    for p = 1:num_workers
        out_temp = 0.0
        for i = bdy_indexes[p]
            for k = 2:j
                tmp = (j-k+1)*du_array[(i-1)*j+p]/j
                for l = k : (j+1)
                    tmp += (k-1)*du_array[(i-1)*j+l+p-1]/(l-1)
                end
                out_temp += dt * ( (K * lnr_der[(i-1)*j+k] + g + 1) - exp( -tmp ) )^2 / (2.0*g*K) -
                                dt * tmp * ( K*(lnr_der[(i-1)*j+k]-lnr_der[(i-1)*j+k-1])/dt ) / g
            end
        end
        out_slow += out_temp
    end
    
     out_slow -= ( T * (2.0+g) / (4.0*K) + (N/2.0) * (theta[1]-theta[2]) -
                (0.5*theta[1]^2+theta[1]*(log(T)-5.0)) - (theta[2]^2+theta[2]) )

end


@spawnat first_worker ( 

    lastind  = length(bdy_indexes[1]);

    localpart(V_slow)[1] = quote

        K  = T*exp(theta[1])
        g  =   exp(theta[2])

        # Potential, boundary terms:
        # ---------------------------
    
        out_slow = ( - exp(- localpart(du)[1]) - localpart(du)[1]*(K * localpart(dlnr_der)[1]   + g + 1) )/g
    
        # Potential, bulk terms:
        # ---------------------------
    
        for i = 1:($lastind)
            out_slow +=  dt * ( (K * localpart(dlnr_der)[i*j+1] + g + 1) - exp( -localpart(du)[i*j+1] ) )^2 / (2.0*g*K) -
                            dt * localpart(du)[i*j+1] * ( K*(localpart(dlnr_der)[i*j+1]-localpart(dlnr_der)[i*j])/dt ) / g
        end

        out_temp = 0.0
        for i = 1:($lastind)
            for k = 2:j
                tmp = (j-k+1)*localpart(du)[(i-1)*j+1]/j
                for l = k : (j+1)
                    tmp += (k-1)*localpart(du)[(i-1)*j+l]/(l-1)
                end
                out_temp += dt * ( (K * localpart(dlnr_der)[(i-1)*j+k] + g + 1) - exp( -tmp ) )^2 / (2.0*g*K) -
                                dt * tmp * ( K*(localpart(dlnr_der)[(i-1)*j+k]-localpart(dlnr_der)[(i-1)*j+k-1])/dt ) / g
            end
        end
        out_slow += out_temp
    
    end
)
@spawnat last_worker (

    lastind  = length(bdy_indexes[end]);
    lastel   = chk_sizes[end];

    localpart(V_slow)[1] = quote

        K  = T*exp(theta[1])
        g  = exp(theta[2])

        # Potential, boundary terms:
        # ---------------------------
        
        out_slow = ( exp(- localpart(du)[$lastel]) + localpart(du)[$lastel]*(K * localpart(dlnr_der)[$lastel-1] + g + 1) )/g
        
        # Potential, bulk terms:
        # ---------------------------
        
        for i = 1:(($lastind)-1)
            out_slow +=  dt * ( (K * localpart(dlnr_der)[i*j+1] + g + 1) - exp( -localpart(du)[i*j+1] ) )^2 / (2.0*g*K) -
                            dt * localpart(du)[i*j+1] * ( K*(localpart(dlnr_der)[i*j+1]-localpart(dlnr_der)[i*j])/dt ) / g
        end

        out_temp = 0.0
        for i = 1:($lastind)
            for k = 2:j
                tmp = (j-k+1)*localpart(du)[(i-1)*j+1]/j
                for l = k : (j+1)
                    tmp += (k-1)*localpart(du)[(i-1)*j+l]/(l-1)
                end
                out_temp += dt * ( (K * localpart(dlnr_der)[(i-1)*j+k] + g + 1) - exp( -tmp ) )^2 / (2.0*g*K) -
                                dt * tmp * ( K*(localpart(dlnr_der)[(i-1)*j+k]-localpart(dlnr_der)[(i-1)*j+k-1])/dt ) / g
            end
        end
        out_slow += out_temp

        out_slow -= ( T * (2.0+g) / (4.0*K) + (N/2.0) * (theta[1]-theta[2]) -
                (0.5*theta[1]^2+theta[1]*(log(T)-5.0)) - (theta[2]^2+theta[2]) )

    end
)

for p = 2:(num_workers-1)

    @spawnat workers()[p] (

        lastind  = length(bdy_indexes[p]);
        
        localpart(V_slow)[1] = quote
            
            K  = T*exp(theta[1])
            g  = exp(theta[2])

            
            # Potential, bulk terms:
            # ---------------------------

            out_slow = 0.0
            
            for i = 1:($lastind)
                out_slow +=  dt * ( (K * localpart(dlnr_der)[i*j+1] + g + 1) - exp( -localpart(du)[i*j+1] ) )^2 / (2.0*g*K) -
                                dt * localpart(du)[i*j+1] * ( K*(localpart(dlnr_der)[i*j+1]-localpart(dlnr_der)[i*j])/dt ) / g
            end

            out_temp = 0.0
            for i = 1:($lastind)
                for k = 2:j
                    tmp = (j-k+1)*localpart(du)[(i-1)*j+1]/j
                    for l = k : (j+1)
                        tmp += (k-1)*localpart(du)[(i-1)*j+l]/(l-1)
                    end
                    out_temp += dt * ( (K * localpart(dlnr_der)[(i-1)*j+k] + g + 1) - exp( -tmp ) )^2 / (2.0*g*K) -
                                    dt * tmp * ( K*(localpart(dlnr_der)[(i-1)*j+k]-localpart(dlnr_der)[(i-1)*j+k-1])/dt ) / g
                end
            end
            out_slow += out_temp
    
        end 
    )
end

function V_slow_eval()
    return sum([@fetchfrom workers()[p] (eval(localpart(V_slow)[1])) for p = 1:num_workers])
end

# @time(begin  # This is faster
#     for jj = 1:100
#         eval(V_slow_proc1)
#     end
# end)
# @time(begin
#     for jj = 1:100
#         V_slow_eval()
#     end
# end)
# @time(begin  # This is much faster!
#     for jj = 1:100
#         V_slow_fun_proc1(theta,du_array)
#     end
# end)

##########################################################################
## Derivatives of V_slow
##########################################################################

## The potential is re-written in the variables u and lnr_der instead of the corresponding local parts du, dlnr_der
## Local parts will be re-introduced after applying rdiff()
## ## Be careful! ##, the indeces are untouched and therefore already set to be used again with localpart()  

for p = 1:num_workers
    @spawnat workers()[p] (
            temp  = string(localpart(V_slow)[1]);
            temp1 = replace(temp,  "localpart(du)", "u");
            temp2 = replace(temp1, "localpart(dlnr_der)", "lnr_der");
            localpart(V_slow_temp)[1] = parse(temp2);
    )
end

for p = 1:num_workers
    @spawnat workers()[p] (
        begin
            localpart(V_slow_der_u_expr)[1] = rdiff(localpart(V_slow_temp)[1], outsym=:out_slow, u     = ones(Float64, N));
        end)
end

## The expressions of the derivatives are re-parsed to re-introduce localparts()
for p = 1:num_workers
    @spawnat workers()[p] (
            temp = string(localpart(V_slow_der_u_expr)[1]);
            temp1 = replace(temp,  "u", "localpart(du)");
            temp2 = replace(temp1, "lnr_der", "localpart(dlnr_der)");
            localpart(V_slow_der_u_expr)[1] = parse(temp2);
    )
end

## Derivatives in {u}
for p = 1:num_workers
    @spawnat workers()[p] @eval V_slow_der_u() = $(localpart(V_slow_der_u_expr)[1])[2][1]
end

## Derivatives in theta
res_slow = rdiff(V_slow_proc1, outsym=:out_slow, theta = ones(Float64, s))
@eval V_slow_der_theta() = $(res_slow)[2][1] 


@everywhere function slow_force_borders()
    return [(@fetchfrom workers()[1] -V_slow_der_u()[1]);
    [(@fetchfrom workers()[p] -V_slow_der_u()[end]) + (@fetchfrom workers()[p+1] -V_slow_der_u()[1]) for p = 1:(num_workers-1)];
    (@fetchfrom workers()[end] -V_slow_der_u()[end])]
end

@everywhere function slow_force_eval(darray_force::DArray)
    force_borders_temp = slow_force_borders()
    for p = 1:num_workers
        @spawnat workers()[p] (begin
            localpart(darray_force)[1]         = force_borders_temp[p]
            localpart(darray_force)[2:(end-1)] = -V_slow_der_u()[2:(end-1)]
            localpart(darray_force)[end]       = force_borders_temp[p+1]
        end)
    end
end

##
##
##########################################################################
## RESPA
##########################################################################

function RESPA(theta, ptheta, du, dpu, force, force_old, force_new, dmu, dtau, nn, deltatau, counter)     # (Tuckerman et al., JCP 97(3), 1992 )
    
    # First long-range step (dtau):
    # ---------------------------
    time1 = time()
   
    ptheta += (dtau/2)*(-V_slow_der_theta());
    slow_force_eval(force)
    @time(begin
        for pid in workers()
            @spawnat pid (localpart(dpu)[:] += (dtau / 2.0) * localpart(force)[:])
        end
    end)
    apu = convert(Array,dpu)
    aforce = convert(Array, force)
    @time(begin
        for i=1:du_size 
            apu[i] += (dtau/2.0)* aforce[i]
        end 
    end)

    time_respa_s[counter-nsample_burnin] += (time()-time1) 
    # Short-range steps (nn*deltatau):
    # ---------------------------

    for counter_fast_respa = 1:nn                       # Verlet integrator 
        timefd = time()
        
        force_old_theta = -V_fast_der_theta() 
        fast_force_eval(force_old)

        # force_old = -V_fast_der(theta,u)               # (Tuckerman et al., JCP 97(3), 1992, eq. 2.17)
        
        time_respa_f_d[counter-nsample_burnin] += (time()-timefd) 
        timefu = time()

        theta += deltatau * ( ptheta  + (deltatau/2) * force_old_theta ) ./ mtheta
        time0 = time()
        for p = 1:num_workers
            @spawnat workers()[p] (begin
                localpart(du)[:] += deltatau * ( localpart(dpu)[:] + (deltatau/2) * localpart(force_old)[:] ) ./ localpart(dmu)[:]
            end)
        end 

      
        du_array = float64(convert(Array,du))

        time_respa_f_u[counter-nsample_burnin] += (time()-timefu)
        timeglobu = time()
        # for pid in workers()
        #     remotecall(pid, x->(global u; u = x; nothing), u);
        #     remotecall(pid, x->(global theta; theta = x; nothing), theta)
        # end
        time_respa_globu[counter-nsample_burnin] += (time()-timeglobu)
        
        timefd = time()

        force_new_theta = -V_fast_der_theta()
        fast_force_eval(force_new)
        # force_new = -V_fast_der(theta,u)

        time_respa_f_d[counter-nsample_burnin] += (time()-timefd)
        timefu = time()
        ptheta += (deltatau/2)*( force_old_theta + force_new_theta )
        for p = 1:num_workers
            @spawnat workers()[p] (begin
                localpart(dpu)[:] += (deltatau/2)*(localpart(force_old)[:]+localpart(force_new)[:])
            end)
        end
        time_respa_f_u[counter-nsample_burnin] += (time()-timefu)
    
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
    ptheta += (dtau/2)*(-V_slow_der_theta())
    slow_force_eval(force)
    for p = 1:num_workers
        @spawnat workers()[p] (begin
            localpart(dpu)[:] += (dtau/2)*localpart(force)[:]
        end)
    end
    time_respa_s[counter-nsample_burnin] += (time()-time2)
end

##########################################################################
##########################################################################
## UNUSED STUFF 
##########################################################################
##########################################################################
#=
##########################################################################
## Fast potential as ... many functions
##########################################################################

@everywhere function V_fast_dat(u, y, r, i::Int64)      # i = 0:n

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

@everywhere function V_fast_ser(theta, u) 
    V_fast_sum = V_fast_dat(u, 0)
    for i = 1:n
        V_fast_sum += V_fast_bdy(theta, u, i) + V_fast_dat(u, y, r, i)
        for k = 2:j 
            V_fast_sum += V_fast_stg(theta, u, i, k)
        end
    end
    return V_fast_sum
end

V_fast_out = dzeros(num_workers)

function V_fast_par()
    @sync begin
        @spawnat first_worker (
                lastind = length(bdy_indexes[1]);
                V_fast_sum = V_fast_dat(localpart(du), localpart(dy), localpart(dr), 0);
                for i = 1:lastind
                    V_fast_sum += V_fast_bdy(theta, localpart(du), i) +
                                    V_fast_dat(localpart(du), localpart(dy), localpart(dr), i)
                    for k = 2:j 
                        V_fast_sum += V_fast_stg(theta, localpart(du), i, k)
                    end
                end;
                localpart(V_fast_out)[1] = V_fast_sum
        )
        for p = 2:num_workers
            @spawnat workers()[p] (
                lastind = length(bdy_indexes[p]);
                V_fast_sum = 0.0;
                for i = 1:lastind
                    V_fast_sum += V_fast_bdy(theta, localpart(du), i) +
                                    V_fast_dat(localpart(du), localpart(dy), localpart(dr), i)
                    for k = 2:j 
                        V_fast_sum += V_fast_stg(theta, localpart(du), i, k)
                    end
                end;
                localpart(V_fast_out)[1] = V_fast_sum
            )
        end
    end
    return sum(V_fast_out)
end

V_fast_par()


@everywhere V_fast_der_dat = rdiff(V_fast_dat,(ones(N), ones(n+1), ones(N), 1))
@everywhere V_fast_der_bdy = rdiff(V_fast_bdy,(ones(2), ones(N),            1))
@everywhere V_fast_der_stg = rdiff(V_fast_stg,(ones(2), ones(N),            1, 1))

@everywhere function dV_fast_local_sum(lastind)
    out_sum = sum([V_fast_der_dat(localpart(du), localpart(dy), localpart(dr), i)[2][1] for i = 1:lastind]) +
                            sum([V_fast_der_bdy(theta, localpart(du), i)[2][2] for i = 1:lastind]) +
                                sum([V_fast_der_stg(theta, localpart(du), i, k)[2][2] for i = 1:lastind, k = 2:j])
    return out_sum
end

@everywhere function dV_fast_local_sum_term0(lastind)
    out_zero = V_fast_der_dat(localpart(du), localpart(dy), localpart(dr), 0)[2][1]
    return out_zero
end

@everywhere function dV_fast_local_sum_first(lastind)
    out_sum = dV_fast_local_sum(lastind) + dV_fast_local_sum_term0(lastind)
    return out_sum
end

dV_fast_out=Array(Any, num_workers)
for p = 1:num_workers
    dV_fast_out[p] = zeros(chk_sizes[p])
end
dV_fast_out = distribute(dV_fast_out)

function dV_fast_u()
   @sync begin
        @async begin
            @spawnat first_worker (
                lastind = length(bdy_indexes[1]);
                localpart(dV_fast_out)[1] = dV_fast_local_sum_first(lastind);
            )
        end
        @async begin
            for p = 2:num_workers
                @spawnat workers()[p] (
                    lastind = length(bdy_indexes[p]);
                    localpart(dV_fast_out)[1] = dV_fast_local_sum(lastind);   
                )
            end
        end
    end

    # @sync begin
    #     for p = 1:(num_workers-1)
    #         @spawnat workers()[p] ( localpart(dV_fast_out)[2][end] += dV_fast_out[2*(p+1)][1])
    #     end
    # end

    # @sync begin
    #     for p = 2:(num_workers)
    #         @spawnat workers()[p] ( localpart(dV_fast_out)[2][1] = dV_fast_out[2*(p-1)][end])
    #     end
    # end

    return dV_fast_out

end

@everywhere function dV_fast_local_sum_theta(lastind)
    out_sum = sum([V_fast_der_bdy(theta, localpart(du), i)[2][1] for i = 1:lastind]) +
                sum([V_fast_der_stg(theta, localpart(du), i, k)[2][1] for i = 1:lastind, k = 2:j])
    return out_sum
end

dV_fast_out_theta=Array(Any, num_workers)
for p = 1:num_workers
    dV_fast_out_theta[p] = zeros(s)
end
dV_fast_out_theta = distribute(dV_fast_out_theta)

function dV_fast_t()
   @sync begin
        @async begin
            for p = 1:num_workers
                @spawnat workers()[p] (
                    lastind = length(bdy_indexes[p]);
                    localpart(dV_fast_out_theta)[1] = dV_fast_local_sum_theta(lastind)       
                )
            end
        end
    end

    return sum(dV_fast_out_theta)
end
=#