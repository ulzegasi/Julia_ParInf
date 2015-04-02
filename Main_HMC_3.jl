## ============================================================================================
## Project: Parameter Inference for a Simple Stochastic Hydrological Model
## File: Main_HMC_1.jl
##
## Description: Single bucket, scale invariant noise, no evaporation yet
##              dS(t)/dt = r(t) - S(t)/k + sqrt(gamma/k) S(t) eta(t)
## The {q} variable set is used for slow (long-range) components. 
## The {u} variable set is used for fast (short-range) components.
## 
## The derivative of the potential V(theta,u/q) w.r.t. theta can be calculated explicitely
## or by automated differentiation (ForwardDiff package). 
## 
## Different masses are used for the burn-in phase.
##
## NEW IN VERSION 2. THE PRIOR FOR PARAMETER K IS A LOG-NORMAL DISTRIBUTION.
## We use the tentative function f(K) = (3 Exp(-(Log(K)-4)^2/2))/K 
##
## NEW IN VERSION 3. PARALLELIZED.
##             
## GitHub repository: https://github.com/ulzegasi/Julia_ParInf.git
##
## Authors: Simone Ulzega, Carlo Albert
##
## simone.ulzega@eawag.ch
## carlo.albert@eawag.ch
## ============================================================================================
##
##
## ============================================================================================
## Preliminaries:
## ============================================================================================
##
##

tstart = time()

## Add worker processes
nw = 2
addprocs(nw)

dir  = "C:/Users/ulzegasi/Julia_files/ParInf_HMC";       # Main project directory   
dir2 = "C:/Users/ulzegasi/SWITCHdrive/JuliaTemp/Data";   # Secondary directory


fname= string("_testing")
## Output files name
## It can be an empty string          

##
##
## ============================================================================================
## System parameters:
## ============================================================================================
##
##
range = 2:5002;

# Time points t.dat
t  = float64(readdlm("$dir/t.dat")[range,2]);

const N = 301;                      # Total number of discrete time points
const j = 30;                       # j-1 = number of staging beads per segment
const n = int64((N-1)/j);           # n + 1 = number of "end point" beads (see Tuckerman '93)
                                    # n = number of segments
                                    # n (j-1) = total number of staging beads
if  ((N-1)%j) != 0 
    error("Be careful, the number of staging points j must fulfill (N-1)/j = integer !")
end
    
t = linspace(t[1], t[end], N);      # Time points
T  = t[end]-t[1];                   # Time interval
dt = T/(N-1);                       # Time step
ty  = iround(linspace(1, N, n+1));  # Indeces of "end point" beads (= "measurement" points)

@everywhere (
    const s = 2;        # s = number of system parameters (k, gamma)
    using ForwardDiff;
)    

##
##
## ============================================================================================
## Parameters to be inferred:
## ============================================================================================
##
##
true_K     = 50.  ;                               # Retention time
true_gamma = 0.2  ;                               # Dimensionless noise parameter
sigma      = 0.05 ;                               # Measurement noise
true_theta = [log(true_K/T),log(true_gamma)]      # Parameters to be inferred (beta, tau)
##
##
## ============================================================================================
## Rain input:
## ============================================================================================
##
##
r  = Array(Float64,N);
r  = sin(t/100.).*sin(t/100.) + 0.1;              # Simple sinusoidal model

# r=float64(readdlm("$dir/r.dat")[range,2])

lnr_der = zeros(N)                                # Calculate logarithmic derivative of rain input
                                                  # Will be used below
for i=1:(N-1)                                     
    lnr_der[i] = ( log(r[i+1]) - log(r[i]) ) / dt     
end
##
##
## ============================================================================================
## System realization (Ito convention):
## ============================================================================================
##
##
for p = procs()                    
    @spawnat p srand(p*18072011)   # Seeding the RNG provides reproducibility
    # @spawnat p srand(p*time())   # Seeding with actual time provides randomness
end                                 
                                 
S = Array(Float64,N);
S[1] = true_K * r[1];   # -------------------->   # Unperturbed steady state (with constant rain input)
for i = 2 : N
    S[i] = S[i-1] + dt * ( r[i-1] - S[i-1]/true_K ) +
                sqrt(dt) * sqrt(true_gamma/true_K) * S[i-1] * randn()
end;
y  = max(0.01,(1/true_K)*S[ty]+sigma*randn(n+1))  # Generating measurement points
##
##
## ============================================================================================
## Hamiltonian Monte Carlo:
## ============================================================================================
##
require("$dir/ParInf_Fun_AD_3.jl")
##
## --------------------------------------------------------------------------------------------
## Local and distributed containers for HMC:
## --------------------------------------------------------------------------------------------
##
nsample_burnin = 0
nsample_eff    = 1
nsample        = nsample_eff + nsample_burnin
theta_sample   = Array(Float64,nsample+1,s)
u_sample       = Array(Float64,nsample+1,N)
u_sample_par   = dzeros((nsample,N), workers(), [1,nworkers()])
energies       = Array(Float64,nsample)

q              = Array(Float64,N)      # Initialization of a q array on local process (1)
u1             = Array(Float64,N)      # Initialization of a u array on local process (1)

mp_theta        = Array(Float64, s)     # Bead masses (parameters)
mp1            = Array(Float64, N)     # Bead masses (u coordinates)

p_theta        = Array(Float64, s)     # Array of momenta (only parameters)
p              = dzeros((N,), workers(), [nworkers()]) # DArray of momenta (u coordinates)
theta_save     = Array(Float64, s)
u_save         = dzeros((N,), workers(), [nworkers()]) 

# To store the output of the derivative of V_fast w.r.t. u
der_out        = dzeros((N,), workers(), [nworkers()])  
##
## --------------------------------------------------------------------------------------------
## MD parameters:
## --------------------------------------------------------------------------------------------
##
dtau    = 0.02                          # Long range MD time step
nn      = 16		                    # To define short-range time steps in RESPA
n_respa = 10                            # Number of large time steps
##
## --------------------------------------------------------------------------------------------
## Initial state:
## --------------------------------------------------------------------------------------------
##
K = 200.0
gamma = 0.5
theta = [log(K/T),log(gamma)]

S_init = Array(Float64,N);
for i = 1:n
    S_init[ty[i]:ty[i+1]] = K*linspace(y[i],y[i+1],ty[i+1]-ty[i]+1)
end

q = log(S_init./(K*r))

## Alternatively, q could be directly initialized as a DArray.
## ----------------------------------------------------------
## q = DArray((N,), workers(), [nworkers()]) do I
##     [log(S_init[i]/(K*r[i])) for i in I[1]]
## end
## ----------------------------------------------------------
## However, the body of the init function -do I ... end- would require
## every object used therein to be previously defined @everywhere

## Transformations q -> u 
## (eqs 2.16-17 in Tuckerman et al., JCP 99 (4), 2796, 1993)

for i=0:(n-1)
    u1[i*j+1] = q[i*j+1]
    for k=2:j
        u1[i*j+k] = q[i*j+k] - ( (k-1)*q[i*j+k+1] + q[i*j+1] )/k
    end
end
u1[N] = q[N]

theta_sample[1,:] = theta
u_sample[1,:]     = u1
##
## --------------------------------------------------------------------------------------------
## Masses (burn-in):
## --------------------------------------------------------------------------------------------
##
M_bdy   = 6.0*K*N/(gamma*T)
M_theta = M_bdy*0.02
M_stage = M_bdy*0.01

mp_theta[1:s] = M_theta

for i=0:(n-1)  
    mp1[i*j+1] = M_bdy                       # End point beads
    for k=2:j 
        mp1[i*j+k] = M_stage*k/(k-1)         # Staging beads
    end
end
mp1[N] = M_bdy                               # Last end point bead
##
## --------------------------------------------------------------------------------------------
## Going parallel:
## --------------------------------------------------------------------------------------------
##
## q  = distribute(q1)    # Distributed q array
u  = distribute(u1)    # Distributed u array
mp = distribute(mp1)   # Distributed array of bead masses
##
## --------------------------------------------------------------------------------------------
## Quick security check:
## --------------------------------------------------------------------------------------------
##
## if  (q.indexes) != (u.indexes) 
##     error("Error, the indexes of q and u are not distributed correctly !")
## end

if  (u.indexes) != [((u_sample_par.indexes)[1][2],),((u_sample_par.indexes)[2][2],)] 
    error("Error, the indexes of u and u_sample_par are not distributed correctly !")
end

if  (u.indexes) != (der_out.indexes) 
    error("Error, the indexes of u and der_out are not distributed correctly !")
end

if  (u.indexes) != (mp.indexes) 
    error("Error, the indexes of u and mp are not distributed correctly !")
end

if  (mp.indexes) != (p.indexes) 
    error("Error, the indexes of mp and p are not distributed correctly !")
end
##
## --------------------------------------------------------------------------------------------
## HMC loops
## --------------------------------------------------------------------------------------------
##
reject_counter_burnin = 0

println(string("\nStarting HMC loops (burn-in)...\n---------------------------------\n"))
tinit=time()

for counter = 1:nsample_burnin

    # Sample momenta and save current state:
    # -------------------

    p_theta = sqrt(mp_theta).*randn(s)           # on local process (1)
    theta_save = theta    

    for wp = workers()                          # parallel on worker processes 
        @spawnat wp (
            for jwp = 1:size(localpart(u),1) 
                localpart(p)[jwp]      = sqrt(localpart(mp)[jwp])*randn();
                localpart(u_save)[jwp] = localpart(u)[jwp];
                # in the serial version: for i = 1:N   u_save[i] = u[i] end
            end
        )
    end
    
    # Calculate energy:
    # -------------------

    H_old = V_fast(theta,u) + V_slow(theta,q) + sum((p_theta.^2) ./ (2*mp_theta)) +
                sum({fetch(@spawnat wp sum((localpart(p).^2)./(2*localpart(mp)))) for wp=workers()})

    energies[counter] = H_old

    # MD Integration:
    # -------------------

    for counter_respa = 1:n_respa 
        RESPA(theta,u,q,p,mp,dtau,nn) 
    end 

    # Calculate energy of proposal state:
    # -------------------

    H_new = V_fast(theta,u) + V_slow(theta,q) + sum((p_theta.^2) ./ (2*mp_theta)) +
                sum({fetch(@spawnat wp sum((localpart(p).^2)./(2*localpart(mp)))) for wp=workers()})

    # Metropolis step:
    # -------------------

    accept_prob = min(1,exp(H_old-H_new))
    if rand() > accept_prob       
        theta = theta_save              # on local process (1) 
        for wp = workers()              # parallel on worker processes 
            @spawnat wp (
                for jwp = 1:size(localpart(u),1) 
                    localpart(u)[jwp] = localpart(u_save)[jwp]
                end
            )
        end
        reject_counter_burnin += 1
    end
    
    theta_sample[counter+1,:] = theta   # on local process (1)
    for wp = workers()                  # parallel on worker processes 
        @spawnat wp (
            for jwp = 1:size(localpart(u),1) 
                localpart(u_sample_par)[counter,jwp] = localpart(u)[jwp]
            end
        )
    end
   
    println(string(counter, " loops completed in ", round(time()-tinit,1), " seconds \n"))   
    #if (counter%100 == 0)
    #    println(string(counter, " loops completed in ", round(time()-tinit,1), " seconds \n"))
    #end

end

## --------------------------------------------------------------------------------------------
## Redefinition of masses (effective HMC loops):
## --------------------------------------------------------------------------------------------

M_bdy_burnin   = M_bdy
M_theta_burnin = M_theta
M_stage_burnin = M_stage

M_bdy   = 6.0*K*N/(gamma*T)
M_theta = M_bdy*[0.03, 0.03]	          # Masses for K and gamma
M_stage = M_bdy*0.02

mp_theta = M_theta

mp1 = zeros(Float64,N)
for i=0:(n-1)  
    mp1[i*j+1] = M_bdy                    # End point beads
    for k=2:j 
        mp1[i*j+k] = M_stage*k/(k-1)      # Staging beads
    end
end
mp1[N] = M_bdy                            # Last end point bead

mp = distribute(mp1)

reject_counter = 0

# Initialization of vectors to store execution times
time_Hold       = Array(Float64,nsample_eff)
time_respa      = Array(Float64,nsample_eff)
time_Hnew       = Array(Float64,nsample_eff)
time_metropolis = Array(Float64,nsample_eff) 

time_respa_f        = Array(Float64,n_respa)
time_respa_f1       = Array(Float64,n_respa,nn)
time_respa_f2       = Array(Float64,n_respa,nn)
time_respa_f3       = Array(Float64,n_respa,nn)
time_respa_f4       = Array(Float64,n_respa,nn)

println(string("\nStarting effective HMC loops...\n---------------------------------\n"))

for counter = (nsample_burnin + 1):nsample

    t0=time()

    # println(string("Starting loop ", counter, " \n"))

    # Sample momenta and save current state:
    # -------------------

    p_theta = sqrt(mp_theta).*randn(s)          # on local process (1)
    theta_save = theta    

    for wp = workers()                          # parallel on worker processes 
        @spawnat wp (
            for jwp = 1:size(localpart(u),1) 
                localpart(p)[jwp]      = sqrt(localpart(mp)[jwp])*randn();
                localpart(u_save)[jwp] = localpart(u)[jwp];
                # in the serial version: for i = 1:N   u_save[i] = u[i] end
            end
        )
    end

    # println(string("Samplin momenta loop ", counter, " \n"))

    # Calculate energy:
    # -------------------

    H_old = V_fast(theta,u) + V_slow(theta,q) + sum((p_theta.^2) ./ (2*mp_theta)) +
                sum({fetch(@spawnat wp sum((localpart(p).^2)./(2*localpart(mp)))) for wp=workers()})

    energies[counter] = H_old
    
    time_Hold[counter-nsample_burnin] = time()-t0

    # println(string("H_old loop ", counter, " \n"))

    # MD Integration:
    # -------------------
    t1=time()
    
    for counter_respa = 1:n_respa 
        RESPA(theta,u,q,p,mp,dtau,nn,counter_respa,time_respa_f,time_respa_f1,time_respa_f2,time_respa_f3,time_respa_f4,der_out) 
    end 

    time_respa[counter-nsample_burnin] = time()-t1

    # println(string("RESPA loop", counter, " \n"))

    # Calculate energy of proposal state:
    # -------------------

    t2=time()

    H_new = V_fast(theta,u) + V_slow(theta,q) + sum((p_theta.^2) ./ (2*mp_theta)) +
                sum({fetch(@spawnat wp sum((localpart(p).^2)./(2*localpart(mp)))) for wp=workers()})

    time_Hnew[counter-nsample_burnin] = time()-t2

    # println(string("H_new loop", counter, " \n"))

    # Metropolis step:
    # -------------------

    t3=time()

    accept_prob = min(1,exp(H_old-H_new))
    if rand() > accept_prob
        theta = theta_save              # on local process (1) 
        for wp = workers()              # parallel on worker processes 
            @spawnat wp (
                for jwp = 1:size(localpart(u),1) 
                    localpart(u)[jwp] = localpart(u_save)[jwp]
                end
            )
        end
        reject_counter += 1
    end

    theta_sample[counter+1,:] = theta   # on local process (1)
    for wp = workers()                  # parallel on worker processes 
        @spawnat wp (
            for jwp = 1:size(localpart(u),1) 
                localpart(u_sample_par)[counter,jwp] = localpart(u)[jwp]
            end
        )
    end

    time_metropolis[counter-nsample_burnin] = time()-t3

    # println(string("Metropolis and end loop", counter, " \n"))

    # println(string(counter, " loops completed in ", round(time()-tinit,1), " seconds \n"))
    if (counter%10 == 0)
        println(string(counter, " loops completed in ", round(time()-tinit,1), " seconds \n"))
    end

end
##
## --------------------------------------------------------------------------------------------
## End of HMC loop
## --------------------------------------------------------------------------------------------
##
## 
## --------------------------------------------------------------------------------------------
## Back transformations (u -> q):
## (Tuckerman et al., JCP 99 (1993), eq. 2.19)
## --------------------------------------------------------------------------------------------
##
u_sample[2:end,:] = convert(Array,u_sample_par)

qs = Array(Float64, nsample+1, N)
for index = 1:(nsample+1)
    ss=n-1
    while ss>=0
        k=j+1
        qs[index,ss*j+k] = u_sample[index,ss*j+k]
        while k>2
           k -= 1
           qs[index,ss*j+k] = u_sample[index,ss*j+k] +
                              u_sample[index,ss*j+1]/k +
                              (k-1)*qs[index,ss*j+k+1]/k
         end
        ss -= 1
    end
    qs[index,1] = u_sample[index,1]
end

predy = Array(Float64,nsample+1,N)                              # Expected outputs
for index = 1:(nsample+1)
    for l=1:N
        predy[index,l] = r[l]*exp(qs[index,l])
    end
end

println(string("\nRun completed in ", round(time()-tinit,2), " seconds \n"))
println(string("Calculation of H_old in ", round(sum(time_Hold),2), " seconds \n"))
println(string("RESPA cycles in ", round(sum(time_respa),2), " seconds \n"))
println(string("Calculation of H_new in ", round(sum(time_Hnew),2), " seconds \n"))
println(string("Metropolis steps in ", round(sum(time_metropolis),2), " seconds \n"))

##
##
## ============================================================================================
## Saving parameters and results
## ============================================================================================
##
##

param_names  = vcat("N", "j", "n", "t[1]", "dt", "s", "true_K", "true_gamma", "sigma", "K", "gamma", 
	"nsample_burnin", "nsample_eff", "M_bdy_burnin", "M_bdy", "M_theta_burnin", "M_theta_K", 
	"M_theta_gamma", "M_stage_burnin", "M_stage", "dtau", "nn", "n_respa")
param_values = vcat(N, j, n, t[1], dt, s, true_K, true_gamma, sigma, K, gamma, 
	nsample_burnin, nsample_eff, M_bdy_burnin, M_bdy, M_theta_burnin, M_theta, 
	M_stage_burnin, M_stage, dtau, nn, n_respa)
writedlm("$dir2/params$fname.dat", hcat(param_names, param_values))

last_qs = qs[1:(nsample+1),N]          
writedlm("$dir2/last_qs$fname.dat", last_qs)

red_range = 200; red_predy = Array(Float64, red_range, N); row_ind = 1
for ip in iround(linspace(1,nsample+1,min(size(predy,1),200)))
	red_predy[row_ind,1:N] = predy[ip,1:N]
	row_ind += 1
end
writedlm("$dir2/predy$fname.dat", red_predy)

thetas = Array(Float64,nsample+1,s)
thetas[:,1] = exp(theta_sample[:,1])*T                  # K
thetas[:,2] = exp(theta_sample[:,2])                    # gamma
writedlm("$dir2/thetas$fname.dat", thetas)

writedlm("$dir2/energies$fname.dat", energies)

reject_rates = [reject_counter_burnin/nsample_burnin*100,reject_counter/nsample_eff*100]
writedlm("$dir2/reject$fname.dat", reject_rates)

S_first  = Array(Any, N+1); S_first  = vcat("System_realization", S)
rain_in  = Array(Any, N+1); rain_in  = vcat("Rain_input", r)
flow_out = Array(Any, n+2); flow_out = vcat("Output_flow", y)
io_data  = Array(Any, 2*N+n+4); io_data  = vcat(S_first, rain_in, flow_out)
writedlm("$dir2/iodata$fname.dat", io_data)

##
##
## ============================================================================================
## Close active worker processes
## ============================================================================================
##
##

for i = 2:nprocs()
    rmprocs(procs()[i])
end

