## ============================================================================================
## Project: Parameter Inference for a Simple Stochastic Hydrological Model
## File: Main_HMC_1.jl
##
## Description: Single bucket, scale invariant noise, no evaporation yet
##              dS(t)/dt = r(t) - S(t)/k + sqrt(gamma/k) S(t) eta(t)
## The {q} variable set is used for slow (long-range) components. (NOT ANY MORE TRUE!)
## The {u} variable set is used for fast (short-range) components.
## 
## The derivative of the potential V(theta,u/q) w.r.t. theta can be calculated explicitely
## or by automated differentiation (ForwardDiff package). 
## 
## Different masses are used for the burn-in phase.
##
## NEW IN VERSION 2. THE PRIOR FOR PARAMETER K IS A LOG-NORMAL DISTRIBUTION.
##                   We use the tentative function f(K) = (3 Exp(-(Log(K)-4)^2/2))/K 
## NEW IN VERSION 3. We use the reverse diff package 
##             
## GitHub repository: https://github.com/ulzegasi/Julia_ParInf.git
##
## Authors: Simone Ulzega, Carlo Albert
##
## simone.ulzega@eawag.ch
## carlo.albert@eawag.ch
## ============================================================================================
dir   = "C:/Users/ulzegasi/Julia_files/ParInf_HMC"       # Main project directory, local repository   
dir2  = "C:/Users/ulzegasi/SWITCHdrive/JuliaTemp/Data"   # Secondary directory (e.g., data storage)
fname = string("_testing")                               # Output files
using ReverseDiffSource, ForwardDiff                     # Automated differentiation package
##
## ============================================================================================
## System parameters
## ============================================================================================
##
range = 2:5002
t  = float64(readdlm("$dir/t.dat")[range,2]) # Time points t.dat

const n = 10             # n+1 -> "end point" beads (see Tuckerman '93), n -> number of segments
const j = 30             # n(j-1) -> total number of staging beads, j-1 -> staging beads per segment
const N = int64(n*j+1)   # Total number of discrete time points
                         # IMPORTANT: (N-1)/j = integer = n (measurement points)
if  ((N-1)%j) != 0 
    error("Be careful, the number of staging points j must fulfill (N-1)/j = integer !")
end

t = linspace(t[1], t[end], N)                 # Time points      
T  = t[end]-t[1]                              # Total time interval
dt = T/(N-1)                                  # Time step
ty  = iround(linspace(1, N, n+1))             # Indeces of "end point" beads (= "measurement" points)

const s = 2                                   # Number of system parameters (k, gamma)     

nsample_burnin  = 0                           # Number of system realizations (burnin and effective) 
nsample_eff     = 100
nsample         = nsample_eff + nsample_burnin

dtau    = 0.02                                # Long range MD time step
nn      = 16  ; deltatau = dtau / nn          # Short-range time steps in RESPA
n_respa = 10                                  # Number of large RESPA time steps

true_K     = 50.                              # Retention time
true_gamma = 0.2                              # Dimensionless noise parameter
sigma      = 0.05                             # Measurement noise
true_theta = [log(true_K/T),log(true_gamma)]  # Parameters to be inferred (beta, tau)

K     = 200.0                                 # Initial state
gamma = 0.5                                   
theta = [log(K/T),log(gamma)]
##                 
## ============================================================================================
## Initialization of vectors and containers
## ============================================================================================
##
q            = zeros(N)             # Coordinates vector
u            = zeros(N)             # Coordinates vector
r            = zeros(N)             # Rain input vector
y            = zeros(n+1)           # Measurements
lnr_der      = zeros(N)             # Logarithmic derivative of rain input

S            = zeros(N)             # System realization (with true parameters, used to create measurements y)
S_init       = zeros(N)             # Initial state (with starting parameter values)

theta_sample = zeros(nsample+1, s)  # Containers for sampled parameters
u_sample     = zeros(nsample+1, N)  # Containers for sampled coordinates
energies     = zeros(nsample)       # Containers for sampled energies
mp           = zeros(s+N)           # Masses for the p momenta
p            = zeros(s+N)           # Momenta
theta_save   = zeros(s)             # To store temporary parameters
u_save       = zeros(N)             # To store temporary coordinates

##
## ============================================================================================
## Rain input:
## ============================================================================================
##
r = sin(t/100.).*sin(t/100.) + 0.1;      # Simple sinusoidal model
## r=float64(readdlm("$dir/r.dat")[range,2])                                                   

for i = 1:(N-1)                          # Calculate logarithmic derivative of rain input                
    lnr_der[i] = ( log(r[i+1]) - log(r[i]) ) / dt
end
##
## ============================================================================================
## System realization (Ito convention) and measurement points:
## ============================================================================================
##
srand(18072011)                 # Seeding the RNG provides reproducibility
#srand(time())                  # Seeding the RNG with actual time provides randomness

S[1] = true_K * r[1]            # Unperturbed steady state (with constant rain input)
for i = 2 : N
    S[i] = S[i-1] + dt * ( r[i-1] - S[i-1]/true_K ) +
                   sqrt(dt) * sqrt(true_gamma/true_K) * S[i-1] * randn()
end
y  = max(0.01,(1/true_K)*S[ty] + sigma*randn(n+1))      # Generation of measurement points
##
## ============================================================================================
## Hamiltonian Monte Carlo:
## ============================================================================================

## --------------------------------------------------------------------------------------------
## Initial state:
## --------------------------------------------------------------------------------------------

for a = 1:n
    S_init[ty[a]:ty[a+1]] = K*linspace(y[a],y[a+1],ty[a+1]-ty[a]+1)
end

q = log(S_init./(K*r))

## Transformations q -> u 
## (eqs 2.16-17 in Tuckerman et al., JCP 99 (4), 2796, 1993)

for i=0:(n-1)
    u[i*j+1] = q[i*j+1]
    for k=2:j
        u[i*j+k] = q[i*j+k] - ( (k-1)*q[i*j+k+1] + q[i*j+1] )/k
    end
end
u[N] = q[N]

theta_sample[1,:] = theta          # Store initial parameter values
u_sample[1,:]     = u              # Store initial coordinates

println(string("\nLoading potentials, derivatives and functions...\n---------------------------------"))
reload("$dir/ParInf_Fun_AD_3_2.jl")

## --------------------------------------------------------------------------------------------
## Masses (burn-in):
## --------------------------------------------------------------------------------------------

M_bdy   = 6.0*K*N/(gamma*T)
M_theta = M_bdy*0.02
M_stage = M_bdy*0.01

mp[1:s] = M_theta
for i=0:(n-1)  
    mp[s+1+i*j] = M_bdy                                 # End point beads
    for k=2:j 
        mp[s+i*j+k] = M_stage*k/(k-1)                   # Staging beads
    end
end
mp[s+N] = M_bdy                                         # Last end point bead

## --------------------------------------------------------------------------------------------
## HMC loops
## --------------------------------------------------------------------------------------------

reject_counter_burnin = 0

println(string("\nStarting HMC loops (burn-in)...\n---------------------------------\n"))
tinit=time()

for counter = 1:nsample_burnin

    # Sample momenta:
    # -------------------

    p = sqrt(mp).*randn(s+N)
    
    # Calculate energy:
    # -------------------

    H_old = V_fast_fun(theta,u) + V_slow_fun(theta,u) + sum((p .* p) ./ (2*mp))
    energies[counter] = H_old
    
    # Save current state:
    # -------------------

    for i = 1:s   theta_save[i] = theta[i]    end
    for i = 1:N   u_save[i]     = u[i]        end

    # MD Integration:
    # -------------------

    for counter_respa = 1:n_respa 
        RESPA(theta, u, p, mp, dtau, nn, deltatau) 
    end 

    # Calculate energy of proposal state:
    # -------------------

    H_new = V_fast_fun(theta,u) + V_slow_fun(theta,u) + sum(p .* p ./ (2*mp))

    # Metropolis step:
    # -------------------

    accept_prob = min(1,exp(H_old-H_new))
    if rand() > accept_prob
        for i=1:s theta[i] = theta_save[i] end
        for i=1:N u[i]     = u_save[i] end
        reject_counter_burnin += 1
    end
    
    theta_sample[counter+1,:] = theta
    u_sample[counter+1,:] = u 
   
    if (counter%100 == 0)
        println(string(counter, " loops completed in ", round(time()-tinit,1), " seconds \n"))
    end

end

## --------------------------------------------------------------------------------------------
## Redefinition of masses (effective HMC loops):
## --------------------------------------------------------------------------------------------

M_bdy_burnin   = M_bdy
M_theta_burnin = M_theta
M_stage_burnin = M_stage

M_bdy   = 6.0*K*N/(gamma*T)
M_theta = M_bdy*[0.03, 0.03]						    # Masses for K and gamma
M_stage = M_bdy*0.02

mp = zeros(Float64, s+N)
mp[1:s] = M_theta
for i=0:(n-1)  
    mp[s+1+i*j] = M_bdy                                 # End point beads
    for k=2:j 
        mp[s+i*j+k] = M_stage*k/(k-1)                   # Staging beads
    end
end
mp[s+N] = M_bdy                                         # Last end point bead

reject_counter = 0

# Initialization of vectors to store execution times

time_Hold       = Array(Float64,nsample_eff)
time_respa      = Array(Float64,nsample_eff)
time_Hnew       = Array(Float64,nsample_eff)
time_metropolis = Array(Float64,nsample_eff) 

println(string("\nStarting effective HMC loops...\n---------------------------------\n"))

for counter = (nsample_burnin + 1):nsample

    # Sample momenta:
    # -------------------
    t0=time()
    
    p = sqrt(mp).*randn(s+N)

    # Calculate energy:
    # -------------------

    H_old = V_fast_fun(theta,u) + V_slow_fun(theta,u) + sum((p .* p) ./ (2*mp))
    energies[counter] = H_old
    
    time_Hold[counter-nsample_burnin] = time()-t0
    
    # Save current state:
    # -------------------

    for i = 1:s   theta_save[i] = theta[i]    end
    for i = 1:N   u_save[i]     = u[i]        end

    # MD Integration:
    # -------------------
    t1=time()
    
    for counter_respa = 1:n_respa 
        RESPA(theta, u, p, mp, dtau, nn, deltatau)  
    end 
    time_respa[counter-nsample_burnin] = time()-t1

    # Calculate energy of proposal state:
    # -------------------

    t2=time()

    H_new = V_fast_fun(theta,u) + V_slow_fun(theta,u) + sum(p .* p ./ (2*mp))

    time_Hnew[counter-nsample_burnin] = time()-t2

    # Metropolis step:
    # -------------------

    t3=time()

    accept_prob = min(1,exp(H_old-H_new))
    if rand() > accept_prob
        for i=1:s theta[i] = theta_save[i] end
        for i=1:N u[i]     = u_save[i] end
        reject_counter += 1
    end
    
    theta_sample[counter+1,:] = theta
    u_sample[counter+1,:] = u 

    time_metropolis[counter-nsample_burnin] = time()-t3
    if (counter%100 == 0)
        println(string(counter, " loops completed in ", round(time()-tinit,1), " seconds \n"))
    end

end

## --------------------------------------------------------------------------------------------
## End of HMC loop
## --------------------------------------------------------------------------------------------

## --------------------------------------------------------------------------------------------
## Back transformations (u -> q):
## (Tuckerman et al., JCP 99 (1993), eq. 2.19)
## --------------------------------------------------------------------------------------------

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

println(string("\nRun completed in ", round(time()-tinit,6), " seconds \n"))
println(string("Calculation of H_old in ", round(sum(time_Hold),6), " seconds \n"))
println(string("RESPA cycles in ", round(sum(time_respa),6), " seconds \n"))
println(string("Calculation of H_new in ", round(sum(time_Hnew),6), " seconds \n"))
println(string("Metropolis steps in ", round(sum(time_metropolis),6), " seconds \n"))

##
## ============================================================================================
## Saving parameters and results
## ============================================================================================
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
