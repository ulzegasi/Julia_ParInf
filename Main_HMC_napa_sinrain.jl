## ============================================================================================
## Project: Parameter Inference for a Simple Stochastic Hydrological Model
##
## Description: Single bucket, scale invariant noise
##              dS(t)/dt = r(t) - S(t)/k + sqrt(gamma/k) S(t) eta(t) 
##
## NAPA version. The fast time evolution of the staging beads is solved analytically (innermost 
##               cycle of the old RESPA version). All other degrees of freedom are propagated 
##               by Velocity Verlet Integrators.   
##
## GitHub repository: https://github.com/ulzegasi/Julia_ParInf.git
##
## Authors: Simone Ulzega, Carlo Albert
##
## simone.ulzega@eawag.ch
## carlo.albert@eawag.ch
## ============================================================================================
dir   = "C:/Users/ulzegasi/Julia_files/ParInf_HMC"  # Main project directory, local repository   
dir2  = "$dir/temp_data"                            # Secondary directory (e.g., output storage)
dir3  = "$dir/input_data"                           # Input data directory
fname = string("_n10_K50_g02_s10_sinr")                        # Output files
using ReverseDiffSource, PyPlot, Distributions #, ForwardDiff      # Automated differentiation package
pygui(true)
srand(18072011)                                     # Seeding the RNG provides reproducibility
# srand(time())                                     # Seeding the RNG with actual time provides randomness
##
## ============================================================================================
## Load data
## ============================================================================================
##
St     = readdlm("$dir3/St$fname.dat")
S      = St[:,1]
long_t = St[:,2]   
y      = vec(readdlm("$dir3/y$fname.dat"))
##
## ============================================================================================
## System parameters
## ============================================================================================
##
range = 2:5002
tdat  = float64(readdlm("$dir/t.dat")[range,2]) # Time points t.dat

const nparams = 2              # Number of system parameters (k, gamma)
const n = 10                  # n+1 -> "end point" beads (see Tuckerman '93), n -> number of segments
const j = 30                  # n(j-1) -> total number of staging beads, j-1 -> staging beads per segment
const N = int64(n*j+1)        # Total number of discrete time points = n*j + 1
                              # IMPORTANT: (N-1)/j = integer = n (measurement points)
if  ((N-1)%j) != 0 
    error("Be careful, the number of staging points j must fulfill (N-1)/j = integer !")
end

const t  = linspace(tdat[1], tdat[end], N) # Generate N time points in [t[1],t[end]]
const T  = t[end]-t[1]                     # Total time interval
const dt = T/(N-1)                         # Time step
const ty = iround(linspace(1, N, n+1))     # Indeces of "boundary" beads (= measurement points)     

const nsample_burnin = 0                   # Number of points in the MCMC
const nsample_eff    = 20000
const nsample        = nsample_eff + nsample_burnin

const dtau           = 0.25                # MD time step
const n_napa         = 3                   # Number of NAPA time steps

const true_K     = 50.                    # Retention time
const true_gam   = 0.2                     # Dimensionless noise parameter
const true_bet   = sqrt(T*true_gam/true_K)
const sigma      = 0.10                    # Measurement noise
const true_theta = [true_bet,true_gam]     # Parameters to be inferred (beta, tau)

K     = 200.0                              # Initial state
gam   = 0.5
bet   = sqrt(T*gam/K)                                   
theta = [bet, gam]

# Logistic prior parameters
# const s_gam = 4
# const gam_0 = 3
# const s_k   = 0.04
# const k_0   = 500

# Parameter limits
K_min   = 0.0
gam_min = 0.0
K_max   = 1000.0
gam_max = 5.0
##                 
## ============================================================================================
## Initialization of arrays
## ============================================================================================
##
## Data
lnr_der = Array(Float64,N)      # Log derivative of (smoothed) rain input
## System coordinates/modes
q       = Array(Float64,N)       
u       = Array(Float64,N)       
## Masses and momenta for HMC
mp      = Array(Float64,nparams+N)     
p       = Array(Float64,nparams+N)
## Containers for sampled objects
theta_sample = Array(Float64,nsample+1,nparams) # Container for sampled parameters
u_sample     = Array(Float64,nsample+1,N)      # Container for sampled coordinates
energies     = Array(Float64,nsample)          # Container for sampled energies
## To store temporary values
theta_save   = Array(Float64,nparams)           # To store temporary parameters
u_save       = Array(Float64,N)                # To store temporary coordinates
## results
qs           = Array(Float64,nsample+1,N)      # Sampled coordinates
predy        = Array(Float64,nsample+1,N)      # Sampled outputs
##
## ============================================================================================
## Synthetic data
## ============================================================================================
##
#=
S[1] = true_K * r[1]  # Define first point (--> <S*Peq(S)> with constant rain input)
for i = 2:N              
    S[i] = S[i-1] + dt*( r[i-1] - S[i-1]/true_K ) + sqrt(dt)*sqrt(true_gam/true_K)*S[i-1]*randn()
    if  S[i] < 0 
        error("Negative volume!")
    end
end
=#
##
## ============================================================================================
## Data
## ============================================================================================
##
r = Array(Float64,N); r = sin(t/100.).*sin(t/100.) + 0.1
for i = 1:(N-1)  
	lnr_der[i] = (log(r[i+1])-log(r[i]))/dt   # Log-derivative of the rain
end

q_init = (1/true_bet)*log(y./r[ty])
bq     = true_bet*q_init
for s = 1:n 								  # ... and their linear interpolation (staging points)
    q[ty[s]:ty[s+1]] = linspace(q_init[s],q_init[s+1],ty[s+1]-ty[s]+1)
end

############## It could be interesting to plot a few figures to visualize input data ##############
#=
plt.figure(figsize=(12.5, 4.5))
axes()[:set_ylim]([0,2])
axes()[:set_xlim]([-5,840])
plt.xlabel("time")
plt.ylabel("r")
plt.plot(t, r, "bo", markersize = 6)
plt.plot(t, r, "b", linewidth=1)

plt.figure(figsize=(12.5, 4.5))
axes()[:set_ylim]([0,2])
axes()[:set_xlim]([-5,840])
plt.xlabel("time")
plt.ylabel("S/K, y")
plt.plot(long_t, S/true_K, "r", linewidth = 2, color = (1,0.65,0))
yerr=2*sigma*y
plt.errorbar(t[ty], y, yerr=(yerr,yerr), fmt="o", color = (1,0,0), markersize = 10, capsize=8, elinewidth=4)

plt.figure(figsize=(12.5, 4.5))
axes()[:set_ylim]([0,2])
axes()[:set_xlim]([-5,840])
plt.xlabel("time")
plt.ylabel("S/K, y")
plt.plot(long_t, S/true_K, "r", linewidth = 2, color = (1,0.65,0))
plt.plot(t, r, "bo", markersize = 6)
plt.plot(t, r, "b", linewidth=1)
yerr=2*sigma*y
plt.errorbar(t[ty], y, yerr=(yerr,yerr), fmt="o", color = (1,0,0), markersize = 10, capsize=8, elinewidth=4)

plt.figure(figsize=(12.5, 4.5))
axes()[:set_ylim]([-1,1])
axes()[:set_xlim]([-5,840])
plt.plot(t, q, "g-",linewidth=1)
plt.plot(t, q, "go",markersize=4)
plt.plot(t[ty], q[ty], "go",markersize=12)
=#                       

##
##
## ============================================================================================
## Hamiltonian Monte Carlo
## ============================================================================================
##
## Transformations q -> u 
## (eqs 2.16-17 in Tuckerman et al., JCP 99 (4), 2796, 1993)
##
for s = 0:(n-1)
    u[s*j+1] = q[s*j+1]
    for k = 2:j
        u[s*j+k] = q[s*j+k] - ( (k-1)*q[s*j+k+1] + q[s*j+1] )/k
    end
end
u[N] = q[N]

theta_sample[1,:] = theta          # Store initial parameters
u_sample[1,:]     = u              # Store initial coordinates

## --------------------------------------------------------------------------------------------
## Masses (burn-in) and frequencies:
## --------------------------------------------------------------------------------------------

m_bdy    = 10.0                       # m = m_q / dt
m_stg    = 1.0                        # we assume m_q prop. to dt ==> m = costant     
m_theta  = 1.0

mp[1:nparams] = m_theta
for s = 1:n  
    mp[nparams+(s-1)*j+1] = m_bdy                       
    for k = 2:j 
        mp[nparams+(s-1)*j+k] = m_stg                   
    end
end
mp[nparams+N] = m_bdy                                       

## --------------------------------------------------------------------------------------------
## Loading potentials, derivatives and Napa 
## --------------------------------------------------------------------------------------------

println(string("\nLoading potentials, derivatives and napa...\n---------------------------"))
reload("$dir/ParInf_Fun_AD_napa.jl")

## --------------------------------------------------------------------------------------------
## HMC loops
## --------------------------------------------------------------------------------------------

reject_counter_burnin = 0

println(string("\nStarting HMC loops (burn-in)...\n---------------------------------\n"))
tinit=time()

for counter = 1:nsample_burnin

    # Sample momenta:
    # -------------------
    p = sqrt(mp).*randn(nparams+N)
    
    # Calculate energy:
    # -------------------
    H_old = V_N_fun(u) + V_n_fun(theta,u) + V_1_fun(theta,u) + sum((p.^2)./(2*mp))
    energies[counter] = H_old
    if  isnan(H_old) 
        error(string("Iteration ", counter, " --> energy values diverged..."))
    end
    
    # Save current state:
    # -------------------
    for i = 1:nparams theta_save[i] = theta[i] end
    for i = 1:N      u_save[i]     = u[i]     end

    # MD Integration:
    # -------------------
    for counter_napa = 1:n_napa 
        napa(theta, u, counter) 
    end 

    # Metropolis step and energy of proposal state
    # -------------------
    if (theta[2] <= gam_min) || (theta[2] >= gam_max) || ((T*theta[2]/(theta[1])^2) <= K_min) || ((T*theta[2]/(theta[1])^2) >= K_max) 
        for i = 1:nparams theta[i] = theta_save[i] end
        for i = 1:N      u[i]     = u_save[i]     end
        reject_counter_burnin += 1
    elseif any(u .> 10.0) || any(u .< -10.0)
        for i = 1:nparams theta[i] = theta_save[i] end
        for i = 1:N      u[i]     = u_save[i]     end
        reject_counter_burnin += 1
    else
        H_new = V_N_fun(u) + V_n_fun(theta,u) + V_1_fun(theta,u) + sum((p.^2)./(2*mp))
        accept_prob = min(1,exp(H_old-H_new))
        if rand() > accept_prob 
            for i = 1:nparams theta[i] = theta_save[i] end
            for i = 1:N      u[i]     = u_save[i]     end
            reject_counter_burnin += 1
        end
    end
    
    theta_sample[counter+1,:] = theta
    u_sample[counter+1,:]     = u 
   
    if (counter%100 == 0)
        println(string(counter, " loops completed in ", round(time()-tinit,1), " seconds \n"))
    end

end

## --------------------------------------------------------------------------------------------
## Redefinition of masses (effective HMC loops):
## --------------------------------------------------------------------------------------------

m_bdy_burnin   = m_bdy
m_stg_burnin   = m_stg
m_theta_burnin = m_theta

m_bdy    = 800                      # m = m_q / dt
m_stg    = 120                      # we assume m_q prop. to dt ==> m = costant     
m_theta  = [130, 130]

mp[1:nparams] = m_theta
for s = 1:n  
    mp[nparams+(s-1)*j+1] = m_bdy                       
    for k = 2:j 
        mp[nparams+(s-1)*j+k] = m_stg                   
    end
end
mp[nparams+N] = m_bdy          

reject_counter = 0

# Initialization of vectors to store execution times

time_respa      = Array(Float64,nsample_eff)
time_respa_s    = zeros(nsample_eff)
time_respa_f    = zeros(nsample_eff)
time_respa_f    = zeros(nsample_eff)

println(string("\nStarting effective HMC loops...\n---------------------------------\n"))

for counter = (nsample_burnin + 1):nsample

    # Sample momenta:
    # -------------------
    t0=time()
    p = sqrt(mp).*randn(nparams+N)

    # Calculate energy:
    # -------------------
    H_old = V_N_fun(u) + V_n_fun(theta,u) + V_1_fun(theta,u) + sum((p.^2)./(2*mp))
    energies[counter] = H_old
    if  isnan(H_old) 
        error(string("Iteration ", counter, " --> energy values diverged..."))
    end
    
    # Save current state:
    # -------------------
    for i = 1:nparams theta_save[i] = theta[i] end
    for i = 1:N      u_save[i]     = u[i]     end

    # MD Integration:
    # -------------------
    t1=time()
    for counter_napa = 1:n_napa 
        napa(theta, u, counter) 
    end 
    time_respa[counter-nsample_burnin] = time()-t1

    # Calculate energy of proposal state:
    # -------------------
    
    # println(H_new-H_old)

    # Metropolis step:
    # -------------------
    
    if (theta[2] <= gam_min) || (theta[2] >= gam_max) || ((T*theta[2]/(theta[1])^2) <= K_min) || ((T*theta[2]/(theta[1])^2) >= K_max) 
        for i = 1:nparams theta[i] = theta_save[i] end
        for i = 1:N      u[i]     = u_save[i]     end
        reject_counter += 1
    elseif any(u .> 10.0) || any(u .< -10.0)
        for i = 1:nparams theta[i] = theta_save[i] end
        for i = 1:N      u[i]     = u_save[i]     end
        reject_counter += 1
    else
        H_new = V_N_fun(u) + V_n_fun(theta,u) + V_1_fun(theta,u) + sum((p.^2)./(2*mp))
        accept_prob = min(1,exp(H_old-H_new))
        if rand() > accept_prob 
            for i = 1:nparams theta[i] = theta_save[i] end
            for i = 1:N      u[i]     = u_save[i]     end
            reject_counter += 1
        end
    end
    
    theta_sample[counter+1,:] = theta
    u_sample[counter+1,:]     = u 
   
    if (counter%100 == 0)
        println(string(counter, " loops completed in ", round(time()-tinit,1), " seconds \n"))
    end

end

## --------------------------------------------------------------------------------------------
## End of HMC loop
## --------------------------------------------------------------------------------------------
println(string("\nRun completed in ", round(time()-tinit,6), " seconds \n"))
println(string("RESPA cycles in ", round(sum(time_respa),6), " seconds \n"))
println(string("Slow RESPA in ", round(sum(time_respa_s),6), " seconds \n"))
println(string("Fast RESPA in ", round(sum(time_respa_f),6), " seconds \n"))
## --------------------------------------------------------------------------------------------
## Back transformations (u -> q):
## (Tuckerman et al., JCP 99 (1993), eq. 2.19)
## --------------------------------------------------------------------------------------------

for sample_ind = 1:(nsample+1)
    for s = 1:n
        qs[sample_ind, (s-1)*j+1] = u_sample[sample_ind, (s-1)*j+1]
        for k = 2:j 
            qs[sample_ind, (s-1)*j+k] = (j-k+1)*u_sample[sample_ind, (s-1)*j+1]/j
            for l = k:(j+1)
                qs[sample_ind, (s-1)*j+k] += (k-1)*u_sample[sample_ind, (s-1)*j+l]/(l-1)
            end
        end
    end
    qs[sample_ind, N] = u_sample[sample_ind, N]
end
# q = (1/bet)*log(S/(K*r))

for sample_ind = 1:(nsample+1)
    for ind = 1:N
        predy[sample_ind, ind] = r[ind]*exp(theta_sample[sample_ind, 1]*qs[sample_ind, ind])
    end
end
# y = r*exp(bet*q)

##
## ============================================================================================
## Save parameters and results
## ============================================================================================
##
fname = string("_sinr800_120") 

param_names  = vcat("N", "j", "n", "t[1]", "dt", "nparams", "true_K", "true_gam", "sigma", "K", "gam", 
	"nsample_burnin", "nsample_eff", "m_bdy_burnin", "m_bdy", "m_theta_burnin", "m_theta_bet", 
	"m_theta_gam", "m_stg_burnin", "m_stg", "dtau", "n_napa")
param_values = vcat(N, j, n, t[1], dt, nparams, true_K, true_gam, sigma, K, gam, 
	nsample_burnin, nsample_eff, m_bdy_burnin, m_bdy, m_theta_burnin, m_theta, 
	m_stg_burnin, m_stg, dtau, n_napa)
writedlm("$dir2/params$fname.dat", hcat(param_names, param_values))

last_qs = qs[1:(nsample+1),N]          
writedlm("$dir2/last_qs$fname.dat", last_qs)

u_chains          = Array(Float64, nsample+2, 4)
u_chains[1,1]     = ty[int64(floor(size(ty,1)/2))]
u_chains[2:end,1] = u_sample[:, ty[int64(floor(size(ty,1)/2))]]
u_chains[1,2]     = ty[int64(floor(size(ty,1)/2))] + int64(floor((j-1)/2))
u_chains[2:end,2] = u_sample[:, ty[int64(floor(size(ty,1)/2))] + int64(floor((j-1)/2))]
u_chains[1,3]     = ty[int64(floor(size(ty,1)/2))] + j
u_chains[2:end,3] = u_sample[:, ty[int64(floor(size(ty,1)/2))] + j]
u_chains[1,4]     = ty[int64(floor(size(ty,1)/2))] + j + int64(floor((j-1)/2))
u_chains[2:end,4] = u_sample[:, ty[int64(floor(size(ty,1)/2))] + j + int64(floor((j-1)/2))]
writedlm("$dir2/u_chains$fname.dat", u_chains)

#=
red_range = 200; red_predy = Array(Float64, red_range, N); row_ind = 1
if red_range < (nsample+1)
    for ip in iround(linspace(1, nsample+1, red_range))
        red_predy[row_ind,1:N] = predy[ip,1:N]
        row_ind += 1
    end
else
    red_predy = predy
end
=#

red_range = min(200,nsample+1); red_predy = Array(Float64, red_range, N); row_ind = 1
for ip in iround(linspace(1, nsample+1, red_range))
    red_predy[row_ind,1:N] = predy[ip,1:N]
    row_ind += 1
end
writedlm("$dir2/predy$fname.dat", red_predy)

red_qs = Array(Float64, red_range, N); row_ind = 1
for ip in iround(linspace(1, nsample+1, red_range))
    red_qs[row_ind,1:N] = qs[ip,1:N]
    row_ind += 1
end
writedlm("$dir2/predq$fname.dat", red_qs)

thetas      = Array(Float64,nsample+1,nparams)
thetas[:,1] = T*theta_sample[:,2]./(theta_sample[:,1]).^2  # K
thetas[:,2] = theta_sample[:,2]                            # gamma
writedlm("$dir2/thetas$fname.dat", thetas)

writedlm("$dir2/energies$fname.dat", energies)

reject_rates = [reject_counter_burnin/nsample_burnin*100,reject_counter/nsample_eff*100]
writedlm("$dir2/reject$fname.dat", reject_rates)

S_red = S[iround(linspace(1, size(S,1), N))]


S_first  = Array(Any, N+1); S_first  = vcat("System_realization", S_red)
#=S_first = vcat("System_realization", zeros(N)) =# 
rain_in  = Array(Any, N+1); rain_in  = vcat("Rain_input", r)
flow_out = Array(Any, n+2); flow_out = vcat("Output_flow", y)
io_data  = Array(Any, 2*N+n+4); io_data  = vcat(S_first, rain_in, flow_out)
writedlm("$dir2/iodata$fname.dat", io_data)