## ============================================================================================
## Project: Parameter Inference for a Simple Stochastic Hydrological Model
## Real data - San Marco catchment
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
dir2  = "$dir/temp_data/sanmarco"                   # Secondary directory (e.g., output storage)
## dir3  = "$dir/input_data"                           # Input data directory
fname = string("_testing")                          # Output files
using ReverseDiffSource, PyPlot #, ForwardDiff      # Automated differentiation package
pygui(true)
srand(18072011)                                     # Seeding the RNG provides reproducibility
# srand(time())                                     # Seeding the RNG with actual time provides randomness
##
## ============================================================================================
## Load real data
## ============================================================================================
##
y   = vec(readdlm("$dir/disc.dat"))  # Discharge data
t_r = readdlm("$dir/rain.dat")       # Time and rain data
const t = t_r[:,1]  # Time points
const r = t_r[:,2]  # Rain
##
## ============================================================================================
## System parameters
## ============================================================================================
##
const params = 2         # Number of system parameters (k, gamma)
const n = 730            # n+1 -> "end point" beads (see Tuckerman '93), n -> number of segments
const j = 10             # n(j-1) -> total number of staging beads, j-1 -> staging beads per segment
const N = int64(n*j+1)   # Total number of discrete time points = n*j + 1
                         # IMPORTANT: (N-1)/j = integer = n (measurement points)
if  N != size(t,1) 
    error("Be careful, N must be equal to size(t) (or size(r)) !")
end

const T  = t[end]-t[1]                     # Total time interval
const dt = T/(N-1)                         # Time step
const ty = iround(linspace(1, N, n+1))     # Indeces of "boundary" beads (= measurement points)     

const nsample_burnin = 0                   # Number of points in the MCMC
const nsample_eff    = 5000
const nsample        = nsample_eff + nsample_burnin

const dtau           = 0.020               # MD time step
const n_napa         = 5                   # Number of NAPA time steps

const true_K     = 0.0                     # Retention time
const true_gam   = 0.0                     # Dimensionless noise parameter
# const true_bet   = sqrt(T*true_gam/true_K)
# const sigma      = 0.10                    # Measurement noise
# const true_theta = [true_bet,true_gam]     # Parameters to be inferred (beta, tau)

# Initial state
K     = 100.0           
gam   = 0.5
bet   = sqrt(T*gam/K)                                   
theta = [bet, gam]
sigma = 0.10

# Logistic prior parameters
# const s_gam = 4
# const gam_0 = 3
# const s_k   = 0.04
# const k_0   = 500

# Parameter limits
K_min   = 0.0
gam_min = 0.0
K_max   = 2000.0
gam_max = 10.0

#=
plt.figure()
plt.xlabel("time")
plt.ylabel("r, y")
plt.title("System I/O")
plt.plot(t, r, "b", linewidth=2)
plt.plot(t[ty], y, "ro", markersize = 8, linewidth=5)
=#
##                 
## ============================================================================================
## Initialization of containers
## ============================================================================================
##
mp           = zeros(params+N)          # Masses
p            = zeros(params+N)          # Momenta

theta_sample = zeros(nsample+1, params) # Container for sampled parameters
u_sample     = zeros(nsample+1, N)      # Container for sampled coordinates
energies     = zeros(nsample)           # Container for sampled energies

theta_save   = zeros(params)            # To store temporary parameters
u_save       = zeros(N)                 # To store temporary coordinates

qs           = zeros(nsample+1, N)      # Sampled coordinates
predy        = zeros(nsample+1, N)      # Sampled outputs
##
## ============================================================================================
## Rain input, log derivative
## ============================================================================================
##    
lnr_der = Array(Float64,N)                      
for i = 1:(N-1)  lnr_der[i] = ( log(r[i+1]) - log(r[i]) ) / dt  end
##
##
## ============================================================================================
## Hamiltonian Monte Carlo
## ============================================================================================
##
## --------------------------------------------------------------------------------------------
## Initial state 
## -------------
## bq = ln(y/r) = beta*q
## N.B.: using the definitions S = (T g r / b^2)*exp(bq) and b = sqrt(Tg/K)
## ==> beta*q = ln(S/(kr)) = ln(y/r)
## --------------------------------------------------------------------------------------------
bq = Array(Float64, n+1)
q  = Array(Float64, N) 

bq = log(y./r[ty])
for s = 1:n
    q[ty[s]:ty[s+1]] = (1/bet)*linspace(bq[s],bq[s+1],ty[s+1]-ty[s]+1)
end

# An initial system realization
S = K*r[ty].*exp(bet*q[ty])
#=
plt.plot(t[ty], S, "r", linewidth = 2)
plt.plot(t, r, "b", linewidth=2)
=#

## Transformations q -> u 
## (eqs 2.16-17 in Tuckerman et al., JCP 99 (4), 2796, 1993)
u  = Array(Float64, N) 
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

m_bdy    = 700                      # m = m_q / dt
m_stg    = 300                      # we assume m_q prop. to dt ==> m = costant     
m_theta  = [150, 150]

mp[1:params] = m_theta
for s = 1:n  
    mp[params+(s-1)*j+1] = m_bdy                       
    for k = 2:j 
        mp[params+(s-1)*j+k] = m_stg                   
    end
end
mp[params+N] = m_bdy                                  

## --------------------------------------------------------------------------------------------
## Loading potentials, derivatives and Napa 
## --------------------------------------------------------------------------------------------

println(string("\nLoading potentials, derivatives and napa algorithm...\n---------------------------"))
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
    p = sqrt(mp).*randn(params+N)
    
    # Calculate energy:
    # -------------------
    H_old = V_N_fun(u) + V_n_fun(theta,u) + V_1_fun(theta,u) + sum((p.^2)./(2*mp))
    energies[counter] = H_old
    if  isnan(H_old) 
        error(string("Iteration ", counter, " --> energy values diverged..."))
    end
    
    # Save current state:
    # -------------------
    for i = 1:params theta_save[i] = theta[i] end
    for i = 1:N      u_save[i]     = u[i]     end

    # MD Integration:
    # -------------------
    for counter_napa = 1:n_napa 
        napa(theta, u, counter) 
    end 

    # Metropolis step and energy of proposal state
    # -------------------
    if (theta[2] <= gam_min) || (theta[2] >= gam_max) || ((T*theta[2]/(theta[1])^2) <= K_min) || ((T*theta[2]/(theta[1])^2) >= K_max) 
        for i = 1:params theta[i] = theta_save[i] end
        for i = 1:N      u[i]     = u_save[i]     end
        reject_counter_burnin += 1
    elseif any(u .> 10.0) || any(u .< -10.0)
        for i = 1:params theta[i] = theta_save[i] end
        for i = 1:N      u[i]     = u_save[i]     end
        reject_counter_burnin += 1
    else
        H_new = V_N_fun(u) + V_n_fun(theta,u) + V_1_fun(theta,u) + sum((p.^2)./(2*mp))
        accept_prob = min(1,exp(H_old-H_new))
        if rand() > accept_prob 
            for i = 1:params theta[i] = theta_save[i] end
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

m_bdy    = 700                      # m = m_q / dt
m_stg    = 300                      # we assume m_q prop. to dt ==> m = costant     
m_theta  = [150, 150]

mp[1:params] = m_theta
for s = 1:n  
    mp[params+(s-1)*j+1] = m_bdy                       
    for k = 2:j 
        mp[params+(s-1)*j+k] = m_stg                   
    end
end
mp[params+N] = m_bdy          

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
    p = sqrt(mp).*randn(params+N)

    # Calculate energy:
    # -------------------
    H_old = V_N_fun(u) + V_n_fun(theta,u) + V_1_fun(theta,u) + sum((p.^2)./(2*mp))
    energies[counter] = H_old
    if  isnan(H_old) 
        error(string("Iteration ", counter, " --> energy values diverged..."))
    end
    
    # Save current state:
    # -------------------
    for i = 1:params theta_save[i] = theta[i] end
    for i = 1:N      u_save[i]     = u[i]     end

    # MD Integration:
    # -------------------
    t1=time()
    for counter_napa = 1:n_napa 
        napa(theta, u, counter) 
    end 
    time_respa[counter-nsample_burnin] = time()-t1


    # Metropolis step and energy of proposal state
    # -------------------
    if (theta[2] <= gam_min) || (theta[2] >= gam_max) || ((T*theta[2]/(theta[1])^2) <= K_min) || ((T*theta[2]/(theta[1])^2) >= K_max) 
        for i = 1:params theta[i] = theta_save[i] end
        for i = 1:N      u[i]     = u_save[i]     end
        reject_counter += 1
    elseif any(u .> 10.0) || any(u .< -10.0)
        for i = 1:params theta[i] = theta_save[i] end
        for i = 1:N      u[i]     = u_save[i]     end
        reject_counter += 1
    else
        H_new = V_N_fun(u) + V_n_fun(theta,u) + V_1_fun(theta,u) + sum((p.^2)./(2*mp))
        accept_prob = min(1,exp(H_old-H_new))
        if rand() > accept_prob 
            for i = 1:params theta[i] = theta_save[i] end
            for i = 1:N      u[i]     = u_save[i]     end
            reject_counter += 1
        end
    end
    
    theta_sample[counter+1,:] = theta
    u_sample[counter+1,:]     = u 
   
    # if (counter%100 == 0)
        println(string(counter, " loops completed in ", round(time()-tinit,1), " seconds \n"))
    # end

end

# writedlm("C:/Users/ulzegasi/respatimeser.dat",      time_respa)
# writedlm("C:/Users/ulzegasi/respaslowser.dat",      time_respa_s)
# writedlm("C:/Users/ulzegasi/respafastforceser.dat", time_respa_f_d)
# writedlm("C:/Users/ulzegasi/respafastupser.dat",    time_respa_f_u)

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
param_names  = vcat("N", "j", "n", "t[1]", "dt", "params", "true_K", "true_gam", "sigma", "K", "gam", 
	"nsample_burnin", "nsample_eff", "m_bdy_burnin", "m_bdy", "m_theta_burnin", "m_theta_bet", 
	"m_theta_gam", "m_stg_burnin", "m_stg", "dtau", "n_napa")
param_values = vcat(N, j, n, t[1], dt, params, true_K, true_gam, sigma, K, gam, 
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

red_range = 200; red_predy = Array(Float64, red_range, N); row_ind = 1
if red_range < (nsample+1)
    for ip in iround(linspace(1, nsample+1, red_range))
        red_predy[row_ind,1:N] = predy[ip,1:N]
        row_ind += 1
    end
else
    red_predy = predy
end
writedlm("$dir2/predy$fname.dat", red_predy)

thetas      = Array(Float64,nsample+1,params)
thetas[:,1] = T*theta_sample[:,2]./(theta_sample[:,1]).^2  # K
thetas[:,2] = theta_sample[:,2]                            # gamma
writedlm("$dir2/thetas$fname.dat", thetas)

writedlm("$dir2/energies$fname.dat", energies)

reject_rates = [reject_counter_burnin/nsample_burnin*100,reject_counter/nsample_eff*100]
writedlm("$dir2/reject$fname.dat", reject_rates)

S_first  = Array(Any, N+1); S_first = vcat("System_realization", zeros(N)) # S_first  = vcat("System_realization", S)
rain_in  = Array(Any, N+1); rain_in  = vcat("Rain_input", r)
flow_out = Array(Any, n+2); flow_out = vcat("Output_flow", y)
io_data  = Array(Any, 2*N+n+4); io_data  = vcat(S_first, rain_in, flow_out)
writedlm("$dir2/iodata$fname.dat", io_data)
