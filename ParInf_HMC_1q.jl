## ============================================================================================
## Project: Parameter Inference for a Simple Stochastic Hydrological Model
## File: ParInf_HMC_1q.jl
##
## Description: Single bucket, scale invariant noise, no evaporation yet
##              dS(t)/dt = r(t) - S(t)/k + sqrt(gamma/k) S(t) eta(t)
## 
## Only the initial coordinate set {q} is used. 
## Transformations {q} <-> {u} are not used here.
## (Tuckerman et al., JCP 99 (4), 2796, 1993)
##
##				
## GitHub repository: https://github.com/ulzegasi/Julia_1.git
##
## Authors: Simone Ulzega, Carlo Albert
##
## simone.ulzega@eawag.ch
## carlo.albert@eawag.ch
## ============================================================================================
##
##
##
##
## ============================================================================================
## System Variables:
## ============================================================================================
##
##

dir = "C:/Users/ulzegasi/Julia_files/Julia_1"           # Main project directory
dir2 = "C:/Users/ulzegasi/Julia_files/Julia_1"          # Files could be saved in a different directory

range = 2:5002

# Time points t.dat
t  = float64(readdlm("$dir/t.dat")[range,2])

N  = length(t)                         # Total number of discrete time points
T  = t[end]-t[1]                       # Time interval
dt = T/(N-1)                           # Time step

const j = 10                           # j-1 = number of staging beads per segment (see Tuckerman '93)

if  ((N-1)%j) != 0 
    error("Be careful, the number of staging points j must fulfill (N-1)/j = integer !")
end

const n = int64((N-1)/j)               # n + 1 = number of "end point" beads (see Tuckerman '93)
                                       # n = number of segments
                                       # n (j-1) = total number of staging beads

ty  = iround(linspace(1, N, n+1))      # Indeces of "end point" beads (= "measurement" points)

const s = 2                            # Number of system parameters (k, gamma)     

##
## ============================================================================================
## Parameters to be inferred:
## ============================================================================================
##
##

true_K     = 200.                                       # Retention time
true_gamma = 0.2                                        # Dimensionless noise parameter
sigma      = 0.01                                       # Measurement noise

true_theta = [log(true_K/T),log(true_gamma)]            # Parameters to be inferred (beta, tau)

##
##
## ============================================================================================
## File names:
## ============================================================================================
##
##

fname= string("_K$true_K","_G$true_gamma","_S$sigma","_Rsin")          # This will be displayed 
                                                                       # in the saved file names
# fname= string("")                                                    # Alternative: empty string                                       

##
##
## ============================================================================================
## Loading libraries:
## ============================================================================================
##
##

using PyPlot, Winston, ForwardDiff, Distributions, KernelDensity
pygui(true)
reload("$dir/1D_functions2_2_SU.jl")

##
##
## ============================================================================================
## Rain input:
## ============================================================================================
##
##

r   = Array(Float64,N);
r   = sin(t/100.).*sin(t/100.) + 0.1;                   # Simple sinusoidal model

# r=float64(readdlm("$dir/r.dat")[range,2])

lnr_der = zeros(N)                                      # Calculate logarithmic derivative of rain input
for i=1:(N-1)                                           # for later use
    lnr_der[i] = ( log(r[i+1]) - log(r[i]) ) / dt     
end

##
##
## ============================================================================================
## System realization (Ito convention):
## ============================================================================================
##
##

srand(18072011)                                         # Seeding the RNG provides reproducibility
#srand(time())                                          # Seeding the RNG with actual time provides randomness
S = Array(Float64,N)
S[1] = true_K * r[1]   # --------------------------->   # Unperturbed steady state (with constant rain input)
for i = 2 : N
    S[i] = S[i-1] + dt * ( r[i-1] - S[i-1]/true_K ) +
                   sqrt(dt) * sqrt(true_gamma/true_K) * S[i-1] * randn()
end


Sr = zeros(N)                                           # System realization WITHOUT noise, 
constr = 0.5                                            # WITH CONSTANT rain. Just for fun.
Sr[1] = 0
for i = 2 : N
    Sr[i] = Sr[i-1] + dt * ( constr - Sr[i-1]/true_K )
end

y  = max(0.01,(1/true_K)*S[ty] + sigma*randn(n+1))      # Generation of measurement points


plt.figure(1)
plt.subplots_adjust(hspace=0.5)
plt.subplot(211)
plt.ylabel("S[t]")
plt.title("System realization")
plt.plot(t, S,"r")
plt.subplot(212)
plt.xlabel("time")
plt.ylabel("r[t]")
plt.title("Rain input")
plt.plot(t,r,"b")
plt.savefig("$dir2/figure1$fname.png",transparent=true)

plt.figure(2)
#axes()[:set_ylim]([0,1.5])
plt.xlabel("time")
plt.ylabel("r[t], S[t]/K")
plt.title("System I/O")
plt.plot(t,r,"b",label="Rain input")
plt.plot(t,S/true_K,"g",label="Output")
plt.legend(loc="upper right",fancybox="true")
plt.savefig("$dir2/figure2$fname.png",transparent=true)

plt.figure(3)
axes()[:set_ylim]([0,120])
plt.xlabel("time")
plt.ylabel("S[t]")
plt.title("Unperturbed system realization - constant rain")
plt.plot(t,Sr,"r")
plt.hlines(true_K*constr, 0, t[end], colors = "g", linestyle = "dotted", label = "Steady state")
plt.text(10,102,"Steady state")
plt.legend(loc="lower right",fancybox="true")
plt.savefig("$dir2/figure3$fname.png",transparent=true)

plt.figure(4)
#axes()[:set_ylim]([0,120])
plt.xlabel("time")
plt.ylabel("S[t]")
plt.title("Measurements")
plt.plot(t[ty],y,"bo",label="data points")
plt.plot(t,S/true_K,"r",label="system realization")
plt.legend(loc="lower right",fancybox="true")
plt.savefig("$dir2/figure4$fname.png",transparent=true)

##
##
## ============================================================================================
## Transformations (eqs 2.16-17 in Tuckerman et al., JCP 99 (4), 2796, 1993):
## ============================================================================================
##
##
##
##q = log(S./(true_K*r))
##
##u = zeros(N)
##for ss=0:(n-1)
##    u[ss*j+1] = q[ss*j+1]
##    for k=2:j
##        u[ss*j+k] = q[ss*j+k] - ( (k-1)*q[ss*j+k+1] + q[ss*j+1] )/k
##    end
##end
##u[N] = q[N]                                           # The trasformations in Tuckerman 93
                                                        # do NOT have the last "end point" bead


##
##
## ============================================================================================
## Hamiltonian Monte Carlo:
## ============================================================================================
##
##

## Containers for HMC:
## -------------------

nsample         = 200
theta_sample    = Array(Float64,nsample,s)
q_sample        = Array(Float64,nsample,N)
energies        = Array(Float64,nsample)

## Masses:
## -------------------

M_theta = 100.
M_bdy   = 1.*true_K*(N-1) / (true_gamma*T)
M_stage = M_bdy*0.01

mp = Array(Float64, s+N)
mp[1:s] = M_theta
for i=0:(n-1)  
    mp[s+1+i*j] = M_bdy                                 # End point beads
    for k=2:j 
        mp[s+i*j+k] = M_stage*k/(k-1)                   # Staging beads
    end
end
mp[s+N] = M_bdy                                         # Last end point bead

## Parameters for MD:
## -------------------

dtau    = .02 #.04                                      # Long range MD time step
nn      = 15		                                    # To define short-range time steps in RESPA
n_respa = 10                                            # Number of large time steps

## Initial state:
## -------------------

K     = 100.0
gamma = 0.3
sigma = 0.01
theta = [log(K/T),log(gamma)]

S_init = Array(Float64,N)
for a=1:n
    S_init[ty[a]:ty[a+1]] = K*linspace(y[a],y[a+1],ty[a+1]-ty[a]+1)
end

plt.figure(5)
plt.xlabel("time")
plt.ylabel("output")
plt.title("Initial state")
plt.plot(t[ty],y,"bo",label="data points")
plt.plot(t,S_init/K,"r",label="initial state",linewidth=2.0)
plt.plot(t,S/true_K,"g",label="system realization")
plt.legend(loc="lower right",fancybox="true")
plt.savefig("$dir2/figure5$fname.png",transparent=true)

## Transformations q -> u 
## (eqs 2.16-17 in Tuckerman et al., JCP 99 (4), 2796, 1993):
## -------------------

q = Array(Float64,N)
q = log(S_init./(K*r))

#u = Array(Float64,N)
#for i=0:(n-1)
#    u[i*j+1] = q[i*j+1]
#    for k=2:j
#        u[i*j+k] = q[i*j+k] - ( (k-1)*q[i*j+k+1] + q[i*j+1] )/k
#    end
#end
#u[N] = q[N]

## HMC loop:
## ------------------------------------------------------
## ------------------------------------------------------

reject_counter = 0

p          = Array(Float64,s+N)
theta_save = Array(Float64,s)
q_save     = Array(Float64,N)

println("\nStarting HMC loops... \n")
println("------------------------\n")
t1=time()

for counter = 1:nsample

    # Sample momenta:
    # -------------------

    p = sqrt(mp).*randn(s+N)
    
    # Calculate energy:
    # -------------------

    H_old = V_fast(theta,q) + V_slow(theta,q) + sum((p .* p) ./ (2*mp))
    energies[counter] = H_old
    
    # Save current state:
    # -------------------

    for i = 1:s   theta_save[i] = theta[i]    end
    for i = 1:N   q_save[i]     = q[i]        end
    
    # MD Integration:
    # -------------------

    for counter_respa = 1:n_respa 
        RESPA(theta,q,p,mp,dtau,nn) 
    end 
    
    # Calculate energy of proposal state:
    # -------------------

    H_new = V_fast(theta,q) + V_slow(theta,q) + sum(p .* p ./ (2*mp))

    # Metropolis step:
    # -------------------

    accept_prob = min(1,exp(H_old-H_new))
    if rand() > accept_prob
        for i=1:s theta[i] = theta_save[i] end
        for i=1:N q[i]     = q_save[i] end
        reject_counter += 1
    end
   
    theta_sample[counter,:] = theta
    q_sample[counter,:] = q 
    
    if (counter%10 == 0)
        println(string(counter, " loops completed in ", round(time()-t1,1), " seconds \n"))
    end

end

#t2=time()
#t_hmc = t2-t1

## End of HMC loop
## ------------------------------------------------------
## ------------------------------------------------------


## Back transformations (u -> q):
## (Tuckerman et al., JCP 99 (1993), eq. 2.19)
## ---------------------

###qs = Array(Float64, nsample-1, N)
###for index = 1:(nsample-1)
###    ss=n-1
###    while ss>=0
###        k=j+1
###        qs[index,ss*j+k] = u_sample[index,ss*j+k]
###        while k>2
###           k -= 1
###           qs[index,ss*j+k] = u_sample[index,ss*j+k] +
###                              u_sample[index,ss*j+1]/k +
###                              (k-1)*qs[index,ss*j+k+1]/k
###         end
###        ss -= 1
###    end
###    qs[index,1] = u_sample[index,1]
###end

ypred = Array(Float64,nsample-1,N)                              # Expected outputs
for index = 1:(nsample-1)
    for l=1:N
        ypred[index,l] = r[l]*exp(q_sample[index,l])
    end
end

println(string("\nRun completed in ", time()-t1, " seconds"))

##
##
## ============================================================================================
## Final plots
## ============================================================================================
##
##

## Spaghetti plot:
## ---------------------

range    = 1:N
redrange = 1:501

maxy=Array(Float64,N)
miny=Array(Float64,N)
maxy=vec(maximum(predy,1))                              # For each time point, max and min
miny=vec(minimum(predy,1))                              # of all the sampled system realizations

plt.figure(6)
#for a in iround( linspace(1, nsample, 100) )
#    plt.plot(range,vec(predy[a,range]),"g",alpha=0.3)
#end
plt.xlabel("N")
plt.ylabel("Output")
plt.title("Sampled system dynamics vs. Data points")
plt.plot(range, maxy, label="Sampled system realizations", color="g", linewidth=1)
plt.plot(range, miny, color="g", linewidth=1)
plt.fill_between(range,maxy,miny,color="g",alpha=0.5)
yerr=2*sigma*ones(length(redrange))
plt.errorbar(ty[redrange], y[redrange], yerr=(yerr,yerr), fmt="o", color="r", 
                                            capsize=4, elinewidth=2,label="Data points")
plt.legend(loc="lower right",fancybox="true")
plt.savefig("$dir2/figure6$fname.png",transparent=true,dpi=300)


## Plot final state (and QQ plot):
## ---------------------


plt.figure(7)
plt.subplots_adjust(hspace=0.5)
plt.subplot(211)
plt.xlabel("exp(q)")
plt.ylabel("Probability density")
plt.title("Final state")
plt.hist(exp(q_sample[:,N]), 50, normed = 1, color = "y")
kd = KernelDensity.kde(exp(q_sample[:,N]))
plt.plot(kd.x, kd.density, color = "r", linewidth = 3, label = "Kernel density kde(exp(u))")
mu = mean(exp(q_sample[:,N]))
sig = std(exp(q_sample[:,N]))
xx = linspace(mu-4sig, mu+4sig, 100)
plt.plot(xx, pdf(Normal(mu,sig),xx), "b--", linewidth = 3, label = "Normal distribution")
plt.legend(loc="upper right",fancybox="true")



plt.subplot(212)

plt.xlabel("exp(q)")
plt.ylabel("Normal(mu,sig)")
plt.title("Final state QQ plot")
res = qqbuild(exp(q_sample[:,N]),Normal(mu,sig))
plt.plot(res.qy, res.qx, linewidth = 3)
xy = [x for x in linspace(minimum([res.qx,res.qy]),maximum([res.qx,res.qy]),4)]
plt.plot(xy, xy, "r--", linewidth = 3)
plt.savefig("$dir2/figure7$fname.png",transparent=true)



## Plot parameters:
## ---------------------


thetas = Array(Float64,nsample,s)
thetas[:,1] = exp(theta_sample[:,1])*T                  # K
thetas[:,2] = exp(theta_sample[:,2])                    # gamma

plt.figure(8)
plt.subplots_adjust(hspace=0.5)
plt.subplot(211)
plt.ylabel("K")
plt.title("K chain")
plt.plot(thetas[:,1],"ro")
plt.subplot(212)
plt.xlabel("N")
plt.ylabel("Gamma")
plt.title("Gamma chain")
plt.plot(thetas[:,2],"bo")
plt.savefig("$dir2/figure8$fname.png",transparent=true)


plt.figure(9)
plt.xlabel("K")
plt.ylabel("Gamma")
plt.title("Evolution in K-Gamma space")
grid  = iround(linspace(1, size(theta_sample)[1], 4000))    
grid2 = grid[length(grid)-49:length(grid)]
plt.plot(thetas[grid,1], thetas[grid,2], linewidth = 0.5)
plt.plot(thetas[grid,1], thetas[grid,2], "b.")
plt.plot(thetas[1,1], thetas[1,2], "go", markersize = 10, label = "Initial state")
plt.plot(thetas[grid2,1], thetas[grid2,2], "yo", markersize = 10, label = "End points")
plt.plot(thetas[size(thetas)[1],1], 
                   thetas[size(thetas)[1],2], "ro", markersize = 10, label = "Last state")
plt.legend(loc="lower right",fancybox="true")
plt.savefig("$dir2/figure9$fname.png",transparent=true,dpi=300)


lowlim = 1000

plt.figure(10)
plt.subplots_adjust(hspace=0.5)
plt.subplot(211)
plt.xlabel("K")
plt.ylabel("Probability density")
plt.title("PDF")
plt.hist(thetas[lowlim:end,1], 50, normed = 1, color = "y")
kd_K = KernelDensity.kde(thetas[lowlim:end,1])
plt.plot(kd_K.x, kd_K.density, color = "r", linewidth = 3, label = "Kernel density kde(K)")
mu_K   = mean(thetas[lowlim:end,1])
sig_K  = std(thetas[lowlim:end,1])
xx_K   = linspace(mu_K-4*sig_K, mu_K+4*sig_K, 100)
plt.plot(xx_K, pdf(Normal(mu_K,sig_K),xx_K), "b--", linewidth = 3, label = "Normal distribution")
plt.axvline(x=mu_K, linewidth=2, color="r")
plt.axvspan(mu_K-sig_K, mu_K+sig_K, facecolor="y", alpha=0.3)

plt.subplot(212)
plt.xlabel("Gamma")
plt.ylabel("Probability density")
plt.title("PDF")
plt.hist(thetas[lowlim:end,2], 50, normed = 1, color = "y")
kd_G = KernelDensity.kde(thetas[lowlim:end,2])
plt.plot(kd_G.x, kd_G.density, color = "r", linewidth = 3, label = "Kernel density kde(Gamma)")
mu_G   = mean(thetas[lowlim:end,2])
sig_G  = std(thetas[lowlim:end,2])
xx_G   = linspace(mu_G-4*sig_G, mu_G+4*sig_G, 100)
plt.plot(xx_G, pdf(Normal(mu_G,sig_G),xx_G), "b--", linewidth = 3, label = "Normal distribution")
plt.axvline(x=mu_G, linewidth=2, color="r")
plt.axvspan(mu_G-sig_G, mu_G+sig_G, facecolor="y", alpha=0.3)

plt.savefig("$dir2/figure10$fname.png",transparent=true)

# Plot energies:
# ---------------------

plt.figure(11)
plt.xlabel("N")
plt.ylabel("Energy")
plt.title("Energies")
plt.plot(energies)
plt.savefig("$dir2/figure11$fname.png",transparent=true)
