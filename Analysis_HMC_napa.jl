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
dir      = "C:/Users/ulzegasi/Julia_files/ParInf_HMC/temp_data"
fname    = string("_testing")
##
##
## ============================================================================================
## Get stored parameters and data
## ============================================================================================
##
##
params   = readdlm("$dir/params$fname.dat")

N              = int64(params[findin(params[:,1],["N"]),2])[1]
j              = int64(params[findin(params[:,1],["j"]),2])[1]
n              = int64(params[findin(params[:,1],["n"]),2])[1] 
params 		   = int64(params[findin(params[:,1],["params"]),2])[1] 
nsample_burnin = int64(params[findin(params[:,1],["nsample_burnin"]),2])[1]
nsample_eff    = int64(params[findin(params[:,1],["nsample_eff"]),2])[1]  
n_napa         = int64(params[findin(params[:,1],["n_napa"]),2])[1] 

nsample = nsample_burnin + nsample_eff

dt 			   = params[findin(params[:,1],["dt"]),2][1]
true_K		   = params[findin(params[:,1],["true_K"]),2][1]
true_gam       = params[findin(params[:,1],["true_gam"]),2][1]
sigma          = params[findin(params[:,1],["sigma"]),2][1]
K              = params[findin(params[:,1],["K"]),2][1]
gam            = params[findin(params[:,1],["gam"]),2][1]
dtau           = params[findin(params[:,1],["dtau"]),2][1]

m_bdy_burnin   = params[findin(params[:,1],["M_bdy_burnin"]),2][1]
m_bdy          = params[findin(params[:,1],["M_bdy"]),2][1]
m_theta_burnin = params[findin(params[:,1],["M_theta_burnin"]),2][1]
m_theta        = [params[findin(params[:,1],["m_theta_bet"]),2][1], params[findin(params[:,1],["M_theta_gam"]),2][1]] 
m_stg_burnin   = params[findin(params[:,1],["m_stg_burnin"]),2][1]
m_stg          = params[findin(params[:,1],["m_stg"]),2][1]

t    = Array(Float64, N)
t[1] = params[findin(params[:,1],["t[1]"]),2][1]
T = dt*(N-1)
for i = 2:N
	t[i] = t[i-1] + dt 
end
ty  = iround(linspace(1, N, n+1))   # Indexes of measurement points

qs       = readdlm("$dir/last_qs$fname.dat")
predy    = readdlm("$dir/predy$fname.dat")
thetas   = readdlm("$dir/thetas$fname.dat")
energies = readdlm("$dir/energies$fname.dat")
reject   = readdlm("$dir/reject$fname.dat")

io_data  = readdlm("$dir/iodata$fname.dat")
S = io_data[2:N+1]
r = io_data[(N+3):(2*N+2)]
y = io_data[(2*N+4):end]

using PyPlot, Winston, ForwardDiff, Distributions, KernelDensity
pygui(true)

##
##
## ============================================================================================
## Get ready to plot!
## ============================================================================================
##
##

## --------------------------------------------------------------------------------------------
## System plots:
## --------------------------------------------------------------------------------------------

plt.figure(1)
plt.subplots_adjust(hspace=0.5)
plt.subplot(211)
plt.ylabel("S[t]")
plt.title("Water volume")
plt.plot(t, S,"r")
plt.subplot(212)
plt.xlabel("time")
plt.ylabel("r[t]")
plt.title("Rain input")
plt.plot(t,r,"b")
plt.savefig("$dir/figure1$fname.png",transparent=true)

plt.figure(2)
#axes()[:set_ylim]([0,1.5])
plt.xlabel("time")
plt.ylabel("r, S/K, y")
plt.title("System I/O")
plt.plot(t,r,"b--",label="Rain input",linewidth=2)
plt.plot(t,S/true_K,"g",label="Output",linewidth=2)
yerr=2*sigma*ones(length(y))
plt.errorbar(t[ty], y, yerr=(yerr,yerr), fmt="o", color="b", 
                                            capsize=4, elinewidth=2,label="Data points")
# plt.plot(t[ty],y,"bo",label="Data points")
# plt.legend(loc="upper right",fancybox="true")
plt.savefig("$dir/figure2$fname.png",transparent=true)

## --------------------------------------------------------------------------------------------
## Spaghetti plot:
## --------------------------------------------------------------------------------------------

range    = 1:N
redrange = iround(linspace(1,n+1,min(n+1,101)))         # Number of data points that will be shown
                                                        # in the plot
maxy=Array(Float64,N)
miny=Array(Float64,N)
maxy=vec(maximum(predy,1))                              # For each time point, max and min
miny=vec(minimum(predy,1))                              # of all the sampled system realizations

plt.figure(3)
plt.subplots_adjust(hspace=0.5)
plt.subplot(211)
plt.xlabel("N")
plt.ylabel("Output")
plt.title("Sampled system dynamics vs. Data points")
plt.plot(1:N, maxy, label="Sampled system realizations", color="g", linewidth=1)
plt.plot(1:N, miny, color="g", linewidth=1)
plt.fill_between(1:N,maxy,miny,color="g",alpha=0.5)
yerr=2*sigma*ones(length(redrange))
plt.errorbar(ty[redrange], y[redrange], yerr=(yerr,yerr), fmt="o", color="r", 
                                            capsize=4, elinewidth=2,label="Data points")
# plt.legend(loc="lower right",fancybox="true")

plt.subplot(212)
plt.xlabel("N")
plt.ylabel("Output")
plt.title("Spaghetti style")
for ip = 1:size(predy,1)
	plt.plot(1:N, vec(predy[ip,1:N]), label="Sampled system realizations", color="g", linewidth=1)
end
yerr=2*sigma*ones(length(redrange))
plt.errorbar(ty[redrange], y[redrange], yerr=(yerr,yerr), fmt="o", color="r", 
                                            capsize=4, elinewidth=2,label="Data points")
plt.savefig("$dir/figure3$fname.png",transparent=true,dpi=300)

## --------------------------------------------------------------------------------------------
## Plot final state (and QQ plot):
## --------------------------------------------------------------------------------------------

range = (nsample_burnin+2):(nsample+1)
plt.figure(4)
plt.subplots_adjust(hspace=0.5)
plt.subplot(211)
plt.xlabel("S / K = r Exp(q)")
plt.ylabel("Probability density")
plt.title("Final state S / K")
plt.hist(r[N]*exp(qs[range]), 50, normed = 1, color = "y")
kd = KernelDensity.kde(r[N]*exp(qs[range]))
plt.plot(kd.x, kd.density, color = "r", linewidth = 3, label = "Kernel density kde(exp(u))")
mu  = mean(r[N]*exp(qs[range]))
sig = std(r[N]*exp(qs[range]))
xx  = linspace(mu-4sig, mu+4sig, 100)
plt.plot(xx, pdf(Normal(mu,sig),xx), "b--", linewidth = 3, label = "Normal distribution")
# plt.legend(loc="upper right",fancybox="true")

plt.subplot(212)
plt.xlabel("exp(u)")
plt.ylabel("Normal(mu,sig)")
plt.title("Final state QQ plot")
res = qqbuild(r[N]*exp(qs[range]),Normal(mu,sig))
plt.plot(res.qy, res.qx, linewidth = 3)
xy = [x for x in linspace(minimum([res.qx,res.qy]),maximum([res.qx,res.qy]),4)]
plt.plot(xy, xy, "r--", linewidth = 3)
plt.savefig("$dir/figure4$fname.png",transparent=true)

## --------------------------------------------------------------------------------------------
## Plot parameters:
## --------------------------------------------------------------------------------------------


plt.figure(5)
plt.subplots_adjust(hspace=0.5)
plt.subplot(211)
plt.ylabel("K")
plt.title("K chain")
plt.plot(thetas[:,1],"r",linewidth=3)
plt.plot(thetas[:,1],"ro",markersize=4)
plt.subplot(212)
plt.xlabel("N")
plt.ylabel("Gamma")
plt.title("Gamma chain")
plt.plot(thetas[:,2],"b",linewidth=3)
plt.plot(thetas[:,2],"bo",markersize=4)
plt.savefig("$dir/figure5$fname.png",transparent=true)


plt.figure(6)
plt.xlabel("K")
plt.ylabel("Gamma")
plt.title("Evolution in K-Gamma space")
grid = iround(linspace(1, size(thetas)[1], 1000))       # Number of points in K-gamma space to be plotted
grid2 = grid[length(grid)-49:length(grid)]              # Last points in K-gamma space to be highlighted
plt.plot(thetas[grid,1], thetas[grid,2], linewidth = 0.5)
plt.plot(thetas[grid,1], thetas[grid,2], "b.")
plt.plot(thetas[1,1], thetas[1,2], "go", markersize = 10, label = "Initial state")
plt.plot(thetas[grid2,1], thetas[grid2,2], "yo", markersize = 10, label = "Last points")
plt.plot(thetas[end,1], thetas[end,2], "ro", markersize = 10, label = "Last state")
plt.plot(true_K, true_gamma, "rs", markersize = 20, label = "True value", linewidth = 2)
# plt.legend(loc="lower right",fancybox="true")
plt.savefig("$dir/figure6$fname.png",transparent=true,dpi=300)


lowlim = nsample_burnin                                 # To possibly discard initial 
                                                        # system realizations (burn-in)
plt.figure(7)
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

plt.savefig("$dir/figure7$fname.png",transparent=true)

## --------------------------------------------------------------------------------------------
## Plot energies:
## --------------------------------------------------------------------------------------------

plt.figure(8)
plt.xlabel("N")
plt.ylabel("Energy")
plt.title("Energies")
plt.plot(energies)
plt.savefig("$dir/figure8$fname.png",transparent=true)

reject
