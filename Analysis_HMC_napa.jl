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
dir   = "C:/Users/ulzegasi/Julia_files/ParInf_HMC/temp_data"
dir2  = "C:/Users/ulzegasi/Julia_files/ParInf_HMC/input_data"
fname = string("_sinr")
##
##
## ============================================================================================
## Get stored parameters and data
## ============================================================================================
##
##
paramslist     = readdlm("$dir/params$fname.dat")

N              = int64(paramslist[findin(paramslist[:,1],["N"]),2])[1]
j              = int64(paramslist[findin(paramslist[:,1],["j"]),2])[1]
n              = int64(paramslist[findin(paramslist[:,1],["n"]),2])[1] 
nparams		   = int64(paramslist[findin(paramslist[:,1],["params"]),2])[1] 
nsample_burnin = int64(paramslist[findin(paramslist[:,1],["nsample_burnin"]),2])[1]
nsample_eff    = int64(paramslist[findin(paramslist[:,1],["nsample_eff"]),2])[1]  
n_napa         = int64(paramslist[findin(paramslist[:,1],["n_napa"]),2])[1] 

nsample = nsample_burnin + nsample_eff

dt 			   = paramslist[findin(paramslist[:,1],["dt"]),2][1]
true_K		   = paramslist[findin(paramslist[:,1],["true_K"]),2][1]
true_gam       = paramslist[findin(paramslist[:,1],["true_gam"]),2][1]
sigma          = paramslist[findin(paramslist[:,1],["sigma"]),2][1]
K              = paramslist[findin(paramslist[:,1],["K"]),2][1]
gam            = paramslist[findin(paramslist[:,1],["gam"]),2][1]
dtau           = paramslist[findin(paramslist[:,1],["dtau"]),2][1]

m_bdy_burnin   = paramslist[findin(paramslist[:,1],["m_bdy_burnin"]),2][1]
m_bdy          = paramslist[findin(paramslist[:,1],["m_bdy"]),2][1]
m_theta_burnin = paramslist[findin(paramslist[:,1],["m_theta_burnin"]),2][1]
m_theta        = [paramslist[findin(paramslist[:,1],["m_theta_bet"]),2][1], paramslist[findin(paramslist[:,1],["m_theta_gam"]),2][1]] 
m_stg_burnin   = paramslist[findin(paramslist[:,1],["m_stg_burnin"]),2][1]
m_stg          = paramslist[findin(paramslist[:,1],["m_stg"]),2][1]

t    = Array(Float64, N)
t[1] = paramslist[findin(paramslist[:,1],["t[1]"]),2][1]
T = dt*(N-1)
for i = 2:N
	t[i] = t[i-1] + dt 
end
ty  = iround(linspace(1, N, n+1))   # Indexes of measurement points

qs       = readdlm("$dir/last_qs$fname.dat")
predy    = readdlm("$dir/predy$fname.dat")
predq    = readdlm("$dir/predq$fname.dat")

thetas   = readdlm("$dir/thetas$fname.dat")
bet      = sqrt(T*thetas[:,2]./thetas[:,1])

energies = readdlm("$dir/energies$fname.dat")
reject   = readdlm("$dir/reject$fname.dat")

io_data  = readdlm("$dir/iodata$fname.dat")
S = io_data[2:N+1]
r = io_data[(N+3):(2*N+2)]
y = io_data[(2*N+4):end]

#=St = readdlm("$dir2/St$fname.dat")
S2  = St[:,1]
tS = St[:,2]=#

u_chains = readdlm("$dir/u_chains$fname.dat")
u1       = int64(u_chains[1,1])
u2       = int64(u_chains[1,2])
u3       = int64(u_chains[1,3])
u4       = int64(u_chains[1,4])
u_chains = u_chains[2:end,:]

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

plt.figure(1,figsize=(12.5, 8.5))
plt.subplots_adjust(hspace=0.5)
plt.subplot(211)[:set_xlim]([-5,605])
plt.grid("on")
plt.ylabel("S[t]")
plt.title("Water volume")
plt.plot(t, S,"r")
plt.subplot(212)[:set_xlim]([-5,605])
plt.grid("on")
plt.xlabel("time")
plt.ylabel("r[t]")
plt.title("Rain input")
plt.plot(t,r,"b")
plt.savefig("$dir/figure1$fname.png",transparent=true)

plt.figure(2,figsize=(12.5, 8.5))
plt.subplots_adjust(hspace=0.5)
plt.subplot(211)[:set_xlim]([-5,605])
plt.grid("on")
plt.ylabel("y[t]")
plt.title("Runoff")
plt.plot(t, S/true_K, color = (1,0.65,0), linewidth=2)
yerr=2*sigma*y
plt.errorbar(t[ty], y, yerr=(yerr,yerr), fmt="o", markersize = 12, color="r", capsize=10, elinewidth=4)
plt.tick_params(length=5, width=3)
plt.subplot(212)[:set_xlim]([-5,605])
plt.grid("on")
plt.xlabel("time")
plt.ylabel("r[t]")
plt.title("Rain input")
plt.plot(t,r,"b")
plt.tick_params(length=5, width=3)
plt.savefig("$dir/figure2$fname.png",transparent=true)

## --------------------------------------------------------------------------------------------
## Spaghetti plot:
## --------------------------------------------------------------------------------------------

range    = 1:N
redrange = iround(linspace(1, n+1, min(n+1,101))) # Number of data points that will be shown
                                                  # in the plot
#=
maxy=Array(Float64,N)
miny=Array(Float64,N)
maxy=vec(maximum(predy,1))                        # For each time point, max and min
miny=vec(minimum(predy,1))                        # of all the sampled system realizations
=#

#=
plt.figure(3)
plt.subplots_adjust(hspace=0.5)
plt.subplot(211)
plt.xlabel("N")
plt.ylabel("Output")
plt.title("Sampled system dynamics vs. Data points")
plt.plot(range, maxy, label="Sampled system realizations", color="g", linewidth=1)
plt.plot(range, miny, color="g", linewidth=1)
plt.fill_between(range,maxy,miny,color="g",alpha=0.5)
yerr=2*sigma*y[redrange]
plt.errorbar(ty[redrange], y[redrange], yerr=(yerr,yerr), fmt="o", color="r", 
                                            capsize=4, elinewidth=2,label="Data points")
# plt.legend(loc="lower right",fancybox="true")

plt.subplot(212)
plt.xlabel("N")
plt.ylabel("Output")
plt.title("Spaghetti style")
for ip = 1:size(predy,1)
	plt.plot(range, vec(predy[ip,range]), label="Sampled system realizations", color="g", linewidth=1)
end
yerr=2*sigma*y[redrange]
plt.errorbar(ty[redrange], y[redrange], yerr=(yerr,yerr), fmt="o", color="r", 
                                            capsize=4, elinewidth=2,label="Data points")
plt.savefig("$dir/figure3$fname.png",transparent=true,dpi=300)
=#

plt.figure(3,figsize=(12.5, 8.5))
plt.subplots_adjust(hspace=0.5)
plt.subplot(211)[:set_xlim]([-5,605])
plt.xticks(size="15")
plt.yticks(size="15") 
plt.tick_params(length=5, width=2)
plt.xlabel("N")
plt.ylabel("Output")
plt.title("Spaghetti plot - q")
for ip = 1:size(predq,1)
	plt.plot(range, vec(predq[ip,range]), label="Sampled system realizations", color="g", linewidth=1)
end
qerr=2*sigma*ones(length(redrange))
plt.errorbar(ty[redrange], vec(predq[1,ty[redrange]]), yerr=(qerr,qerr), fmt="o", color="r", markersize = 12, capsize=8, elinewidth=4)
plt.subplot(212)[:set_xlim]([-5,605])
plt.xticks(size="15")
plt.yticks(size="15") 
plt.tick_params(length=5, width=2)
plt.xlabel("N")
plt.ylabel("Output")
plt.title("Spaghetti plot - y")
for ip = 1:size(predy,1)
	plt.plot(range, vec(predy[ip,range]), label="Sampled system realizations", color="g", linewidth=1)
end
yerr=2*sigma*y[redrange]
plt.errorbar(ty[redrange], y[redrange], yerr=(yerr,yerr), fmt="o", color="r", markersize = 12, capsize=8, elinewidth=4)
plt.savefig("$dir/figure3$fname.png",transparent=true,dpi=300)

## --------------------------------------------------------------------------------------------
## Plot final state (and QQ plot):
## --------------------------------------------------------------------------------------------

range = (nsample_burnin+2):(nsample+1)
plt.figure(4)
plt.subplots_adjust(hspace=0.5)
plt.subplot(211)
plt.xlabel("S / K = r Exp(bq)")
plt.ylabel("Probability density")
plt.title("Final state S / K")
plt.hist(r[N]*exp(bet[range].*qs[range]), 50, normed = 1, color = "y")
kd = KernelDensity.kde(r[N]*exp(bet[range].*qs[range]))
plt.plot(kd.x, kd.density, color = "r", linewidth = 3, label = "Kernel density kde(exp(u))")
mu  = mean(r[N]*exp(bet[range].*qs[range]))
sig = std(r[N]*exp(bet[range].*qs[range]))
xx  = linspace(mu-4sig, mu+4sig, 100)
plt.plot(xx, pdf(Normal(mu,sig),xx), "b--", linewidth = 3, label = "Normal distribution")
# plt.legend(loc="upper right",fancybox="true")

plt.subplot(212)
# plt.xlabel("exp(u)")
plt.ylabel("Normal(mu,sig)")
plt.title("Final state QQ plot")
res = qqbuild(r[N]*exp(bet[range].*qs[range]),Normal(mu,sig))
plt.plot(res.qy, res.qx, linewidth = 3)
xy = [x for x in linspace(minimum([res.qx,res.qy]),maximum([res.qx,res.qy]),4)]
plt.plot(xy, xy, "r--", linewidth = 3)
plt.savefig("$dir/figure4$fname.png",transparent=true)


#=
range = (nsample_burnin+2):(nsample+1)
plt.figure(4)
plt.subplots_adjust(hspace=0.5)
plt.subplot(211)
plt.xlabel("q")
plt.ylabel("Probability density")
plt.title("Final state q")
plt.hist(qs[range], 50, normed = 1, color = "y")
kd = KernelDensity.kde(qs[range])
plt.plot(kd.x, kd.density, color = "r", linewidth = 3, label = "Kernel density kde(exp(u))")
mu  = mean(qs[range])
sig = std(qs[range])
xx  = linspace(mu-4sig, mu+4sig, 100)
plt.plot(xx, pdf(Normal(mu,sig),xx), "b--", linewidth = 3, label = "Normal distribution")
# plt.legend(loc="upper right",fancybox="true")

plt.subplot(212)
# plt.xlabel("exp(u)")
plt.ylabel("Normal(mu,sig)")
plt.title("Final state QQ plot")
res = qqbuild(qs[range],Normal(mu,sig))
plt.plot(res.qy, res.qx, linewidth = 3)
xy = [x for x in linspace(minimum([res.qx,res.qy]),maximum([res.qx,res.qy]),4)]
plt.plot(xy, xy, "r--", linewidth = 3)
plt.savefig("$dir/figure4$fname.png",transparent=true)
=#

## --------------------------------------------------------------------------------------------
## Plot parameters:
## --------------------------------------------------------------------------------------------


plt.figure(5)
plt.subplots_adjust(hspace=0.5)
plt.subplot(211)
plt.ylabel("K")
plt.title("K chain")
plt.plot(thetas[:,1], "r",  linewidth=2)
plt.plot(thetas[:,1], "ro", markersize=3)
plt.subplot(212)
plt.xlabel("t")
plt.ylabel("Gamma")
plt.title("Gamma chain")
plt.plot(thetas[:,2], "b",  linewidth=2)
plt.plot(thetas[:,2], "bo" ,markersize=3)
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
plt.plot(true_K, true_gam, "rs", markersize = 20, label = "True value", linewidth = 2)
# plt.legend(loc="lower right",fancybox="true")
plt.savefig("$dir/figure6$fname.png",transparent=true,dpi=300)


lowlim = maximum([2,nsample_burnin+1])                  # To possibly discard initial 
                                                        # system realizations (burn-in)
plt.figure(7)
plt.subplots_adjust(hspace=0.5)
plt.subplot(211)
plt.xlabel("K")
plt.ylabel("Probability density")
plt.title("PDF")
kd_K = KernelDensity.kde(thetas[lowlim:end,1])
plt.axis([-0.01,kd_K.x[end],0,maximum(kd_K.density)])
plt.hist(thetas[lowlim:end,1], 100, normed = 1, color = "y")
plt.plot(kd_K.x, kd_K.density, color = "r", linewidth = 3, label = "Kernel density kde(K)")
# mu_K   = mean(thetas[lowlim:end,1])
# sig_K  = std(thetas[lowlim:end,1])
xx_K   = linspace(0.0001, 400, 500)
plt.plot(xx_K, pdf(LogNormal(4,1),xx_K), "b", linewidth = 1, label = "LogNormal distribution")
# plt.axvline(x=mu_K, linewidth=2, color="r")
# plt.axvspan(mu_K-sig_K, mu_K+sig_K, facecolor="y", alpha=0.3)

plt.subplot(212)
plt.xlabel("Gamma")
plt.ylabel("Probability density")
plt.title("PDF")
kd_G = KernelDensity.kde(thetas[lowlim:end,2])
plt.axis([-0.01,kd_G.x[end],0,maximum(kd_G.density)])
plt.hist(thetas[lowlim:end,2], 100, normed = 1, color = "y")
plt.plot(kd_G.x, kd_G.density, color = "r", linewidth = 3, label = "Kernel density kde(Gamma)")
# mu_G   = mean(thetas[lowlim:end,2])
# sig_G  = std(thetas[lowlim:end,2])
xx_G   = linspace(0.0001, 1, 500)
plt.plot(xx_G, pdf(LogNormal(0,2),xx_G), "b", linewidth = 1, label = "LogNormal distribution")
# plt.axvline(x=mu_G, linewidth=2, color="r")
# plt.axvspan(mu_G-sig_G, mu_G+sig_G, facecolor="y", alpha=0.3)

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

plt.figure(9)
plt.subplots_adjust(hspace=0.2)
plt.subplots_adjust(wspace=0.5)
plt.subplot(221)
plt.ylabel("u[$u1]")
plt.title("u_bdy chain")
plt.plot(u_chains[:,1], "r",  linewidth=3)
plt.plot(u_chains[:,1], "ro", markersize=4)
plt.subplot(222)
plt.ylabel("u[$u2]")
plt.title("u_stg chain")
plt.plot(u_chains[:,2], "b",  linewidth=3)
plt.plot(u_chains[:,2], "bo" ,markersize=4)
plt.subplot(223)
plt.xlabel("t")
plt.ylabel("u[$u3]")
plt.plot(u_chains[:,3], "r",  linewidth=3)
plt.plot(u_chains[:,3], "ro" ,markersize=4)
plt.subplot(224)
plt.xlabel("t")
plt.ylabel("u[$u4]")
plt.plot(u_chains[:,4], "b",  linewidth=3)
plt.plot(u_chains[:,4], "bo" ,markersize=4)
plt.savefig("$dir/figure9$fname.png",transparent=true)

plt.figure(10)
plt.xlabel("t")
plt.ylabel("Beta")
plt.title("Beta chain")
plt.axis([0, nsample + 1, 0, 5])
plt.plot(bet, "g",  linewidth=2)
plt.plot(bet, "go", markersize=3)
plt.savefig("$dir/figure10$fname.png",transparent=true)

reject

############################################################################################
############################################################################################

#=
plt.figure(figsize=(12.5, 2.5))
axes()[:set_ylim]([0,1.8])
axes()[:set_xlim]([-10,850])
plt.axis("off")
# plt.xticks(size="16")
# plt.yticks(size="16") 
plt.plot(t[ty], y, "ro", markersize = 10)
# plt.plot(t,r,"bo",label="Rain input",markersize = 5, color="b")
# plt.errorbar(t[ty], r[ty], yerr=(2*sigma*r[ty],2*sigma*r[ty]), fmt="o", markersize = 12, color="b", capsize=6, elinewidth=2)
# plt.plot(t,S/true_K,"g",label="Output",linewidth=2)
# plt.plot(t,S_ex,"g",label="Output",linewidth=2)
# plt.errorbar(t[ty], y, yerr=(2*sigma*y,2*sigma*y), fmt="o", markersize = 10, color="r", capsize=10, elinewidth=3)
# plt.tick_params(length=6, width=2)
=#

#=
S_ex=vec(predy[1,:])
=#
#=
plt.figure(1)
plt.figure(figsize=(10.24, 2.56))
axes()[:set_ylim]([0,1.8])
plt.xticks(size="16")
plt.yticks(size="16") 
plt.plot(t,r,"bo",label="Rain input",markersize = 3, color="b")
# plt.errorbar(t[ty], r[ty], yerr=(2*sigma*r[ty],2*sigma*r[ty]), fmt="o", markersize = 12, color="b", capsize=6, elinewidth=2)
# plt.plot(t,S/true_K,"g",label="Output",linewidth=2)
plt.plot(t,S_ex,"g",label="Output",linewidth=2)
plt.errorbar(t[ty], y, yerr=(2*sigma*y,2*sigma*y), fmt="o", markersize = 10, color="r", capsize=10, elinewidth=3)
plt.tick_params(length=6, width=2)

range    = 1:N
redrange = iround(linspace(1, n+1, min(n+1,101))) # Number of data points that will be shown
                                                  # in the plot
plt.figure(2)
plt.xticks(size="15")
plt.yticks(size="15") 
for ip = 1:size(predy,1)
	plt.plot(range, vec(predy[ip,range]), label="Sampled system realizations", color="g", linewidth=1)
end
yerr=2*sigma*y[redrange]
plt.errorbar(ty[redrange], y[redrange], yerr=(yerr,yerr), fmt="o", markersize = 10, color="r", capsize=6, elinewidth=4)
plt.tick_params(length=5, width=2)


plt.figure(3)
axes()[:set_ylim]([0,2])
axes()[:set_xlim]([0,400])
plt.xticks(size="15")
plt.yticks(size="15") 
grid = iround(linspace(1, size(thetas)[1], 500))       
grid2 = grid[length(grid)-49:length(grid)]             
plt.plot(thetas[grid,1], thetas[grid,2], linewidth = 0.5)
plt.plot(thetas[grid,1], thetas[grid,2], "b.")
plt.plot(thetas[1,1], thetas[1,2], "go", markersize = 15, label = "Initial state")
plt.plot(thetas[grid2,1], thetas[grid2,2], "yo", markersize = 10, label = "Last points")
plt.plot(thetas[end,1], thetas[end,2], "ro", markersize = 15, label = "Last state")
plt.plot(true_K, true_gam, "rs", markersize = 20, label = "True value", linewidth = 2)
plt.tick_params(length=5, width=2)

lowlim = maximum([2,nsample_burnin+1])
plt.figure(4)
plt.xticks(size="15")
plt.yticks(size="15") 
kd_K = KernelDensity.kde(thetas[lowlim:end,1])
plt.axis([-0.01,400,0,maximum(kd_K.density)+0.002])
plt.hist(thetas[lowlim:end,1], 200, normed = 1, color = "y", alpha = 0.5)
plt.plot(kd_K.x, kd_K.density, color = "r", linewidth = 3, label = "Kernel density kde(K)")
plt.axvline(x=true_K, linewidth=3, color = "b")
plt.tick_params(length=5, width=2)

plt.figure(5)
plt.xticks(size="15")
plt.yticks(size="15") 
kd_G = KernelDensity.kde(thetas[lowlim:end,2])
plt.axis([-0.01,2,0,3])
plt.hist(thetas[lowlim:end,2], 80, normed = 1, color = "y", alpha = 0.5)
plt.plot(kd_G.x, kd_G.density, color = "r", linewidth = 3, label = "Kernel density kde(Gamma)")
plt.axvline(x=true_gam, linewidth=3, color = "b")
plt.tick_params(length=5, width=2)
=#
############################################################################################
############################################################################################