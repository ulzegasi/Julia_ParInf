dir   = "C:/Users/ulzegasi/Julia_files/ParInf_HMC/good_results"
dir2  = "C:/Users/ulzegasi/Julia_files/ParInf_HMC/input_data"
fname = string("_dtau030_nnapa3_j30")
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

thetas   = readdlm("$dir/thetas$fname.dat")
bet      = sqrt(T*thetas[:,2]./thetas[:,1])

energies = readdlm("$dir/energies$fname.dat")
reject   = readdlm("$dir/reject$fname.dat")

io_data  = readdlm("$dir/iodata$fname.dat")
S = io_data[2:N+1]
r = io_data[(N+3):(2*N+2)]
y = io_data[(2*N+4):end]

# St = readdlm("$dir2/St$fname.dat")
# S  = St[:,1]
# tS = St[:,2]

u_chains = readdlm("$dir/u_chains$fname.dat")
u1       = int64(u_chains[1,1])
u2       = int64(u_chains[1,2])
u3       = int64(u_chains[1,3])
u4       = int64(u_chains[1,4])
u_chains = u_chains[2:end,:]

using PyPlot, Winston, ForwardDiff, Distributions, KernelDensity
using ReverseDiffSource
pygui(true)

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

K_init     = 200.0                              # Initial state
gam_init   = 0.5
bet_init   = sqrt(T*gam_init/K_init)

bq = Array(Float64,n+1)
for i = 1:n+1
	bq[i] = log(y[i]/r[ty[i]])
end
q = Array(Float64,N)
for s = 1:n
    q[ty[s]:ty[s+1]] = (1/bet_init)*linspace(bq[s],bq[s+1],ty[s+1]-ty[s]+1)
end

plt.figure(figsize=(12.5, 2.5))
axes()[:set_ylim]([-1.5,1.5])
axes()[:set_xlim]([-10,850])
plt.axis("off")
plt.plot(t[ty], q[ty], "ro", markersize = 18)

plt.figure(figsize=(12.5, 2.5))
axes()[:set_ylim]([-1.5,1.5])
axes()[:set_xlim]([-10,850])
plt.axis("off")
plt.plot(t, q, "o", markersize = 8, color = (1,0.65,0))
plt.plot(t[ty], q[ty], "ro", markersize = 18)

u = Array(Float64,N)
for s = 0:(n-1)
    u[s*j+1] = q[s*j+1]
    for k = 2:j
        u[s*j+k] = q[s*j+k] - ( (k-1)*q[s*j+k+1] + q[s*j+1] )/k
    end
end
u[N] = q[N]

mp      = Array(Float64,nparams+N)     
p       = Array(Float64,nparams+N)


mp[1:nparams] = m_theta
for s = 1:n  
    mp[nparams+(s-1)*j+1] = m_bdy                       
    for k = 2:j 
        mp[nparams+(s-1)*j+k] = m_stg                   
    end
end
mp[nparams+N] = m_bdy          

p = sqrt(mp).*randn(nparams+N)

params = nparams
theta = [bet_init, gam_init]
r = Array(Float64,N)
r = sin(t/100.).*sin(t/100.) + 0.1;
lnr_der = Array(Float64,N)                      
for i = 1:(N-1)  lnr_der[i] = ( log(r[i+1]) - log(r[i]) ) / dt  end
reload("C:/Users/ulzegasi/Julia_files/ParInf_HMC/ParInf_Fun_AD_napa.jl")

H_old = V_N_fun(u) + V_n_fun(theta,u) + V_1_fun(theta,u) + sum((p.^2)./(2*mp))

time_respa      = Array(Float64,nsample_eff)
time_respa_s    = zeros(nsample_eff)
time_respa_f    = zeros(nsample_eff)
time_respa_f    = zeros(nsample_eff)
counter = 1
for counter_napa = 1:10 
    napa(theta, u, counter) 
end

# Back transformations (u -> q) to update the {q} variables :
    # (Tuckerman et al., JCP 99 (1993), eq. 2.19)
    # ---------------------------
    ss=n-1
    while ss>=0
        k=j+1
        q[ss*j+k] = u[ss*j+k]
        while k>2
           k -= 1
           q[ss*j+k] = u[ss*j+k] + u[ss*j+1]/k + (k-1)*q[ss*j+k+1]/k
        end
        ss -= 1
    end
    q[1] = u[1]





plt.figure(figsize=(12.5, 2.5))
axes()[:set_ylim]([-1.5,1.5])
axes()[:set_xlim]([-10,850])
plt.axis("off")
plt.plot(t, q, linewidth=4, color = (1,0.0,0))
plt.plot(t, q, "o", markersize = 8, color = (1,0.65,0))
plt.plot(t[ty], q[ty], "ro", markersize = 18)



counter = 2
for counter_napa = 1:10 
    napa(theta, u, counter) 
end

# Back transformations (u -> q) to update the {q} variables :
    # (Tuckerman et al., JCP 99 (1993), eq. 2.19)
    # ---------------------------
    ss=n-1
    while ss>=0
        k=j+1
        q[ss*j+k] = u[ss*j+k]
        while k>2
           k -= 1
           q[ss*j+k] = u[ss*j+k] + u[ss*j+1]/k + (k-1)*q[ss*j+k+1]/k
        end
        ss -= 1
    end
    q[1] = u[1]

plt.figure(figsize=(12.5, 2.5))
axes()[:set_ylim]([-1.5,1.5])
axes()[:set_xlim]([-10,850])
plt.axis("off")
plt.plot(t, q, linewidth=4, color = (1,0.0,0))
plt.plot(t, q, "o", markersize = 8, color = (1,0.65,0))
plt.plot(t[ty], q[ty], "ro", markersize = 18)


for counter = 3:100
	for counter_napa = 1:10 
	    napa(theta, u, counter) 
	end

	# Back transformations (u -> q) to update the {q} variables :
	    # (Tuckerman et al., JCP 99 (1993), eq. 2.19)
	    # ---------------------------
	    ss=n-1
	    while ss>=0
	        k=j+1
	        q[ss*j+k] = u[ss*j+k]
	        while k>2
	           k -= 1
	           q[ss*j+k] = u[ss*j+k] + u[ss*j+1]/k + (k-1)*q[ss*j+k+1]/k
	        end
	        ss -= 1
	    end
	    q[1] = u[1]

	plt.figure(figsize=(12.5, 2.5))
	axes()[:set_ylim]([-1.5,1.5])
	axes()[:set_xlim]([-10,850])
	plt.axis("off")
	plt.plot(t, q, linewidth=4, color = (1,0.0,0))
	plt.plot(t, q, "o", markersize = 8, color = (1,0.65,0))
	plt.plot(t[ty], q[ty], "ro", markersize = 18)
	figindex = counter+2
	plt.savefig(string("C:/Temp/figure_",figindex,"q_stg"),transparent=true)
end




plt.figure(figsize=(6.5, 6.5))
axes()[:set_ylim]([0,1.5])
axes()[:set_xlim]([0,250])
grid = fullgrid[1:1]                   
plt.plot(thetas[grid,1], thetas[grid,2], linewidth = 0.5)
plt.plot(thetas[grid,1], thetas[grid,2], "bo", markersize = 6)
plt.plot(thetas[1,1], thetas[1,2], "ro", markersize = 18, label = "Initial state")
plt.plot(true_K, true_gam, "go", markersize = 18, label = "True value", linewidth = 2)
plt.tick_params(length=5, width=2)

fullgrid = iround(linspace(1, size(thetas)[1], 120)) 
for counter = 2:120
	plt.figure(figsize=(6.5, 6.5))
	axes()[:set_ylim]([0,1.5])
	axes()[:set_xlim]([0,250])
	grid = fullgrid[1:counter]                   
	plt.plot(thetas[grid,1], thetas[grid,2], linewidth = 0.5)
	plt.plot(thetas[grid,1], thetas[grid,2], "bo", markersize = 6)
	plt.plot(thetas[1,1], thetas[1,2], "ro", markersize = 18, label = "Initial state")
	plt.plot(true_K, true_gam, "go", markersize = 18, label = "True value", linewidth = 2)
	plt.tick_params(length=5, width=2)
	plt.savefig(string("C:/Temp/figure_",counter,"_2dchain"),transparent=true)
end

fullgrid2 = iround(linspace(1, size(thetas)[1], 5000))
plt.figure(figsize=(6.5, 6.5))
	axes()[:set_ylim]([0,1.5])
	axes()[:set_xlim]([0,250])
	grid = fullgrid2                   
	plt.plot(thetas[grid,1], thetas[grid,2], linewidth = 0.5)
	plt.plot(thetas[grid,1], thetas[grid,2], "bo", markersize = 6)
	plt.plot(thetas[1,1], thetas[1,2], "ro", markersize = 18, label = "Initial state")
	plt.plot(true_K, true_gam, "go", markersize = 18, label = "True value", linewidth = 2)
	plt.tick_params(length=5, width=2)


plt.figure(figsize=(12.5, 3.5))
axes()[:set_ylim]([0,250])
axes()[:set_xlim]([0,50000])
plt.plot(thetas[:,1], "r",  linewidth=2)
plt.plot(thetas[:,1], "ro", markersize=3)

plt.figure(figsize=(12.5, 3.5))
axes()[:set_ylim]([0,1.8])
axes()[:set_xlim]([0,50000])
plt.plot(thetas[:,2], "b",  linewidth=2)
plt.plot(thetas[:,2], "bo", markersize=3)

range    = 1:N
redrange = iround(linspace(1, n+1, min(n+1,101))) # Number of data points that will be shown
                                                  # in the plot
plt.figure(figsize=(12.5, 6.5))
plt.xticks(size="15")
plt.yticks(size="15") 
axes()[:set_ylim]([0,3])
axes()[:set_xlim]([0,310])
for ip = 1:size(predy,1)
	plt.plot(range, vec(predy[ip,range]), label="Sampled system realizations", color="g", linewidth=1)
end
yerr=2*sigma*y[redrange]
plt.errorbar(ty[redrange], y[redrange], yerr=(yerr,yerr), fmt="o", markersize = 18, color="r", capsize=6, elinewidth=8)
plt.tick_params(length=5, width=2)


lowlim = maximum([2,nsample_burnin+1])
plt.figure(figsize=(8.5, 6.5))
plt.xticks(size="15")
plt.yticks(size="15") 
kd_K = KernelDensity.kde(thetas[lowlim:end,1])
plt.axis([-0.01,200,0,maximum(kd_K.density)+0.002])
plt.hist(thetas[lowlim:end,1], 100, normed = 1, color = "y", alpha = 0.5)
plt.plot(kd_K.x, kd_K.density, color = "r", linewidth = 3, label = "Kernel density kde(K)")
plt.axvline(x=true_K, linewidth=3, color = "b")
plt.tick_params(length=5, width=2)

plt.figure(figsize=(8.5, 6.5))
plt.xticks(size="15")
plt.yticks(size="15") 
kd_G = KernelDensity.kde(thetas[lowlim:end,2])
plt.axis([-0.01,1.2,0,4.0])
plt.hist(thetas[lowlim:end,2], 80, normed = 1, color = "y", alpha = 0.5)
plt.plot(kd_G.x, kd_G.density, color = "r", linewidth = 3, label = "Kernel density kde(Gamma)")
plt.axvline(x=true_gam, linewidth=3, color = "b")
plt.tick_params(length=5, width=2)