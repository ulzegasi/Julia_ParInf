## ============================================================================================
## DATA GENERATOR
##
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

dir  = "C:/Users/ulzegasi/Julia_files/ParInf_HMC" 
dir2 = "C:/Users/ulzegasi/Julia_files/ParInf_HMC/input_data" 


##
## ============================================================================================
## System parameters
## ============================================================================================
##
range = 2:5002
tdat  = float64(readdlm("$dir/t.dat")[range,2]) # Time points t.dat

n_temp = 10               # n+1 -> "end point" beads (see Tuckerman '93), n -> number of segments
N_temp = n_temp*1000 + 1  # Total number of points used to generate a system realization

t_temp  = linspace(tdat[1], tdat[end], N_temp)      # Generate N time points in [t[1],t[end]]      
T_temp  = t_temp[end]-t_temp[1]                     # Total time interval
dt_temp = T_temp/(N_temp-1)                         # Time step
ty_temp = iround(linspace(1, N_temp, n_temp + 1))   # Indeces of "boundary" beads (= measurement points)     

K     = 50.         # Retention time
gam   = 0.2         # Dimensionless noise parameter
sig   = 0.10        # Measurement noise

S = zeros(N_temp)   # System realization (with true parameters)

##
## ============================================================================================
## Rain input
## ============================================================================================
##
r = sin(t_temp/100.).*sin(t_temp/100.) + 0.1;    # Simple sinusoidal model                                                   
## Real rain input data: r=float64(readdlm("$dir/r.dat")[range,2])
##
## ============================================================================================
## Synthetic output discharge data 
## Use a system realization with true parameters (--> Ito) to generate synthetic discharge data 
## Measurement error model --> ln(y/r) = beta*q + sigma*epsilon
## N.B.: using the definitions S = (T g r / b^2)*exp(bq) and b = sqrt(Tg/K)
## ==> beta*q = ln(S/(kr)) = ln(y/r)
## ============================================================================================
##
S[1] = K * r[1]   # Define first point (--> <S*Peq(S)> with constant rain input)
for i = 2:N_temp              
    S[i] = S[i-1] + (dt_temp)*( r[i-1] - S[i-1]/K ) + sqrt(dt_temp)*sqrt(gam/K)*S[i-1]*randn()
    if  S[i] < 0 
        error("Negative volume!")
    end
end
bq = log((1/K)*(S[ty_temp]./r[ty_temp])) + sig*randn(n_temp+1)  
y  = r[ty_temp].*exp(bq)                                   

#=
plt.figure()
plt.xlabel("time")
plt.ylabel("r, S/K, y")
plt.title("System I/O")
plt.plot(t, r, "b--", linewidth=2)
plt.plot(t, S/true_K, "r", linewidth=2)
yerr=2*sigma*y
plt.errorbar(t[ty], y, yerr=(yerr,yerr), fmt="o", color="b", capsize=4, elinewidth=2)
=#

St      = zeros(N_temp,2)
St[:,1] = S
St[:,2] = t_temp

writedlm("$dir2/St_n10_K50_g02_s10_sinr.dat", St)
writedlm("$dir2/bq_n10_K50_g02_s10_sinr.dat", bq)
writedlm("$dir2/y_n10_K50_g02_s10_sinr.dat",  y)