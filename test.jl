
x=1.3
typeof(x)

int64('b')

for i=97:110
	print(string(char(i), "\n"))
end
x=(1,2,3)
typeof(x)

x[1],x[3]

x[4]=44

a,b=x

a=x

a

x,y=a

x

y


[1 2 3 4 5; 11 12 13 14 15; 1 2 3 4 5; 7 8 9 1 2; 2 4 6 8 1]$

q=Array(Int64,2)

rand(2,2)

[[1 2], [3 4]]

 Int64[vec([1 2] [3 4])]

 [x*y for x in 1:2, y in 1:4]

[[1 2] [3 4]]

w=vec([1 2 3; 4 5 6; 7 8 9])

[1 2 3; 4 5 6; 7 8 9]

rt=reshape(w,3,3)

sum(rt)
mean(rt)
std(rt)
std(rt)/mean(rt)*100

count(x->x>=8,rt)

rt

a=[x*y for x in 1:2, y in 1, z in 1:3]

a[7]

rt[2,3]

rt[end]

rt[1:2:end]

asa=reshape(1:9,3,3)

asa[1:2:end]

linspace(1,2,5)

a=linspace(0,1)

a[repmat([true, false],50)]


x=reshape(1:8,2,4)

x[:,2:3]=repmat([1 3],2)

x

repmat([1,3],2)

x=cell(2)

x[1]=ones(2)
x[2]=trues(3)

x

a=x

a

b=copy(x)
c=deepcopy(x)

x[1]="Bang"

x[2][1]=false

x
a

b
c

p=Point(0,0.0,"Origin")

typeof(p)

p.x=1
names(p)

names(Point)

p.x


y={"a"=>1,(2,3)=>true}

y["a"]
y[(2,3)]

haskey(y,"a")


keys(y),values(y)

typeof(keys(y))

repr(123.3)

x=123
asas="$x+3=$(x+3)"
typeof(asas)

x==123 ? a=23 : a=44

i=0;
while true
	i += 1
	if i > 10
		break
	end
	print(string(i,"\n"))
end

for x in 1:10
	if 3 < x < 6
		continue
	end
	println(x)
end


g(x::Int64, y::Bool)=x*y

g(1,true)

asad=true
typeof(asad)

int64(true::Bool),int64(false::Bool)



methods(g)




h(x...)=sum(x)/length(x)-mean(x)/2
h(1,2,3)
x=(2,3)

h(x)

qq(f::Function,x)=2*f(x)
qq(x->2*x,10)

print(asa)

mm=reshape(1:12,3,4)

map(x->x^2,mm)

qqq(10) do x
2*x
end


f() = global t = 1

f(),t

Fs=cell(2)
i=1
j=1
Fs[1]=()->j
i=i+1
j=i
Fs[i]=()->j

Fs[1]()

asa=cell(2)
i=1
while i <= 2
	j=i
	asa[i]=j
	i += 1
end
asa[2]

module M
export x
x=1
y=2
end

using M
x

y


[1 2] .< [2 1]
[1 2] * [4, 5]

 Pkg.update()
using DataFrames
 df=DataFrame(A=1:4,B="a")
df[:C]=repmat([true, false],2)
tail(df)
df[1:2,[:A,:C]]

using PyPlot
pygui(true)
x=linspace(0,2*pi,1000)
y=sin(2*pi*x) .* exp(-1*x)
plot(x,y)
grid(true)
axis(true)



help(linspace)


length(x)
using Distributions
srand(1)

x=randn(10000);
bins=20;
plt.hist(x,50,normed=1,facecolor="y")
points=linspace(bins[1],bins[end])
plot(points,pdf(Normal(),points),"r")


v=[randn(10000)]
nbins=10
fig=figure()
fig=fig[:canvas][:draw]()
ax=axes()
ax[:set_ylim]([0,1]),ax[:set_xlim]([-5,5])
h=PyPlot.plt.hist(v,nbins,normed=1)



typeof(range)
range[1]
length(range)
asa=readdlm("C:/Users/ulzegasi/Julia_files/Julia_1/t.dat")
typeof(asa)
length(asa)
size(asa)
10657*2

typeof(asa)


a=reshape(1:12,3,4)
a[2:3,2]


asa=readdlm("C:/Users/ulzegasi/Julia_files/Julia_1/t.dat")[range,2]


aspect_ratio()
using PyPlot
pygui(true)
plot(asa)

diff=[t[i+1]-t[i] for i in 1:(length(t)-1)] 

dt
N

T


dt
dt2=(t[end]-t[1])/(N-1)
diff[1]
dt2-diff[1]
dt-diff[1]

typeof(dt)

diff-dt

diff[end]-diff[1]



t2 = float64((readdlm("C:/Users/ulzegasi/Julia_files/Julia_1/t.dat")[range,2]))

typeof(t2)

range = 2:5002

typeof(range)

typeof(range[1])

t  = float64(readdlm("C:/Users/ulzegasi/Julia_files/Julia_1/t.dat")[range,2])

length(t)

N  = length(t)                         # Number of discrete time points
T  = t[end]-t[1]                       # Total time
dt = T/(N-1)	                       # Time step


diff=[t[i+1]-t[i] for i in 1:(length(t)-1)]

diffdiff=[diff[1]-diff[i] for i in 1:length(diff)]

T/(N-1)
T/N
(N-1)/j
typeof((N-1)/j)


linspace(1, N, n+1)

asa=linspace(1, 5001, n+1)

isinteger(asa[1])==true

for i in 1:length(asa)
	if isinteger(asa[i])==false
		println(i,"\n")
	end
end
i=5
print(i)
isinteger([1.0, 11.0])

int64(90.1)
asa[10]

for i in 1:length[asa]
	if isinteger(asa[i])==false
		print(i)
	end
end



N

Pkg.add("Winston")
Pkg.update()


r = zeros(N);
typeof(r)

r[1:5] = 22;

asas = sin(t/100).^2 + 0.1

typeof(asas)
typeof(r)

asas==r
plot(r)
plot(r)
 length(t)

using PyPlot
pygui(true)

plt.plot(t,r)

x=linspace(0,1)
y = sin(4*pi*x).*exp(-5*x)
plot(x,y)

lnr_der

plt.plot(lnr_der)

length(lnr_der)

theta

K/T

length(r)

plt.plot(S)

srand(1)
asa = zeros(10)
for i = 1 : 10
    asa[i] = randn()
end

asa

srand(1)
asa2 = zeros(10)
for i = 1 : 10
    asa2[i] = randn()
end

asa2
S[1]
K
r[1]


plt.plot(u)
plt.plot(t, S[1:N])
plt.plot(t,r)

plt.plot([1,2,3,4])
plt.ylabel('some numbers')

t = 0:0.2:5
length(t)
t[20]
plt.figure(1)
plt.subplot(211)
plt.plot(t,"b-",t.^2,"ro",t.^3,"g--")
plt.subplot(212)
plt.plot(t.^2,"ro",t.^4,"g--")


plt.figure(2)
plt.subplot(211)
plt.plot(t,"b-",t.^2,"ro",t.^3,"g--")
plt.subplot(212)
plt.plot(t.^2,"ro",t.^4,"g--")

asa=[11,22,33,45,56,69,73,89,91,108]
pt=[1,3,6,9]
asa[pt]

name = input("What's your name? ");

print("enter year  "); year = int(readline(STDIN))

readline(STDIN)as

println()
print("enter year  ") 
year = int(readline(STDIN))

date = Date(year, month, day)
println(date)

year

qq  = float64(readdlm("$dir/t.dat")[range,2])


kk = 5
ss = 6
testname = string("kk$kk","_ss$ss")
testname

mp[1:12]
mp[13:24]

asa = function(theta) end

theta

asa

fnn            = function(theta) V_fast(theta,u) end

function(theta) V_fast(theta,u) end

theta

fnn            = function(theta) V_fast(theta,u) end
V_fast_bytheta = forwarddiff_gradient(fnn, Float64, fadtype=:dual, n=s)

help(V_fast_bytheta)

V_fast_der(theta,u)

V_fast_bytheta(theta)

f(x) = x[1]^2+x[1]*x[2]^3

# Using forwarddiff_gradient
g = forwarddiff_gradient(f, Float64, fadtype=:dual, n=2)
g([2., 2.])

V_fast(theta,u)
x=2; y=3
qq = [x, y]
z=2
function asafunfun(a,b)
	out = a[1]^2*a[2]+2*a[2]^2*a[1]^3+b*a[1]*a[2]
	return out
end

asafunfun(qq,z)

asa1= qq[1]^2+2*qq[2]+z

asafunfun2 = function(qq) asafunfun(qq,z) end 
qq2=[2,1]

asafunfun2(qq2)

using ForwardDiff

gasa2 = forwarddiff_gradient(asafunfun2, Float64, fadtype=:dual, n=2)

gasa2([3.0, 5.5])

function asafunfun3(qq)
	asafunfun(qq,z)
end

asafunfun3(qq2)

gasa3 = forwarddiff_gradient(asafunfun3, Float64, fadtype=:dual, n=2)

gasa3([3.0, 5.5])

function test(xx)
	out = xx^2
	return out
end

function test2(xx)
	out = xx^3
	return out
end

test(4), test2(4)

srand(28101973)

randn(100)

randn(153)

a1=randn()

srand(28101973)
randn(100)
randn(153)

a2=randn()

a1==a2

srand(28101973)

randn(102)
randn(151)
a3=randn()

a2==a3

datetime(today())

time()

t1=time()

t2=time()

time()/3600/24/365

ty

asa1=time()

asa2=time()

asa2-asa1

for i=1:10
	println(i)
end

2*20000/3600

(0.119-0.11)/0.11

mptest

mptest2 = zeros(2+10)
mptest2[1:2] = 100.

mptest2

for i=1:10 println("asa") end

2.5*20000/3600

0.02/15

function asa(x,y)
	global x = x+2
	global y = -y 
end

x=[1,2]
y=[44,45]

asas(x,y)

x
y

theta

RESPA(theta,u,p,mp,dtau,nn)

theta

for counter_respa = 1:n_respa 
        RESPA(theta,u,p,mp,dtau,nn)
    end 
theta

theta
reload("$dir/1D_functions2_2.jl")
asa(theta)
theta
asa(theta)
theta
typeof(x)
x=[1.,2.]
y=[3.,4.]

asa2(x)
x
a=[1,2]
sd=[1,2]
sd=sd+1
asa(a)
a

t0=time()
for counter=1:100
	if (counter%10 == 0)
		println(counter, " is ", time(), " asa\n")
	end
end

nasa=10000000
v=zeros(nasa)
for i=1:4
	for j=1:nasa 
		v[j]=j^2+1/j^3-j*(j^4-1)
	end
	print(string("-> ", i, " done\n"))
end
ypred[:,1:4]

asa = [1.0, 2.0, 3.0]

asa::Array(Float64)

theta

q

V_slow(theta,q)

V_fast(theta,q) + V_slow(theta,q) + sum((p .* p) ./ (2*mp))
q_sample

u = Array(Float64,N)
for i=0:(n-1)
    u[i*j+1] = q[i*j+1]
    for k=2:j
        u[i*j+k] = q[i*j+k] - ( (k-1)*q[i*j+k+1] + q[i*j+1] )/k
    end
end
u[N] = q[N]

u

tasa=time()

    tmp=0.
    tmp=q[1]^2+q[N]^2-2.0*q[1]*q[2]
    out=0.
    for ii = 2:N-1
        tmp += 2.0*(q[ii]^2-q[s]*q[s+1])
    end
     out += K*tmp/(2.0*gamma*dt)
    
    # Data:
    # ---------------------------

    for ii=0:n
        out += ( y[ii+1] - r[ii*j+1]*exp(q[ii*j+1]) )^2/(2*sigma^2)
    end
tasa2=time()

tasa2-tasa

reload("$dir/1D_functions2_2_SU.jl")
tasa=time()
RESPA(theta, q, p, mp, dtau, nn)
tasa2=time()

tasa2-tasa

reload("$dir/1D_functions2_2_SU.jl")
timeq=zeros(100)
for i=1:100
	tasa=time()
	RESPA(theta, q, p, mp, dtau, nn)
	tasa2=time()
	timeq[i]=tasa2-tasa
end
timeq
mean(timeq)
reload("$dir/1D_functions2_2.jl")
timeu=zeros(100)
for i=1:100
	tasa=time()
	RESPA(theta, u, p, mp, dtau, nn)
	tasa2=time()
	timeu[i]=tasa2-tasa
end
timeu
mean(timeu)

tmp=0.
    tmp=q[1]^2+q[N]^2-2.0*q[1]*q[2]
    out=0.
     for ii=1:n
        tmp1 = 0.
        for k = 2:j
            tmp1 += 2.0*(q[(ii-1)*j+k]^2-q[(ii-1)*j+k]*q[(ii-1)*j+k+1])
        end
        if ii == 1 
            tmp += tmp1
        else 
            tmp += 2.0*( q[(ii-1)*j+1]^2 - q[(ii-1)*j+1]*q[(ii-1)*j+2]) + tmp1
        end        
    end
   


    #for s = 2:N-1
    #    tmp += 2.0*(q[s]^2)#-q[s]*q[s+1])
    #end
    out += K*tmp/(2.0*gamma*dt)
tmp
out
q[1]^2+q[N]^2-2.0*q[1]*q[2]

q[1]

q
for ii=2:(N-1)
	tmp += 2.0*(q[ii]^2-q[ii]*q[ii+1])
end
out += K*tmp/(2.0*gamma*dt)
ypred
writedlm("$dir/ypredu.dat", ypred)
asau=readdlm("$dir/ypredu.dat")
asaq=readdlm("$dir/ypredq.dat")


asadiff=asaq-asau

plt.plot(vec(asaq[39,ty]),"bo")
plt.plot(vec(asau[39,ty]),"r^")


vec(asaq[39,ty])

maximum(asadiff./asau*100)

trest=[1 2 3; 4 5 6; 7 8 9]
trest2=[2 2 2; 3 3 3; 1 1 1]

trest./trest2

asatest=[1 2 3; 4 5 6; 7 8 9]

plt.plot(vec(maximum(asatest,1)))

xvec=randn(1000)

plt.plot(kde(xvec))

X = rand(Normal(),10);
kk = KernelDensity.kde(X)
plt.plot(kk.x, kk.density, color = "y", linewidth = 3.0)

display(plt.plot(kk))

x = 0:0.1:10
y = sin(x)
y2 = sin(2sin(2sin(x)))
plt.plot(x, y, "g^", x, y2, "b-o")

kdeboot = KernelDensity.kde(bootcor)
# plot results
plt.hist(bootcor, 50, normed = 1)
plot(kdeboot.x, kdeboot.density, color = "y", linewidth = 3.0)
axvline(0.5, color = "r", linewidth = 3.0)
savefig("corboot.pdf", format = "pdf") # save results to pdf

using ReverseDiffSource

function testasa(x,y)
	out=x^3+y^2+x*y
	return out 
end
testasa2(x)=x[1]^3+x[2]^2+x[1]*x[2]
asaout=rdiff(testasa2, (ones(2),)) 
testasa(0,2)
rdiff(:(testasa(x,1)),x=2.)

asaout(2)
asa
function funfun(x,y)
	x[1]^2+x[2]^3+y[1]+3*y[2]
end
x=[2,3]
funfun([1,4],[2,2])
asa2=[0,1]

function funtest2(asa)
	funfun(asa,asa2)
end
funtest3(asa) = funfun(asa, asa2)
funtest4(asa) = asa[1]^2+asa[2]^3+asa2[1]+3*asa2[2]
funtest      = function(asa) funfun(asa,asa2) end

funtest_der  = forwarddiff_gradient(funtest, Float64, fadtype=:dual, n=2)

funtest_der([1.0,3.0])
funtest4([2,1])
funtest2
rosenbrock
funtest3
dai = rdiff(funtest4, (ones(2),))
rosenbrock([1,0])
rosenbrock(x) = (1 - x[1])^2 + 100(x[2] - x[1]^2)^2+ hh
hh=87
rosen2 = rdiff(rosenbrock, (ones(2),))

rosen2([1,2])

using ReverseDiffSparse
asa = placeholders(2)
ex = @processNLExpr asa[1]^2+asa[2]^3+asa2[1]+3*asa2[2]
fg = genfgrad_simple(ex)

function asa(x,y)
	x[1]^2+2*x[2]+y 
end

asa([1,2],1)

asared      = function(x) asa(x,y) end 
asaredder  = forwarddiff_gradient(asared, Float64, fadtype=:dual, n=2)
y=3
asaredder([2.,2098.])


asaredder2=forwarddiff_gradient((asared2=function(x) asa(x,y) end), Float64, fadtype=:dual,n=2)



(86-71)/86

20000/30*71/3600

(79-71)/79

predyold=readdlm("$dir2/predyOLD.txt");
predynew=readdlm("$dir2/predyNEW");

predyold-predynew

asa1=[1 2 3; 4 5 6; 7 8 9]
asa2=[1 2 3; 4 5 6; 6 7 5]
asa1-asa2

predydiff=vec(predyold-predynew)
length(predydiff)

ind0=0
non0vec=Array(Float64,length(predydiff))
for i=1:length(predydiff)
	if(predydiff[i] != 0 && predydiff[i] > 0.000000000000001 )
		ind0 += 1
	end
end
ind0

maximum(predyold)
minimum(predyold)

maximum(predynew)
minimum(predynew)
y=[0,0]
function asa2(x,y)
	qw = x[1]^2+2*x[2]+y[1]+2*y[2]
	return qw
end

asa([1,2],[2,3])

asared2 = function(x) asa2(x,y) end
y=[2,3]
asared2([2,0])

y=[5,5]

asared([2,0])

for i=1:5
	y=[i+1,0]
	println(y[1]," ",2*i," ",asared2([0,i]))
end

15*3600*0.1

15*3600-5400

48600/3600

M_theta 
M_bdy   
M_stage 


plt.figure(12)
plt.xlabel("K")
plt.ylabel("Gamma")
plt.title("Evolution in K-Gamma space")
grid  = iround(linspace(1, size(theta_sample)[1], 1000))	
grid2 = grid[length(grid)-49:length(grid)]
#plt.plot(thetas[grid,1], thetas[grid,2], linewidth = 0.5)
#plt.plot(thetas[grid,1], thetas[grid,2], "b.")
#plt.plot(thetas[1,1], thetas[1,2], "go", markersize = 10, label = "Initial state")
plt.plot(thetas[grid2,1], thetas[grid2,2], "yo", markersize = 10, label = "End points")
#plt.plot(thetas[end,1], thetas[end,2], "ro", markersize = 10, label = "Last state")
plt.plot(true_K, true_gamma, "rx", markersize = 15, label = "True value")


asa=[[1,2],[3,4]]

x=Array(Int64,10)
a=0
for i = 1 : 10
	a=i
	x[i]=a
end

x
a
qq=4
function asa(x,y)
	y[1] = y[1] + qq
	y[2] = y[2] + 3
	a = x^2
	return(a)
end

y=[4,3]

asa(2,y)
y

3+
4+

qw=65
typeof(qw)
asa=Array(Float64,2)
asa=[qw,728]

quid = ["a", "b", "c"]
quod = [5.6, 6.7, 7.8]

writedlm("C:/Users/ulzegasi/Julia_files/quidquod.dat", hcat(quid,
	quod))

writedlm("C:/Users/ulzegasi/Julia_files/params", hcat(param_names,param_values))

asa=Array(Float64,6,3)
for r = 1:6
	for c = 1:3
		asa[r,c] = r*c
	end
end

asa2=Array(Float64,3,3)
r2=1
for r = 2:2:6
	asa2[r2,:] = asa[r,:]
	r2 += 1
end
asa2


asas = ["asa1", "asa2"]

axa=repr(123)

asa=["N", "j", "n", "sigma", "dt", "true_K", "true_gamma", "nsample_burnin", "nsample_eff",
 "dtau", "nn", "n_respa", "M_bdy_burnin", "M_theta_burnin", "M_stage_burnin", "M_bdy", "M_theta", "M_stage"]

 match("N", asa)

 str="A simple string"

typeof(str)
typeof(asas[1])
asas[1][1]

for i in (1:length(asa))
	println(i, " ", asa[i], "\n")
	if(asa[i] == 'N')
		println("Found!")
		N = i
	end
end
N


for i = 1:length(asa)
	if(asa[i] == "j")
		j=i
	end
end 

asa[2]=="j"

j=0
for ind = 1 : length(asa)
	println("ind = ", ind, " --- asa[", ind,"] = ", asa[ind], " --- = j? ", (asa[ind]=="j"))
	println("=======******** ")
	if (asa[ind] == "j")
		j = ind
	end 
end

post=findin(asa,["j"])
post

param_names = vcat("N", "j", "n", "sigma", "dt")
param_values = vcat(3, 7, 44, 61, 5)

bamboo = hcat(param_names,param_values)

N = bamboo[findin(bamboo[:,1],["N"]),2]
j = bamboo[findin(bamboo[:,1],["j"]),2]
n = bamboo[findin(bamboo[:,1],["n"]),2]
sig = bamboo[findin(bamboo[:,1],["sigma"]),2]

for i = 1:4; println(i); end

esuna =vcat("a", "b", 
	"c", "d")
mier1 = Array(Int64, 2); mier1 = vcat(1, 2)
mier2 = Array(Float64, 2); mier2 = vcat(3.0, 4.0)
mier  = Array(Any, 4); mier = vcat(mier1,mier2)

asa = vcat([1 2 3])

size(asa,1)

function outf(asa)
	f = open("C:/Users/ulzegasi/Julia_files/hello.txt","a")
	write(f,asa)
	close(f)
end

outf(string("\n","Hallo world!"))

IOStream

f = open("C:/Users/ulzegasi/Julia_files/hello.txt","a")
f
@printf(f,"First line\n")
@printf(f,"Second line\n")
writedlm(f,string("First line\n"),"\n"); flush(f)

writedlm(f,string("\nSecond line\n"),"\n"); flush(f)

flush(f)

open(println(string("\nFirst line\n")),f)

close(f)

counter = 1868
p = sqrt(mp).*randn(s+N)
H_old = V_fast(theta,u) + V_slow(theta,q) + sum((p .* p) ./ (2*mp))
energies
energies[counter] = H_old

for i = 1:s   theta_save[i] = theta[i]    end
for i = 1:N   u_save[i]     = u[i]        end

deltatau = dtau/nn 
force = -V_slow_der(theta,q)
for i=1:(N+s) 
    p[i] = p[i]  + (dtau/2)* force[i]
end
p


mp = Array(Float64, 2+10)

asa=3.0*[1.0, 3.0]

mp[1:2]=asa

mp

mp = Array(Float64, 2+10)

asa = 3

asa -= 5+2

asa  = float64(vec(readdlm("C:/Users/ulzegasi/SWITCHdrive/JuliaTemp/Testing_s5_n10_j30_k20/kln_s5_n10_j30_k20_nn32.dat")))

asa[3]

asa = 1


asa -= 3-(7+1)

beta = 3
0.5*beta^2

asa = 4 + 2 -
7 

asa

reload("C:/Users/ulzegasi/Julia_files/funtest.jl")

asa=[1.0, 2.0, 3.0]

qw = Array(Float64, 3)

qw = asa

qw

ws = 1.5

typeof(ws)

test()

qw

function asa(q,w)
	xx = q[1]
	yy = q[2]
	xx2 = 2*xx
	yy2 = 3*yy
	out=q[1]^2+xx2^2+2*q[2]^3-yy2^2+w^2
	return out
end

q=[0.0,1.0]
asa(q,1)
using ForwardDiff
asaq = function(q) asa(q,w) end
asaqder = forwarddiff_gradient(asaq, Float64, fadtype=:dual, n=2)

methods(asaqder)
typeof(q)
w=1
asaqder([3.0,4.0])


addprocs(4)

wk=workers()

r = remotecall(2,rand,2,2)

fetch(r)
rmprocs(2,3,4,5)
for i=1:20  
   # spawn the process  
   r = @spawn rand(Int64)  

   # print out the results  
   @printf("process: %d init_process: %d id: %d int: %f\n", r.where, r.whence, r.id, fetch(r))  
 end  
3+3
addprocs(4)

nprocs()
nworkers()

# get this process's identifier  
 id = myid()  

 # create variable, store random number  
 rand_num = rand(Int64)  
 @printf("process: %d   int: %d", id, rand_num)  

spawn()
methods(spawn)

wk[2]

fetch(r)

r = @spawn rand(Int64)

rmat=remotecall(5,rand,2,2)

fetch(rmat)

rmat2 = remotecall_fetch(5, rand, 2,2)

rmat2

m1 = @spawnat 4 rand(2,2)
m2 = @spawn rand(3,3) 

fetch(@spawnat 2 rand(5,5))
@fetchfrom(2, rand(5,5))

m1

remotecall_fetch(2,getindex,r,1,1)

asa=remotecall(2,getindex,r,1,1)

fetch(asa)

fetch(r)

r[1,1]

r

message = "This string is constructed locally"

typeof(message)

shouting = @spawnat 1 uppercase(message)

fetch(shouting)

big_matrix = ones(10000,10000)
@spawn big_matrix .+ pi 

fetch(RemoteRef(2,1,65))
pi+0.463795

rand(100,2)
asa=[[1, 2, 3] [4, 5, 6]]

help(rand)

asa=asa.^2
asa[:,1]

sum(asa[:,1]+asa[:,2])

17+29+45

asa=rand(20,2).^2

numinc=sum(asa[:,1]+asa[:,2] .< 1)

for i=1:20
	println(asa[i,1]+asa[i,2] < 1)
end

asa[6,1]+asa[6,2]

tic()
numIter = 200000000
rand_coords = rand(numIter,2).^2
num_in_circle = sum(rand_coords[:,1]+rand_coords[:,2].<1)
ans = 4 * num_in_circle / numIter
toc()
print(ans)

ans=0
rand_coords=0
gc()

nprocs()

addprocs(2)
nprocs()
@everywhere function circle_number(numIter)
	rand_coords = rand(numIter,2).^2
	return(sum(rand_coords[:,1]+rand_coords[:,2].<1))
end
numIter = 200000000
tic()
numIter2=convert(Int64,numIter/4)
M = {numIter2 for i = 1:4}
sum_answer = pmap(circle_number,M)
ans = sum(sum_answer) / numIter2
toc()
print(ans)

sum_answer

rmprocs(6)

procs()
t1=time()
t2=time()

t2-t1

addprocs(6)

for i = 2:nprocs()
	rmprocs(procs()[i])
end

procs()

procs()[2]

rmprocs(8,9)

2:nprocs()

int(sum(rand(2).^2)<1)

int(true)

sum([true, true, true])
nprocs()
procs()


addprocs(2)
using DummyModule
require("C:/Users/ulzegasi/Julia_files/testfun.jl")
@everywhere 
rmprocs(2,3)
f(2)



procs()

workers()

asa=MyType(7)

typeof(asa)

reload("C:/Users/ulzegasi/Julia_files/MyModule.jl")

@everywhere using MyModule

x()

addprocs(2)

remotecall_fetch(2, ()->x())

a = 44
tic()
@parallel for i = 1:1000000000
	a^2+i^3-a*i
end
toc()


reload("C:/Users/ulzegasi/Julia_files/count_heads.jl")

addprocs(2)

procs()

tic()
a = @spawnat 2 count_heads(100000000)
b = @spawnat 3 count_heads(100000000)
c = @spawnat 4 count_heads(100000000)
d = @spawnat 5 count_heads(100000000)
out = (fetch(a)+fetch(b)+fetch(c)+fetch(d))/400000000
println("Result = ",out)
toc();

tic()
out2 = count_heads(400000000)/400000000
println("Result = ",out2)
toc();

tic()
nheads = @parallel (+) for i = 1:400000000
	int(randbool())
end
println("Result = ",nheads/400000000)
toc();

tic()
m = {100000000 for i = 1:4}
sum_ans = pmap(count_heads,m)
out3 = sum(sum_ans)/400000000
println("Result = ",out3)
toc();


q = zeros(10000)
@parallel for i = 1:100
	println(myid())
end

cell(3)

myid()

while true
	println("asa")
end

#######################################################################

addprocs(2)

np = nprocs()

asa = Array(Any,10000000)
for i = 1:length(asa)
		asa[i] = i
end

beba = distribute(asa)

t1 = time()
for i = 1:length(asa)
		asa[i] += sqrt(asa[i])*cbrt(asa[i])
end 
t2 = time()
println(string("asa[3] = ", asa[3], " in ", round(t2-t1,2)))

t3 = time()
@sync begin
	for p = 1:np
		if p != myid() || np == 1
			@async begin
				@spawnat p for j in 1:length(localpart(beba)) localpart(beba)[j] += sqrt(localpart(beba)[j])*cbrt(localpart(beba)[j]) end
			end
		end
	end
end
ans = remotecall_fetch(1,getindex,beba,3)
t4 = time()
println(string("beba[3] = ", ans, " in ", round(t4-t3,2)))



############################################################################

addprocs(2)

np = nprocs()

asa = Array(Any,10)
for i = 1:length(asa)
		asa[i] = i
end
masa = Array(Any,10)
for i = 1:length(asa)
		asa[i] = 2
end

fasa = @parallel [2*i for i in [1:10]]
basa = @parallel [3*i for i in [1:10]]

fasa.indexes
basa.indexes

@fetchfrom 3 localindexes(fasa)

chln=[map(length,fasa.indexes[i]) for i = 1:nworkers()]
 
for p = 2:np
	@spawnat p (
		for j in 1:chln[p-1][1]
			localpart(fasa)[j] = (localpart(fasa)[j])-1
			localpart(basa)[j] = (localpart(basa)[j])-2
		end;
		localpart(fasa)[2] = 1789;
		localpart(basa)[2] = 1848
		)
end

fasa
basa 
fetch(@spawnat 2 localpart(fasa))

remotecall_fetch(3, x->(localpart(x)),fasa)

function sq(x::DArray)
	for p = 2:np
		@spawnat p (for j in 1:length(localpart(x)) localpart(x)[j] = (localpart(x)[j])^2 end)
	end 
end

@spawnat 3 (for j in 1:length(localpart(fasa)) localpart(fasa)[j] = (localpart(fasa)[j])^2 end)
sq(fasa)
fasa
procs(fasa)
vec(map(fetch,sq(fasa)))

hcat([1 2], [3 4])

fasa = DArray(I->[asa[i]/masa[i] for i in I[1]],(length(asa),),workers(),[2])
fasa = @parallel [asa[i]/masa[i] for i = 1:10]



hk = remotecall(2,rand,2,2)
workers()
fetch(hk)
qw=DArray(I->[i^2 for i in I[1]],(10,),workers(),[2])
fetch(@spawnat 2 localpart(qw))
qw
remotecall_fetch(2,getindex,qw,3)
remotecall_fetch(2,rand,2,2)

fetch(@spawnat 2 localindexes(beba))
fetch(@spawnat 3 localindexes(beba))
fetch(@spawnat 2 localpart(beba))
fetch(@spawnat 3 localpart(beba))

I = (1:5,)
map(length,I)
I[1]

typeof(I)
map(length,I)
tic()
pbeba(np)
toc()
I = (1:3,1:4)
map(length,I)



beba2[3]
tic()
@spawnat 2 for j in 1:length(localpart(beba)) localpart(beba)[j] += sqrt(localpart(beba)[j])*cbrt(localpart(beba)[j]) end
println(beba[3])
toc()

tic()
@parallel for j = 1:(length(beba)/nworkers())
	localpart(beba)[j] += sqrt(localpart(beba)[j])*cbrt(localpart(beba)[j])
end
toc()
owner(beba,1)
beba

tic()
@spawnat 2 for j in 1:length(localpart(beba)) localpart(beba)[j] += sqrt(localpart(beba)[j])*cbrt(localpart(beba)[j]) end
fetch(@spawnat 3 for j in 1:length(localpart(beba)) localpart(beba)[j] += sqrt(localpart(beba)[j])*cbrt(localpart(beba)[j]) end)
toc()

beba
asa
workers()
ubauba = DArray(I->[i^2+j^2 for i=I[1],j=I[2]],(10,10),workers()[1:4],[2,2])
workers()

fetch(@spawnat 2 localpart(ubauba))
fetch(@spawnat 3 localpart(ubauba))
qasa = DArray(Int64,10,10)

map(length,(1:10,1:10))

der = (1,2)

der[2]

typeof(der)

(1:10,1:10)[1]

println([i for i=1:10])

for i=1:10
	println(i)
end

rmprocs(2,3,4,5)
procs(4)

workers()

S = SharedArray(Int, (3,3), init = S -> S[localindexes(S)] = myid())

asa = [1 2 3; 4 5 6]
length(asa)
nworkers()
addprocs(2)
da = @parallel [2*i for i = 1:10]
procs(da)
da.chunks
fetch(da.chunks[1])
da.indexes
da[3:5]

a = rand(10)

function maxnum_serial(a,s,e)
	if s == e
		a[s]
	else 
		mid = ifloor((s+e)/2)
		low = maxnum_serial(a,s,mid)
		high = maxnum_serial(a,mid+1,e)
		low > high ? low:high 
	end
end

maxnum_serial(a,1,10)

addprocs(2)

v = @parallel [2*i for i = 1:10]
gg = dzeros((10,),workers(),[nworkers()])

for p = procs(v)
	@spawnat p ( localpart(gg)[j]=localpart(v)[j] for j = 1:length(localpart(v)) )
end
gg
for p = procs(v)
	@spawnat p (
		for j = 1:length(localpart(v)) 
			localpart(gg)[j]=localpart(v)[j] 
		end
	)
end
println("Step 1 ok")
for p = procs(v)
	@spawnat p (
		for j = 1:length(localpart(gg)) 
			localpart(gg)[j] += 99 
		end
	)
end
println("Step 2 ok")
ans = @spawnat 2 localpart(gg)[1] -= 104.9  
println("Step 3 ok /n")
gg

gg
v

v.indexes

fetch(@spawnat 2 localpart(v))
fetch(@spawnat 3 localpart(v))

fetch(v.chunks[1])

require("C:/Users/ulzegasi/Julia_files/funtest.jl")

for p = 2:3
	@spawnat p test(localpart(v))
end 

vadd = [i for i = 1:10]

v

for p = 2:3
	@spawnat p (for j = 1:length(localpart(v)) localpart(v)[j] += vadd[localindexes(v)[1][j]] end)
end

v

workers()

@everywhere 
@everywhere 

@everywhere (
	stest = Array(Int64,10);
	rtest = Array(Int64,10);
	stest = [i^2 for i = 1:10];
	rtest = [i+1 for i = 1:10];
)


@everywhere (
	kk = 4;
	qq = 2
)
ww = 17

stest
rtest
qtest = DArray((10,),workers(),[2]) do I 
	[stest[i]/rtest[i] for i in I[1]]
end

qtest2 = DArray((10,),workers(),[2]) do I 
	[stest[i]/rtest[i]+kk-qq+ww for i in I[1]]
end

qtest2

@spawnat 2 (for j = 1:length(localpart(qtest2)) localpart(qtest2)[j] += ww end)
@spawnat 3 (for j = 1:length(localpart(qtest2)) localpart(qtest2)[j] += ww end)

qtest2

mm = DArray(I->[stest[i]/rtest[i] for i in I[1]],(10,))
procs()
addprocs(2)
for p = procs()
	@spawnat p srand(p*18072011)
end

s1 = rand(2)
s2 = fetch(@spawnat 2 rand(4))
s3 = fetch(@spawnat 3 rand(2))
(s1, s2, s3)

s1 = rand(1)
s2 = fetch(@spawnat 2 rand(1))
s3 = fetch(@spawnat 3 rand(1))
(s1, s2, s3)


@spawnat 2 (
	rand()
)

da = @parallel [2*i for i = 1:11]

for wp = workers()               
    @spawnat wp (for j = 1:size(localpart(da),1) localpart(da)[j] = j end)
end
da
fetch(@spawnat 2 length(localpart(da)))
sum({fetch(@spawnat p sum(localpart(da).^2)) for p=procs(da)})
convert(Array, da)
da
beba = [1, 2, 3]
beba.^2

th = [33.0,44.0]

ths = Array(Float64,2)

ths[1:2] = th[1:2]

th = [12.1,66.7]

th = ths

th[1:2] = 88.9

th

da = DArray(I->[(i,j) for i in I[1],j in I[2]],(4,4),workers(),[nworkers(),1])
size(da)
da.indexes
workers()
fetch(@spawnat 2 localpart(da))
fetch(@spawnat 3 localpart(da))
abba = da

da = DArray((4,4),workers(),[nworkers(),1]) do I
	[(i,j) for i in I[1],j in I[2]]
end
da = 0.0
da

for wp = workers()               
    @spawnat wp (for j = 1:size(localpart(da),2) localpart(da)[1,j] = (11,j^2) end)
end

fetch(@spawnat 2 localpart(da))

da
I = (5:8,5:8)
length(I[1])
a = [2*i for i = 1:10]

da = distribute(a)

b = [i^2 for i = 1:10]

da = distribute(b)

da


addprocs(2)
procs()
da = @parallel [2*i for i = 1:11]
dachunks=map(fetch, da.chunks)
size(dachunks,1)==nworkers()

fetch(@spawnat 2 localindexes(da)[1])
fetch(@spawnat 3 localindexes(da)[1][2])

da
rembra = [i for i = 1:11]
workers()
for wp in workers()
	@spawnat wp (
		indvec = [localindexes(da)[1]];
		for j = 1:size(localpart(da),1)   
			localpart(da)[j] = (indvec[j])^2+rembra[indvec[j]]
		end;
	)
end
da
abba =(1:4,)
indvec=[abba[1]]

typeof(indvec)

sum(da)

da[1:3]

a = [2*i for i = 1:50000]
@time(sum(a))

da = @parallel [2*i for i = 1:50000]
@time(sum(da))

const s = 44

function asa()
	println(s)
	for s = 1:10
		println(s)
	end
end
asa()
s


addprocs(2)
workers()
da = @parallel [i for i = 1:21]
dachunks=map(fetch, da.chunks)


critind = [(i-1)*4+1 for i = 1:6]

for ind = 1:21
	if !in(ind,critind)
		println(ind)
	end
end


trop = 
sum(map(fetch, { @spawnat wp (indvec = [localindexes(da)[1]]; loctmp = 0;
	for j = 1:size(localpart(da),1)   
		if !in(indvec[j],critind)
			loctmp += localpart(da)[j]
		end
	end; 
	loctmp ) for wp in workers() } ))
trop

da
size(da,1)
for j = 1:size(da,1)
	da[j] += 1
end

da[20] += 1

tmp =
        sum( 
            map(fetch, 
                { @spawnat wp
                ( indvec = [localindexes(da)[1]]; loctmp = 0;
                    for locj = 1:size(localpart(da),1)
                        if !in(indvec[locj],critind)
                            loctmp += localpart(da)[locj]
                        end
                    end;
                    loctmp ) for wp in workers() } 
                )
            )


	tuz= 
		sum(map(fetch, 
			{@spawnat wp ( indvec = [localindexes(da)[1]]; loctmp = 0; 
				for locj = 1:size(localpart(da),1)
					if !in(indvec[locj],critind) 
						loctmp += localpart(da)[locj]
					end
				end;
				loctmp ) for wp in workers()}))
tuz

trub = 44
@spawnat 1 (trub += sum(da))
trub
da
critind
tonton=findfirst(x->x>7,critind)
7-critind[tonton-1]+1

addprocs(2)
da = @parallel [i for i = 1:21]
db = @parallel [i+100 for i = 1:21]
a = [44 ,55]

ada = (a,da)

ada[1]
ada[2]

procs(ada)

nprocs()
addprocs(2)
workers()
da = @parallel [i for i = 1:20000001];

t1=time()
for wp in workers()
	@spawnat wp (
		indvec = [localindexes(da)[1]];
		for j = 1:size(localpart(da),1)   
			localpart(da)[j] = sqrt(localpart(da)[j]*localpart(da)[j])-1
		end;
	)
end
fetch(@spawnat 2 (localpart(da)[1],localpart(da)[2]))
#asa(da,100)
println("Time is ",time()-t1)

function asa(da::DArray, n::Int64)
	println(da[n])
end
asa(da,1)
asa(da,2)
rmprocs(2,3)
reload("$dir/ParInf_Fun_AD_3.jl")
 
counter = 1

da[1]

datemp = convert(Array,da)

for i = 1:size(datemp,1)
	datemp[i] = db[i]+100
end


da = distribute(datemp)
sum(time_respa_f)
sum(time_respa_f1)
sum(time_respa_f2)
sum(time_respa_f3)
sum(time_respa_f4)

sum(time_respa_f)
sum(time_respa_f1)+sum(time_respa_f3)

mod(12,3)

workers()

addprocs(2)

workers()

da = DArray((16,8),workers(),[nworkers(),1]) do I 
	[i+j for i in I[1], j in I[2]]
end
fetch(da.chunks[1])
fetch(da.chunks[4])
size(da,1)
asaI = (9:12,1:4)
length(asaI[1])
procs(da)

[(i,j) for i in asaI[1], j in asaI[2]]

first(asaI[1])

mod(-1,8)
procs()
rembra = (@spawnat 7 sum(da[:,1]))

fetch(rembra)

rmprocs(workers(procs()))
procs()
addprocs(2)
der_out = dzeros((5,), workers(), [nworkers()])
map(fetch,der_out.chunks)
function asas()
	for wp in workers()
		@spawnat wp (
			indexes = [localindexes(der_out)[1]];
			for j = 1:size(localpart(der_out),1)   
				localpart(der_out)[j] = localpart(der_out)[j]+indexes[j]
			end;
		)
	end
end

asas()

procs(der_out)

ty[2:end-1]

der_out

findin(der_out,30.0)

findin(ty,9)
ty[findin(ty,9)+1]

der_out
procs()
fetch(@spawnat 9 localpart(der_out)[2]+der_out[5])

rmprocs(workers())
procs()

reload("$dir/ParInf_Fun_AD_3.jl")


workers()
wp = 2

@spawnat wp ( uind = [localindexes(u)[1]];
    for jwp = 1:size(localpart(der_out),1) 
        if !in(uind[jwp],ty)
            k = uind[jwp]-ty[findfirst(x->x>uind[jwp],ty)-1]+1;
            localpart(der_out)[jwp] = (K*k/(gamma*dt*(k-1))) * localpart(u)[jwp]
        elseif in(uind[jwp],ty[2:end-1])
            localpart(der_out)[jwp] = ( K/(gamma*j*dt) ) *
            ( 2*localpart(u)[jwp]-u[ty[findin(ty,uind[jwp])-1]] - u[ty[findin(ty,uind[jwp])+1]] )-
                    r[uind[jwp]] * exp(localpart(u)[jwp])*( y[((uind[jwp]-1)/j)+1] -
                        r[uind[jwp]]*exp(localpart(u)[jwp]) )/sigma^2
        elseif uind[jwp] == 1
            localpart(der_out)[jwp] = K*( localpart(u)[jwp] - u[ty[2]] )/( gamma*j*dt )-
                r[1] * exp(localpart(u)[jwp]) * ( y[1] - r[1]*exp(localpart(u)[jwp]) )/sigma^2
        elseif uind[jwp] == N 
            localpart(der_out)[jwp] = K*( localpart(u)[jwp] - u[ty[end-1]] )/( gamma*j*dt )-
                r[N] * exp(localpart(u)[jwp]) * ( y[n+1] - r[N]*exp(localpart(u)[jwp]) )/sigma^2
        end
    end
)

@spawnat wp ( uind = [localindexes(u)[1]];
		jwp = 151;
        if !in(uind[jwp],ty)
            k = uind[jwp]-ty[findfirst(x->x>uind[jwp],ty)-1]+1;
            localpart(der_out)[jwp] = (K*k/(gamma*dt*(k-1))) * localpart(u)[jwp]
        elseif in(uind[jwp],ty[2:end-1])
            localpart(der_out)[jwp] = ( K/(gamma*j*dt) ) *
            ( 2*localpart(u)[jwp]-u[ty[findin(ty,uind[jwp])-1]] - u[ty[findin(ty,uind[jwp])+1]] )-
                    r[uind[jwp]] * exp(localpart(u)[jwp])*( y[((uind[jwp]-1)/j)+1] -
                        r[uind[jwp]]*exp(localpart(u)[jwp]) )/sigma^2
        elseif uind[jwp] == 1
            localpart(der_out)[jwp] = K*( localpart(u)[jwp] - u[ty[2]] )/( gamma*j*dt )-
                r[1] * exp(localpart(u)[jwp]) * ( y[1] - r[1]*exp(localpart(u)[jwp]) )/sigma^2
        elseif uind[jwp] == N 
            localpart(der_out)[jwp] = K*( localpart(u)[jwp] - u[ty[end-1]] )/( gamma*j*dt )-
                r[N] * exp(localpart(u)[jwp]) * ( y[n+1] - r[N]*exp(localpart(u)[jwp]) )/sigma^2
        end
)
fetch(@spawnat wp (uind = [localindexes(u)[1]]; jwp = 151; 
	localpart(der_out)[jwp] = ( K/(gamma*j*dt) ) *
            ( 2*localpart(u)[jwp]-u[ty[findin(ty,uind[jwp])-1]] - u[ty[findin(ty,uind[jwp])+1]] )-
                    r[uind[jwp]] * exp(localpart(u)[jwp])*( y[((uind[jwp]-1)/j)+1] -
                        r[uind[jwp]]*exp(localpart(u)[jwp]) )/sigma^2)
)

ty
fetch(@spawnat wp (uind = [localindexes(u)[1]]; jwp = 151; 
	localpart(der_out)[jwp] = ( K/(gamma*j*dt) ) *
            ( 2*localpart(u)[jwp]-u[ty[findin(ty,uind[jwp])-1][1]] - u[ty[findin(ty,uind[jwp])+1][1]] ))
)

fetch(@spawnat wp (uind = [localindexes(u)[1]]; uind[151]))

fetch(@spawnat wp (indasa = ty[findin(ty,151)-1][1];
	u[indasa]))
)
indasa=ty[findin(ty,151)-1]
fetch(@spawnat wp u[121])

x = zeros(Int64, 10)
ind = 1
function asa(x,ind)
	x[ind] = ind
	ind += 1
end

for i = 1 :10 
	asa(x)
end
x

V_fast_der_theta(theta) 

V_fast_der_u(u)
theta
u
V_fast(theta,u)

V_fast_theta(theta)


function asafun1(x)
	[3*x[1], 2*x[2]^2]
end

asafun1([4,5])
function asafun2(x)
	sum([3*x[1], 2*x[2]^2])
end

using ReverseDiffSource
pot = quote
	H = 0.0
	rembra = y[1]
	tonton = y[2]
	nenno = rembra^2
	bebo = tonton*2
	for i = 1:2
		H += i*x[i]^i
	end
	H += nenno*x[1]+bebo*x[2] 
end
typeof(pot)
y = [1,1]
x = [5,5]

asadiff = rdiff(pot, outsym = :H, y = ones(Float64,2) ) 
@eval foo(y) = ($asadiff)[2][1]
foo([5,6])

eval(pot)
theta
V_fast_der_theta(theta)
a=[1,2]
b=[3,4]
[a,b]

z = 5

Z = 7

z + Z

x = [1,2,3,4,5]
y = [6,7,8,9,10]
z = Array(Int64,5)
beta = Array(Int64,5)
for i = 1:5
	z[i] = x[i] + y[i]
	beta[i] = z[i] - 2
end

beta
z

out_slow = ( exp(- u[N]) +  u[N]*(K * lnr_der[N-1] + gamma + 1) - exp(- u[1]) - u[1]*(K * lnr_der[1] + gamma + 1) )/gamma
    # Potential, bulk terms:
    # ---------------------------
    
    # for i=2:(N-1)
    #     out_slow +=  dt * ( (K * lnr_der[i] + gamma + 1) - exp( -q[i] ) )^2 / (2.0*gamma*K) -
    #             dt * q[i] * ( K*(lnr_der[i]-lnr_der[i-1])/dt ) / gamma
    # end
    for i=1:(n-1)
        out_slow +=  dt * ( (K * lnr_der[i*j+1] + gamma + 1) - exp( -u[i*j+1] ) )^2 / (2.0*gamma*K) -
                        dt * u[i*j+1] * ( K*(lnr_der[i*j+1]-lnr_der[i*j])/dt ) / gamma
    end

    for i = 0:(n-1)
        for k = 2:j
            out_slow += dt * ( (K * lnr_der[i*j+k] + gamma + 1) - exp( -u[i*j+k] ) )^2 / (2.0*gamma*K) -
                            dt * u[i*j+k] * ( K*(lnr_der[i*j+k]-lnr_der[i*j+k-1])/dt ) / gamma
        end
    end

bim = 0.0
for i = 0:(n-1)
        for k = 2:j
            bim += dt * ( (K * lnr_der[i*j+k] + gamma + 1) - exp( -q[i*j+k] ) )^2 / (2.0*gamma*K) -
                            dt * q[i*j+k] * ( K*(lnr_der[i*j+k]-lnr_der[i*j+k-1])/dt ) / gamma
        end
    end

    bim

for i = 1:2
	println(typeof(i))
end
using ReverseDiffSource

asa = Array(Any,10)
bebo = Array(Float64,10)
for i = 1:10
	bebo[i] = 2
end
for i = 1:10
	asa[i] = :(0)
end
for i = 3:3:9
	tmp1 = :(bebo[$i]^2)
	tmp2 = :(bebo[$i])
	asa[i] = :($tmp1+$tmp2)
end
for i in [3]
	tmp1 = :(bebo[$i]^2)
	tmp2 = quote
		tmp0 = 0
		for jj = ($i-1):($i+1)
			tmp0 += bebo[jj]
		end
		tmp0
	end
	asa[i] = :($tmp1+$tmp2)
end
# for i = 3:3:9
# 	asa[i] += bebo[i]
# end
asa
eval(asa[3])
eval(asa[6])
typeof(asa[3])
ii=3
jj=6
zz=9
ex = Expr(:call, :+, Expr(:call, :+, asa[ii], asa[jj]), asa[zz])
for ijz in [3,6,9]
	ex = Expr[]
eval(ex)
tz = 3
qw = 4
testexpr = quote

	# for i = 3:3:9
	# 	asa[i] = bebo[i]^2
	# end
	# for i = 3:3:9
	# 	asa[i] += bebo[i]
	# end
	# rembra = zeros(10)
	# foo = $(asa[3])+$(asa[6])+$(asa[9])
	foo = 3*(bebo[1]^2+bebo[1])+2*(bebo[2]^2+bebo[2])+3*(bebo[4]^2+bebo[4])+2*(bebo[5]^2+bebo[5])+
	 		3*(bebo[7]^2+bebo[7])+2*(bebo[8]^2+bebo[8])
	# for i in [6,9]
	# 	foo += asa[i]
	# end
	tmp1 = 0.0
	for i = 3:3:9
		tmp1 += bebo[i]^2
	end

	tmp2 = 0.0
	for i = 3:3:9
		tmp3 = bebo[i-1]+3*bebo[i+1]
		for jj = (i-1) : (i+1)
			tmp2 += tz*exp(-bebo[jj]) - qw*(bebo[i])
		end
		tmp2 -= tmp3
		tmp2 += bebo[i]^2
	end
	
	foo += (tmp1+tmp2)
end
eval(testexpr) 
testexpr
beboder=rdiff(testexpr, outsym=:foo, bebo = ones(Float64, 10) ) 
@eval res(bebo) = ($beboder)[2][1]
res(bebo)

 asa = :( :a in $( :(:a + :b) ) )
 typeof(asa)
 typeof((asa.args)[1])

 bebo = begin
 	out = 0.0
 	for i = 1:3
 		out += i
 	end
 	for j = 1:3
 		out += j^2
 	end
 	out -= 21
 end

out 

asa= 5.0

fr = asa/.[2.0,3.0]

addprocs(2)

da_sample = zeros(10,4)
a = [i^2 for i in [1:4]]
da = DArray(I->[i^3 for i in I[1]],(4,),workers(),[nworkers()])
da_sample[1,:] = da[1:end]
size(da)
da.indexes
workers()
fetch(@spawnat 2 localpart(da))
fetch(@spawnat 3 localpart(da))
abba = da

da

function tonton(vec::DArray)
	out = @parallel (+) for i = 1:size(vec,1)
		sum(vec[i])
	end
	return out
end

tonton(da)

da_save = dzeros(4)
a_save = zeros(4)
da_save = da
da_save
a_save = da 

t0 = time()
da = DArray(I->[i for i in I[1]],(20000000,),workers(),[nworkers()])
da_save_1 = dzeros(20000000)
da_save_2 = dzeros(20000000)
da_save_3 = dzeros(20000000)
round(time()-t0,6)
t1 = time()
da_save = da 
round(time()-t1,6)

t2 = time()
@async begin
	for wp in workers()
		@spawnat wp for j = 1:size(localpart(da),1)
			localpart(da_save)[j] = localpart(da)[j]
		end
	end
end
da_save
round(time()-t2,6)
da.indexes == da_save.indexes


asa = [1,2,3]
asa0 = zeros(3)
asa0 = asa
asa0
rmprocs(2,3,4,5)
da
function dsave(da::DArray)
	DArray(size(da),procs(da)) do I 
		[da[i] for i in I[1]]
	end
end
da_save = dsave(da)
da_save_2 = da 

for pid in workers()
    remotecall(pid, x->(global da; da=x; nothing), da)
end

t3 = time()
da_save_3 = DArray(size(da),procs(da)) do I 
	out = zeros(length(I[1]))
	for i = I[1]
   	    out[i-minimum(I[1])+1] = da[i]
   	end
   	out
end
round(time()-t3,6)
da.indexes

tref=time()
da_save_1 = da
round(time()-tref,6)

da_save_1

tref=time()
for np in workers()
	@spawnat np (
		for i = 1:size(localpart(da),1)
			localpart(da_save_2)[i] = localpart(da)[i]
		end
	)
end
round(time()-tref,6)
da_save_2

tref=time()
for np in workers()
	@spawnat np (localpart(da_save_3) = localpart(da))
end
round(time()-tref,6)

fetch(@spawnat 2 da)


a = [randn() for i = 1:50000000]
@time(sum(a.^4-a.^3-a.^2+sqrt(abs(a))))

addprocs(4)
da = distribute(a)

@everywhere function strange_sum(numIter)
	temp_vec = [randn() for i = 1:numIter]
	return sum(temp_vec.^4-temp_vec.^3-temp_vec.^2+sqrt(abs(temp_vec)))
end
@time(strange_sum(50000000))

M={convert(Int64,50000000/4) for i = 1:4}
@time(sum(pmap(strange_sum,M)))

convert(Int64,50000000/4)
@time(sum(map(fetch,{(@spawnat p strange_sum(12500000)) for p in workers()})))

rmprocs(2,3,4,5)

function strange_sum_2(ar::Array)
	sum(ar.^4-ar.^3-ar.^2+sqrt(abs(ar)))
end

@time(sum(map(fetch,{(@spawnat p sum(localpart(da).^4-localpart(da).^3-localpart(da).^2+sqrt(abs(localpart(da))))) for p in workers()})))
b = zeros(size(a,1))
a
b = a
db = dzeros(size(da,1))
@time(for p in workers()
	@spawnat p (for j = 1:size(localpart(da),1)
		localpart(db)[j] = localpart(da)[j]
	end)
end)
db[end]

for p in workers()
	@spawnat p (for j = 1:size(localpart(da),1)
		localpart(da)[j] = localpart(da)[j]^2+1
	end)
end

da
da
aw=[1,2,3]

db = 2.*da.+1

xs = zeros(10,4)
xk = [i for i = 1:4]
procs()
dxs = DArray((4,),workers()[1:2]) do I
	[i for i in I[1]]
end

xs[1,:] = xk

xs[2,:] = dxs[1:4]


addprocs(2)
using ReverseDiffSource

np = nworkers()
first_worker = workers()[1]
funtest = Array(Any, 2)
for funind = 1:np
	funtest[funind] = quote
		out = 0.0
		for i = (($funind-1)*5+1):($funind*5)
			out += (i^2)*x[i]^i
		end
		funout = out
	end
end


@everywhere x = ones(10)
eval(funtest[1])+eval(funtest[2])
funtest_der = [rdiff(funtest[ind], outsym=:funout, x = ones(Float64, 10) ) for ind = 1:2]
ind=1
@eval funtest_der_eval(x) = ($(funtest_der[ind]))[2][1]+($(funtest_der[ind+1]))[2][1]
funtest_der_eval(x)
temp3 = Array(Any,2)
for ii = 1:2
	ind = ii
	temp3[ind] = (@eval ($(funtest_der[ind]))[2][1])
end
temp3
sum(temp3)

@everywhere (for i = 1:size(x,1)
	x[i] = x[i] - 1
end)
x
tonton(x) = sum(map(fetch,{@spawnat p (
	ind = (p - (first_worker-1));
	@eval ($(funtest_der[ind]))[2][1]) for p in workers()} ))

tonton(x)


fetch(@spawnat 3 eval(V_fast[2]))
typeof(proc_indexes[1])
rang = [1:5,6:10]

testexpr = quote
	out = 0.0
	for i = proc_indexes[1]
		out += i
	end
	out
end

eval(testexpr)

rag = 1:5
typeof(rag)
rag[end]

x = 5

for i = 2:3 x -= 1 end
x

for pid in workers()
    remotecall(pid, x->(global u; u = x; nothing), u);
    remotecall(pid, x->(global theta; theta = x; nothing), theta)
end

workers()

bibba = 2

fetch(@spawnat 2 (3+bibba))

bebo = [(@spawnat workers()[i] bibba+i) for i = 1:nworkers()]
fetch(bebo[2])

fetch(@spawnat 2 (fetch(bebo[1])+12))

fetch(bebo[1])
id(bebo[1])

bebo[1]

asa = dzeros(nworkers())
nworkers()
asaqq = DArray((nworkers(),),workers(),[nworkers()]) do I
	


end


asaq.indexes
asaq

for p in workers()
	@spawnat p localpart(asaq)[1] = (localindexes(asaq)[1][1])^2
end

asaq
workers()
fetch(@spawnat 2 localindexes(asaq)[1])

rew = (1:5,)

rew[2]

dasa = quote
	rasa = x[1]^2+2*x[2]^2+3*x[3]^2
end
x =[1,2,3]
eval(dasa)
mm = rdiff(dasa, outsym=:rasa, x = ones(Float64, 3) )
@eval ($mm)[2][1]

mm2 = @eval ($(rdiff(dasa, outsym=:rasa, x = ones(Float64, 3))))[2][1]
bvc = [5, 6, 7, 8]
findin(bvc,7)
qq = [89,78]
qq[findin(workers(),3)]




@everywhere function fib(n)
	if (n < 2) then
		return n 
	else 
		return (fib(n-1)+fib(n-2))
	end
end

@time [fib(i) for i = 1:45]

@everywhere function fib_par(n)
	if (n < 40) then
		return fib(n) 
	else 
		x = @spawn fib_par(n-1)
		y = fib_par(n-2)
		return fetch(x)+y
	end
end

@time [fib_par(i) for i = 1:45];

addprocs(2)

@time [fib_par(i) for i = 1:45];

rem(2,3)

N=16
da = dones(3,N)

for p = procs(da) @spawnat p println(localpart(da)) end 

fetch(@spawnat 2 localindexes(da))

[size((da.indexes)[w][2],1) for w = 1:2]

da.indexes[1][2]

div(7,2)
addprocs(1)
workers()
[distribute(zeros(21)).indexes[i][1] for i = 1:nworkers()]

addprocs(1)




da = @parallel [2*i for i = 1:5]
da.indexes
daa = dzeros(5)

@spawnat p localpart(daa)[i]

@spawnat 2 (a2 = [i for i = 1:5])

du_size
pts_indexes
zut
zut = [size(pts_indexes[i],1) for i = 1:size(pts_indexes,1)]
asa = dzeros((5,),workers()[1:1],[1])
asa = dzeros((5,),workers()[1:1],[1])
asa.indexes
Icustom = [(1:5,),(6:14,),(15:19,),(20:24,)]
bunny = DArray((24,)) do 
	[i for i in I[1]]
	end

bunny = DArray((Icustom) -> [i for i in Icustom[1]],(24,))
bunny.indexes
reshape(bunny,6,4)

iround(linspace(1, 24+1, 4+1))

addprocs(4)
chk_indexes
chunk_ind = [(1:5,),(6:14,),(15:19,),(20:24,)]
chunks = Array(RemoteRef,4)
for i = 1:nworkers()
	chunks[i]=remotecall(workers()[i], I->[i for i = I[1]], chunk_ind[i])
end

map(fetch,chunks)

@spawnat 2 (fetch(chunks[1]))[1] = 23 

function myDArray(init, dims, procs, dist)
    np    	= prod(dist)
    procs 	= procs[1:np]
    cuts    = [[chk_indexes[i][1][1] for i = 1:num_workers] for ii = 1:1]
    idxs    = chk_indexes
    chunks  = Array(RemoteRef, dist...)
    for i = 1:np
        chunks[i] = remotecall(procs[i], init, idxs[i])
    end
    p = 0
    for i = 1:length(procs)
        if procs[i] == myid()
            p = i
        end
    end
    p = max(1 , p)
    A = remotecall_fetch(procs[p], r->typeof(fetch(r)), chunks[p])
    DArray{eltype(A),length(dims),A}(dims, chunks, procs, idxs, cuts)
end

p = 0.5

function defaultdist(sz::Int, nc::Int)
    if sz >= nc
        iround(linspace(1, sz+1, nc+1))
    else
        [[1:(sz+1)], zeros(Int, nc-sz)]
    end
end

cuts = map(defaultdist, 24, 4)
typeof(cuts)
cuts[1]
# compute indexes array for dividing dims into chunks
function chunk_idxs(dims, chunks)
	cuts = Array(Any,1)
    cuts[1] = [chunk_ind[i][1][1] for i = 1:num_workers]
    n = length(dims)
    idxs = [(1:5,),(6:14,),(15:19,),(20:24,)]
    idxs, cuts
end

function localpartindex(pmap::Vector{Int})
    mi = myid()
    for i = 1:length(pmap)
        if pmap[i] == mi
            return i
        end
    end
    return 0
end
localpartindex(workers()[1:nworkers()])
localpartindex(d::DArray) = localpartindex(d.pmap)
chunk_idxs(24,4)
myasa = myDArray(I->[i for i = I[1]],(24,),workers(),[nworkers()])

myasa.chunks
myasa.indexes
for p in workers()
	@spawnat p println(localpart(myasa))
end


chunk_ind[1][1][1]
myasa

pts_indexes[3][1]

n = 5
j = 4

hasa = vcat([1,2,3],[4,5,6])

da = DArray(I->[i for i in I[1]],(10,),workers(),nworkers())
da.indexes 

a = [i^2+1 for i = 1:10]

for p in workers()
	@spawnat p (for i = 1:size(localpart(da),1)
		localpart(da)[i] = a[da.indexes[p-1][1][i]]
	end)
end

dda = distribute(a)
a[1] = 1455
a 

x = [i for i = 1:10]
dx = distribute(x)
@everywhere f1(x) = x[1]^2+2*x[2]^2+3*x[3]^2
@everywhere f2(x) = 6*x[1]^2+7*x[2]^2
for pid in workers()
    remotecall(pid, y->(global x; x = y; nothing), x)
end
expr = quote
	out = 0.0
	for i = 1:10
		out += i*x[i]^2
	end
end

fetch(@spawnat 2 f1(localpart(dx)))
fetch(@spawnat 3 f2(localpart(dx)))


der_expr = rdiff(expr, outsym=:out, x = ones(Float64, 10))
@eval $(der_expr)[2][1]
dexpr = DArray(I->[:($i) for i in I[1]],(nworkers(),),workers(),[nworkers()])
dexpr_der = DArray(I->[:($i) for i in I[1]],(nworkers(),),workers(),[nworkers()])
@spawnat 2 (localpart(dexpr)[1] = quote
	out = 0.0
	for i = 1:5
		out += i*x[i]^2
	end
end)
@spawnat 3 (localpart(dexpr)[1] = quote
	out = 0.0
	for i = 6:10
		out += i*x[i]^2
	end
end)


@spawnat 2 (
    localpart(dexpr_der)[1] = rdiff(localpart(dexpr)[1], outsym=:out, x = ones(Float64, 5));
)
@spawnat 3 (
    localpart(dexpr_der)[1] = rdiff(localpart(dexpr)[1], outsym=:out, x = ones(Float64, 10));
)

sum(
    map(
        fetch, {@spawnat p (
            @eval ($(localpart(dexpr_der)[1]))[2][1]) for p in 2:3}
        )
    )
addprocs(1)
addprocs(3)
da = DArray(I->[i for i in I[1]],(10,),workers(),[nworkers()])
da.indexes

@time(fetch(@spawnat 3 (
	out=0.0;
	for i = 1:1000000
		out += localpart(da)[1]+localpart(da)[2]-localpart(da)[3]+localpart(da)[4]
	end;
	final_out = out)
))

@time(fetch(@spawnat 3 (
	out=0.0;
	for i = 1:1000000
		out += da[4]+da[5]-da[6]+da[7]
	end;
	final_out = out)
))
fetch(@spawnat 2 (localpart(dy)))
fetch(@spawnat 3 (localpart(dy)))


fetch(@spawnat 2 ( 

    localpart(V_fast)[1] = quote
        out_fast = 1000.0
        out_temp = 0.0
        for i = 1:6
        	out_temp += localpart(dy)[i]
        end
        out_fast += out_temp
    end;
))
@spawnat 3 (
	for i = 1:size(localpart(dy),1)
		localpart(dy)[i] = 0.1
	end)
@spawnat 3 ( 
	# y_temp = Array(Float64,6);
 	ytemp = localpart(dy);
  
    localpart(V_fast)[1] = quote
        out_fast = 2000.0
        out_temp = 0.0
        for i = 1:6
        	out_temp += ($ytemp[i])
        end
        out_fast += out_temp
    end;
)

fetch(@spawnat 3 (eval(localpart(V_fast)[1])))

@spawnat 2 (

        asa3 = Array(Int64,3);
        asa3 = localpart(asa);
        localpart(V_fast)[1] = quote
        	out = sum([$asa3[i] for i = 1:3])
        end
)


asa = DArray(I->[i^2 for i in I[1]],(6,),workers(),[nworkers()])
asa.indexes
fetch(@spawnat 2 eval(localpart(V_fast)[1]))
asa2 =[i^2 for i = 1:6]
for pid in workers()
    remotecall(pid, x->(global asa2; asa2 = x; nothing), asa2)
end

x = [i for i = 1:10]
dx = distribute(x)
@everywhere f1(x) = x[1]^2+2*x[2]^2+3*x[3]^2
@everywhere f2(x) = 6*x[1]^2+7*x[2]^2
fetch(@spawnat 2 f1(localpart(dx)))
fetch(@spawnat 3 f2(localpart(dx)))
der = DArray(I->[:($i) for i in I[1]],(nworkers(),),workers(),[nworkers()])
@spawnat 2 (localpart(der)[1] = rdiff(f1, (ones(5),)))
@spawnat 3 (localpart(der)[1] = rdiff(f2, (ones(5),)))

fetch(@spawnat 2 localpart(der)[1](localpart(dx)))[2][1]
fetch(@spawnat 3 localpart(der)[1](localpart(dx)))[2][1]

arr = map(fetch, {@spawnat p localpart(der)[1](localpart(dx))[2][1] for p in 2:3})
vcat(arr[1],arr[2])

############################################################################################
############################################################################################
## rdiff with expressions containg localpart of darrays??
############################################################################################
############################################################################################
addprocs(2)
@everywhere using ReverseDiffSource

V_fast     = DArray(I->[:($i) for i in I[1]],(nworkers(),),workers(),[nworkers()])
V_fast_der = DArray(I->[:($i) for i in I[1]],(nworkers(),),workers(),[nworkers()])

x = [i for i = 1:10]

@spawnat 2 ( 

localpart(V_fast)[1] = quote
	out = 0.0
	out_temp = 0.0
	for i = 1:3
		out_temp += i*x[i]^2
	end
	out = out_temp + x[4]
end
)

@spawnat 3 ( 

localpart(V_fast)[1] = quote
	out = 0.0
	out_temp = 0.0
	for i = 6:8
		out_temp += i*x[i]^2
	end
	out = out_temp + x[9]
end
)

eval(fetch(@spawnat 2 localpart(V_fast)[1]))
eval(fetch(@spawnat 3 localpart(V_fast)[1]))

sum(map(eval,map(fetch,{(@spawnat p (localpart(V_fast)[1])) for p in workers()})))

for i = 1:2
    @spawnat workers()[i] (
        localpart(V_fast_der)[1] = rdiff(localpart(V_fast)[1], outsym=:out, x = ones(Float64, 10));
    )
end

V_fast_der_fun(x) = begin
sum(
    map(
        fetch, {@spawnat p (
            @eval ($(localpart(V_fast_der)[1]))[2][1]) for p in workers()}
        )
    )
end
V_fast_der_fun(x)

# ok this works
# but what about the same with the darray dx??

V_fast_2     = DArray(I->[:($i) for i in I[1]],(nworkers(),),workers(),[nworkers()])
V_fast_der_2 = DArray(I->[:($i) for i in I[1]],(nworkers(),),workers(),[nworkers()])

dx = DArray(I->[i for i in I[1]],(10,),workers(),nworkers())

foo = @spawnat 2 quote
	out = x[1]
	for i = 2:5
		out += x[i]
	end
	out 
end

eval(fetch(foo))


d_foo = @spawnat 2 begin
	dxtemp = zeros()
	quote
		out = localpart(dx)[1]
		for i = 2:5
			out += localpart(dx)[i]
		end
		out 
	end
end

eval(fetch(d_foo))

sum(localpart(dx)))	

@everywhere function math_expr(da::DArray)
	global loc_da = localpart(da)
	expr = quote
		out = 0.0;
		out_temp = 0.0;
		for i = 1:3
			out_temp += i*loc_da[i]^2
		end;
		out = out_temp + loc_da[4]
	end
	return expr
end
dx
eval(fetch(@spawnat 2 (localpart(V_fast_2)[1] = math_expr(dx))))

for i = 1:2
    @spawnat workers()[i] (localpart(V_fast_der_2)[1] = rdiff(math_expr(dx), outsym=:out, dx = ones(Float64, 5));)
end

V_fast_der_2

#### It would work well with functions!

x = [i for i = 1:10]
dx = distribute(x)
@everywhere f1(x) = x[1]^2+2*x[2]^2+3*x[3]^2
@everywhere f2(x) = 6*x[1]^2+7*x[2]^2+8*x[3]^2

fetch(@spawnat 2 f1(localpart(dx)))
fetch(@spawnat 3 f2(localpart(dx)))
der = DArray(I->[:($i) for i in I[1]],(nworkers(),),workers(),[nworkers()])
@spawnat 2 (localpart(der)[1] = rdiff(f1, (ones(5),)))
@spawnat 3 (localpart(der)[1] = rdiff(f2, (ones(5),)))

fetch(@spawnat 2 localpart(der)[1](localpart(dx)))[2][1]
fetch(@spawnat 3 localpart(der)[1](localpart(dx)))[2][1]

arr = map(fetch, {@spawnat p localpart(der)[1](localpart(dx))[2][1] for p in 2:3})
out1 = vcat(arr[1],arr[2])


@everywhere fix(x,i) = 5*(x[i])^2

out = dzeros(2)
for p in workers()
	@spawnat p localpart(out)[1] = sum([fix(localpart(dx),i) for i = 1:3])
end
out

@everywhere derfix = rdiff(fix,(ones(10),1))

derfix(x,6)[2][1]
dx.indexes
arr2 = map(fetch, {@spawnat p sum([derfix(localpart(dx),i)[2][1] for i = 1:5]) for p in 2:3})
vcat(arr2[1],arr2[2])

@everywhere fix2(x,i) = 5*(x[i]+x[i+1])^2

out2 = dzeros(2)
for p in workers()
	@spawnat p localpart(out2)[1] = sum([fix2(localpart(dx),i) for i = 1:3])
end
out2

@everywhere derfix2 = rdiff(fix2,(ones(10),1))

derfix2(x,1)

x2 = [1,2,3,4,5,6,5,6,7,8,9,10]
dx2 = distribute(x2)
dx2.indexes

derfix2(x,1)[2][1]
derfix2(x,2)[2][1]
derfix2(x,3)[2][1]
derfix2(x,4)[2][1]
derfix2(x,5)[2][1]
derfix2(x,6)[2][1]
derfix2(x,7)[2][1]
derfix2(x,8)[2][1]
derfix2(x,9)[2][1]

fetch(@spawnat 3 derfix2(localpart(dx),4)[2][1])


arr2[1] = fetch(@spawnat 2 sum([derfix2(localpart(dx2),i)[2][1] for i = 1:5]))[1:5]
arr2[2] = fetch(@spawnat 3 sum([derfix2(localpart(dx2),i)[2][1] for i = 1:5]))[2:6]
vcat(arr2[1],arr2[2])


############################################################################################
############################################################################################

###########################################################################################
## Could we use expressions from Mathematica????
###########################################################################################
asa = readdlm("C:/Users/ulzegasi/Julia_files/ParInf_HMC/testdiff")
asa[1,2]
l1 = size(asa,1)
l2 = size(asa,2)
out = Array(String,l1)
for jj = 1:l1
	out[jj] = ""
end
for jj = 1:l1 
	for ii = 1:l2
		out[jj] = "$(out[jj])$(asa[jj,ii])"
	end
end
out
for jj = 1:l1 
	out[jj] = replace(out[jj], "E^", "exp")
end
out

################### EXAMPLE #######################
const j = 4
const n = 5
dir   = "C:/Users/ulzegasi/Julia_files/ParInf_HMC"       # Main project directory, local repository   
dir2  = "C:/Users/ulzegasi/SWITCHdrive/JuliaTemp/Data"   # Secondary directory (e.g., data storage)

trange = 2:5002
t  = float64(readdlm("$dir/t.dat")[trange,2]) # Time points t.dat

const N = int64(n*j+1)   # Total number of discrete time points
                         # IMPORTANT: (N-1)/j = integer = n (measurement points)
if  ((N-1)%j) != 0 
    error("Be careful, the number of staging points j must fulfill (N-1)/j = integer !")
end

t = linspace(t[1], t[end], N);                # Time points      
T  = t[end]-t[1]                              # Total time interval

dt = T/(N-1)                                  # Time step
ty  = iround(linspace(1, N, n+1))             # Indeces of "end point" beads (= "measurement" points)

const s = 2                                   # Number of system parameters (k, gamma)     

true_K     = 50.                              # Retention time
true_gamma = 0.2                              # Dimensionless noise parameter
sigma      = 0.05                             # Measurement noise
true_theta = [log(true_K/T),log(true_gamma)]  # Parameters to be inferred (beta, tau)

K     = 200.0                                 # Initial state
gamma = 0.5                                   
theta = [log(K/T),log(gamma)]

K = T*exp(theta[1])
g = exp(theta[2])

r = sin(t/100.).*sin(t/100.) + 0.1;      # Simple sinusoidal model
lnrder      = zeros(N)             # Logarithmic derivative of rain input
for i = 1:(N-1)                          # Calculate logarithmic derivative of rain input                
    lnrder[i] = ( log(r[i+1]) - log(r[i]) ) / dt
end

S            = zeros(N)             # System realization (with true parameters, used to create measurements y)
S_init       = zeros(N)
srand(18072011)                 # Seeding the RNG provides reproducibility
#srand(time())                  # Seeding the RNG with actual time provides randomness

S[1] = true_K * r[1]            # Unperturbed steady state (with constant rain input)
for i = 2 : N
    S[i] = S[i-1] + dt * ( r[i-1] - S[i-1]/true_K ) +
                   sqrt(dt) * sqrt(true_gamma/true_K) * S[i-1] * randn()
end
y  = max(0.01,(1/true_K)*S[ty] + sigma*randn(n+1))      # Generation of measurement points

for a = 1:n
    S_init[ty[a]:ty[a+1]] = K*linspace(y[a],y[a+1],ty[a+1]-ty[a]+1)
end

q            = zeros(N)             # Coordinates vector
u            = zeros(N)             # Coordinates vector
q = log(S_init./(K*r))

## Transformations q -> u 
## (eqs 2.16-17 in Tuckerman et al., JCP 99 (4), 2796, 1993)

for i=0:(n-1)
    u[i*j+1] = q[i*j+1]
    for k=2:j
        u[i*j+k] = q[i*j+k] - ( (k-1)*q[i*j+k+1] + q[i*j+1] )/k
    end
end;
u[N] = q[N]
u
outparsed = Array(Expr,l1)
for jj = 1:l1
	outparsed[jj] = parse(out[jj])
end
outparsed

sum([eval(outparsed[i]) for i = 1:l1])


### AND WHAT IF WE PARALLELIZED?!?! ###############
addprocs(2)

@everywhere outparallel = Array(String,l1)
for jj = 1:l1
	outparallel[jj] = out[jj]
end
for jj = 1:l1 
	outparallel[jj] = replace(outparallel[jj], "u", "localpart(du)")
end
outparallel

outpp = Array(Expr,l1)
for jj = 1:l1
	outpp[jj] = parse(outparallel[jj])
end
outpp

du = distribute(u)
du.indexes

for pid in workers()
    remotecall(pid, x->(global K; K = x; nothing), K)
    remotecall(pid, x->(global g; g = x; nothing), g)
    remotecall(pid, x->(global lnrder; lnrder = x; nothing), lnrder)
    remotecall(pid, x->(global dt; dt = x; nothing), dt)
    remotecall(pid, x->(global du; du = x; nothing), du)
end
eval(outparsed[1])
fetch(@spawnat 2 eval(outpp[1]))

eval(outparsed[2])
fetch(@spawnat 2 eval(outpp[2]))

eval(outparsed[3])
fetch(@spawnat 2 eval(outpp[3]))

eval(outparsed[4])
fetch(@spawnat 2 eval(outpp[4]))

eval(outparsed[5])
fetch(@spawnat 2 eval(outpp[5]))

eval(outparsed[6])
fetch(@spawnat 2 eval(outpp[6]))

eval(outparsed[7])
fetch(@spawnat 2 eval(outpp[7]))

eval(outparsed[8])
fetch(@spawnat 2 eval(outpp[8]))

fetch(@spawnat 2 localpart(du))
fetch(@spawnat 3 localpart(du))

outpp[9]

asa="E^2"
asa=replace(asa,"E^","exp(")

eval(parse(asa))

##############################################################################
##############################################################################

V_fast_par     = DArray(I->[:($i) for i in I[1]],(nworkers(),),workers(),[nworkers()])
V_fast_der_par = DArray(I->[:($i) for i in I[1]],(nworkers(),),workers(),[nworkers()])
V_fast_par.indexes
@spawnat 2 (
localpart(V_fast_par)[1] = quote       # Fast dynamics (harmonic and data terms)
	
    	# Variables:
    	# ---------------------------
    	
    	bet      = theta[1]
    	tau      = theta[2]
    	
    	K        = T*exp(bet)
    	gamm     = exp(tau)
	
    	# Harmonic term:
    	# ---------------------------
    	out_fast = 0.0
    	out_temp = ( y[1] - r[1]*exp(u[1]) )^2/(2*sigma^2)
    	for i=1:2
    	    out_temp += K * ( u[(i-1)*j+1] - u[i*j+1] )^2/( j*dt )   / ( 2*gamm ) +
    	                ( y[i+1] - r[i*j+1]*exp(u[i*j+1]) )^2/(2*sigma^2)
    	end
    	for i=1:2
    	    for k=2:j
    	        out_temp += k*u[(i-1)*j+k]^2/( (k-1)*dt )*K/(2*gamm)
    	    end
    	end
    
    	out_fast = out_temp
    
end;

)           


@spawnat 2 (
	localpart(V_fast_der_par)[1] = rdiff(localpart(V_fast_par)[1], outsym=:out_fast, u = ones(21))
)
@spawnat 2 (begin
	temp = string(localpart(V_fast_der_par)[1]);
	temp = replace(temp, "u","localpart(du)");
	localpart(V_fast_der_par)[1] = parse(temp);
end)

fetch(@spawnat 2 eval(localpart(V_fast_der_par)[1])[2][1])

temp = string(V_fast_der_u)
temp = replace(temp, "u","localpart(du)")
V_fast_der_du = parse(temp)
V_fast_der_du
eval(V_fast_der_u)
fetch(@spawnat 2 eval(V_fast_der_du))

@eval V_fast_der(du) = ($V_fast_der_u)[2][1]
@time(
for jj = 1:1000
    V_fast_der(theta,u)
end
)
addprocs(3)
asa = [i for i = 1:60000000]
dasa = distribute(asa)

@time(
	begin
		asa[:] = asa[:].*2+asa[:].*3
		for p = 1:nworkers()
			@spawnat workers()[p] (localpart(dasa)[:] = localpart(dasa)[:].*2+localpart(dasa)[:].*3)
		end
		#asa[1]+dasa[end]
	end)
dasa
asa[:]=asa[:].*2+asa[:].*3

asa2 = Array(Any,5)
asa2 = [1,2,3,4,5]
sqrt(asa2)

procs()

@fetchfrom 3 ((localpart(dpu)[:]).^2)./ (2*localpart(dmu)[:])

capa = [i for i = 1:10]
asanov = quote
	out = capa[end]+capa[1]
end
eval(asanov)
asanovic = rdiff(asanov, outsym=:out, capa = ones(Float64, 10))

asato = [1,2,3,4,5]
marat = [1,1,1,2,3]

forbe = [6,7,8,9,10]

mensh = [1,2,1,2,1]

2*(asato+marat)



asato[:] += 2*(forbe[:]+marat[:])./mensh[:]
dpu

asato[2:(end-1)]

@time(begin
	for jj = 1:1000
		@spawnat 2 (begin
    		# force_u = -V_slow_der_u()
        	localpart(dpu)[:] += (dtau/2)*(-V_slow_der_u()[:])   # (dtau/2)*force_u[:] 
    	end)
	end
end)
@spawnat 2 (begin
        force_u = -V_slow_der_u()
        localpart(dpu)[:] += (dtau/2)*force_u[:] 
    end)

@spawnat 2 (@everywhere asani = 58)
@fetchfrom 3 5+asani
asani

dasa = DArray(I->[i for i in I[1]],(10,))
darat = DArray(I->[i^2 for i in I[1]],(10,))



@everywhere function testfun(a::Array)
	b = Array(Float64,size(a,1))
	b[:] = a[:] + 1
	return b
end 
@fetchfrom 2 (testfun(localpart(dasa)))

@fetchfrom 2 localpart(dasa)[:] += 2*testfun(localpart(dasa))[:]

ambasa = [i for i = 1:10]

ambasa[[1:4, 6:10]]

asa = [i for i = 1:10000000]
dasa = distribute(asa)
dasa.indexes

@time(begin
	for jj = 1:length(asa)
		asa[jj] += asa[jj]^2+1
	end
end)

@time(begin
	for p = 1:num_workers
		@spawnat workers()[p] (localpart(dasa)[:] += localpart(dasa)[:].^2+1)
	end
	dasa
end)


dpua = convert(Array,dpu)
forcea = convert(Array,force)
@time(begin
	for i=1:(du_size) 
    	dpua[i] += (dtau/2)*forcea[i]
	end
end)
@time(begin for p = 1:num_workers
    @spawnat workers()[p] (begin
        localpart(dpu)[:] += (dtau/2)*localpart(force)[:]
    end)
end 
end)
@time(begin
	compute_something(dpu, force)
end)

asa.fields
@everywhere dtau = 0.02
function compute_something(da::DArray, db::DArray)
	DArray(size(da), procs(da)) do I
		[da[i] + (dtau/2)*db[i] for i in I[1]]
	end
end
force
dpu_save = dpu

compute_something(dasa,basa)


asa = [i for i = 1:10]
basa = distribute(asa*2)

basa = slice(asa,1:5)

2*basa[1:5]
a = b = zeros(size(dpu))
dasa = distribute(asa)
dasa.indexes  
@fetchfrom 2 typeof(dasa)

##############################################################
##############################################################
#######   SERIAL VS PARALLEL: A MATCH WITHOUT WINNER   #######

## TASK-BASED PROBLEM , very large number of iterations 
#######################################################
addprocs(2)

@everywhere function count_heads(n)
	c::Int = 0
	for i = 1:n
		c += randbool()
	end
	return c
end

niter = 1000000000
piter = int64(niter/nworkers())
@time(count_heads(niter)/niter*100)  # Serial version: 3.6 seconds

@time(begin   # 3.6 seconds 
	sum([(@fetchfrom p count_heads(piter)) for p in workers()])/niter*100
end)		  # The command @fetchfrom makes the loop wait   
			  # for the result of one proc before moving to the next.
			  # --> equivalent to the serial code!

@time(begin 							# Parallel version: 2.0 seconds.
	a = @spawnat 2 count_heads(piter)
	b = @spawnat 3 count_heads(piter)
	(fetch(a)+fetch(b))/niter*100
end)

# or equivalently:

@time(begin
	sum(map(fetch, {(@spawnat p count_heads(piter)) for p in workers()}))/niter*100
end)

## TASK-BASED PROBLEM , small number of iterations 
#######################################################

niter = 1000
piter = int64(niter/nworkers())
@time(count_heads(niter)/niter*100)     # Serial version: 2-3 * 10^(-5) seconds

@time(begin 							# Parallel version: 0.015 seconds. 
	a = @spawnat 2 count_heads(piter)   # Much slower than the serial one!
	b = @spawnat 3 count_heads(piter)
	(fetch(a)+fetch(b))/niter*100
end)

## Something like our DATA-BASED PROBLEM , small number of iterations 
#####################################################################

addprocs(2)

n = 300
p = [2*i for i = 1:n]
force = [i for i = 1:n]

time_ser = 0.0
for iter = 1:200
	t0 = time()
	begin
		for i = 1:n 
    		p[i] = p[i] + (1.0/2.0) * force[i]
		end
		p
	end
	time_ser += (time()-t0)
end
time_ser 			        # Serial version: 0.03 sec for 200 iterations

@time(begin
	for i = 1:n 
    	p[i] = p[i] + (1.0/2.0) * force[i]
	end
end)  			            # 0.0002-3 sec for a single iteration

n = 300
p = [2*i for i = 1:n]
force = [i for i = 1:n]
dp = distribute(p)			# Introduce distributed arrays to parallelize
dforce = distribute(force)

time_par = 0.0
for iter = 1:200
	t0 = time()
	@sync begin
		for pid in workers()
	    	@spawnat pid (localpart(dp)[:] += (1.0/2.0) * localpart(dforce)[:])
		end
		dp
	end
	time_par += (time()-t0)
end
time_par 					# 0.27 sec  for 200 iterations
							# --> much slower than the serial version!

@time(begin
	for pid in workers()
	   	@spawnat pid (localpart(dp)[:] += (1.0/2.0) * localpart(dforce)[:])
	end
end)                        # 0.003 sec for a single iteration

## Something like our DATA-BASED PROBLEM , large number of iterations 
#####################################################################

addprocs(2)

n = 30000 ################# 30k
p = [2*i for i = 1:n]
force = [i for i = 1:n]

time_ser = 0.0
for iter = 1:200
	t0 = time()
	begin
		for i = 1:n 
    		p[i] = p[i] + (1.0/2.0) * force[i]
		end
		p
	end
	time_ser += (time()-t0)
end
time_ser 				# Serial version: 2.9 seconds for 200 iterations

@time(begin
	for i = 1:n 
    	p[i] = p[i] + (1.0/2.0) * force[i]
	end
end)  					# 0.02 - 0.015 seconds for a single iteration

p = [2*i for i = 1:n]
force = [i for i = 1:n]
dp = distribute(p)
dforce = distribute(force)

time_par = 0.0
for iter = 1:200
	t0 = time()
	@sync begin
		for pid in workers()
	    	@spawnat pid (localpart(dp)[:] += (1.0/2.0) * localpart(dforce)[:])
		end
	end
	time_par += (time()-t0)
end
time_par 				# Parallel version: usually 2.7 seconds, but with fluctuations (sometimes more than 3 seconds)
						# Generally, about the same as the serial version 

@time(@sync begin
	for pid in workers()
	   	@spawnat pid (localpart(dp)[:] += (1.0/2.0) * localpart(dforce)[:])
	end
end)                    # 0.02 sec

n = 300000           ################# 300k
p = [2*i for i = 1:n]
force = [i for i = 1:n]

time_ser = 0.0
for iter = 1:100
	t0 = time()
	begin
		for i = 1:n 
    		p[i] = p[i] + (1.0/2.0) * force[i]
		end
		p
	end
	time_ser += (time()-t0)
end
time_ser 				# Serial version: 16 seconds for 100 iterations

@time(begin
	for i = 1:n 
    	p[i] = p[i] + (1.0/2.0) * force[i]
	end
end)  					# 0.15 seconds for a single iteration

p = [2*i for i = 1:n];
force = [i for i = 1:n];
dp = distribute(p);
dforce = distribute(force);

time_par = 0.0
for iter = 1:100
	t0 = time()
	@sync begin
		for pid in workers()
	    	@spawnat pid (localpart(dp)[:] += (1.0/2.0) * localpart(dforce)[:])
		end
	end
	time_par += (time()-t0)
end
time_par 				# Parallel version: about 16 seconds
						# Same as the serial version 



n = 3000000           ################# 3000k
p = [2*i for i = 1:n];
force = [i for i = 1:n];

time_ser = 0.0
for iter = 1:10
	t0 = time()
	begin
		for i = 1:n 
    		p[i] = p[i] + (1.0/2.0) * force[i]
		end
		p
	end
	time_ser += (time()-t0)
end
time_ser 				# Serial version: 29 seconds for 10 iterations

@time(begin
	for i = 1:n 
    	p[i] = p[i] + (1.0/2.0) * force[i]
	end
end)  					# 2.9 seconds for a single iteration

p = [2*i for i = 1:n];
force = [i for i = 1:n];
dp = distribute(p);
dforce = distribute(force);

time_par = 0.0
for iter = 1:10
	t0 = time()
	@sync begin
		for pid in workers()
	    	@spawnat pid (localpart(dp)[:] += (1.0/2.0) * localpart(dforce)[:])
		end
	end
	time_par += (time()-t0)
end
time_par 				# Parallel version: 17 - 20 seconds
						# Finally faster than the serial version 

@time(@sync begin
	for pid in workers()
	   	@spawnat pid (localpart(dp)[:] += (1.0/2.0) * localpart(dforce)[:])
	end
end)                    # 1.5 seconds

asa = [1,2,3,4,-1]

if asa[1:end] < 0
	println("Careful!")
end

asa = [1,2,3]
(2).^asa
exp(asa)
exp(3)
xep = zeros(100000000)
asa = [i for i = 1:100000000]
@time(begin
	for i = 1:length(asa)
		xep[i] = asa[i]
	end
end)
@time(begin
	xep = asa
end)

xep[101]

gamma = 5
gamma(44)

Kuppa = 5

kuppa 

u = zeros(301)
u[1] = 0.354
for i = 2:301
	u[i] = u[i-1] + 0.832/i
end
for i = 1:11
	bq[i] = i 
end
for i = 1:301
	lnr_der[i] = i/10.0 
end
theta
[1,2,3]+[4,5,6]
asa = 4

(asa *= 6) + 8.5

cos(pi/2)

require("$dir/function.jl")

kappa = 78

testa(3)

zebra = 5
zebra_old = zebra

zebra += 15

asas = 4+ zebra_old