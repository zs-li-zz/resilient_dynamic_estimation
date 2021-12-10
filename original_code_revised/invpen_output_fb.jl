#=
with oracle feedback, used to quantify the estimation error


=#

using LinearAlgebra, GaussianDistributions, Random
using Plots; pgfplotsx() #pgfplotsx()  PGFPlotsX;
using JLD

include("origin_utility.jl")

m=0.1;M=1;l=1
g=9.8
bx = 0.05
bt = 0.1


Acon=  [0.    1.    0.    0.
        0.   -0.05 -0.98  0.1
        0.    0.    0.    1.
        0.    0.05 10.78 -1.1 ]
Bcon=[0; 1/M; 0; -1/M/l]
Ts=0.02
A_in=exp(Acon*Ts)

B_in= [ 2.00000000e-02  1.99933382e-04 -1.35137313e-06 1.37052059e-07
        0.00000000e+00  1.99900096e-02 -1.94575387e-04 1.85551482e-05
        0.00000000e+00  7.41919822e-08  2.00143371e-02 1.98612613e-04
        0.00000000e+00  9.92731564e-06  2.14097916e-03 1.97958303e-02]*Bcon

C_in=[1 0 0 0
      1 0 0 0
      1 0 0 0
      0 0 1 0]

n=size(A_in,1)
m=size(C_in,1)
Q_in=Ts^2*Diagonal([0.1; 0.1; 0.01; 0.01])
R_in=Ts^2*Diagonal([0.1; 0.1; 0.01; 0.01])
Σ_in=Q_in

K_lqr= [-30 -15 -115 -32]
function LQGcontrol(X)
    u=-K_lqr*X # u is a scalar
    return u[1]
end
Λ, K_km, C, T = preprocess(A_in, B_in, C_in, Q_in, R_in, Σ_in)
# T is the transformation matrix from non-Diagonal to Diagonal
# K is the Kalman fixed gain
x0=[ 0 ; 1 ; 0 ; 1 ]

## get data
time_scale = 401
# w=load("w.jld","w")
# v=load("v.jld","v")
# X=load("X.jld","X")
# Y=load("Y.jld","Y")
# Ya=load("Ya.jld","Ya")

w=zeros(n,time_scale)
v=zeros(m,time_scale)
X=zeros(n,time_scale)
Y=zeros(m,time_scale)
Ya=zeros(m,time_scale)

# a=4*(rand(time_scale).-0.5)
a=0.02.*(1:time_scale)
#
for k=1:time_scale
    w[:,k]=rand(Gaussian(zeros(n),Q_in))
    v[:,k]=rand(Gaussian(zeros(m),R_in))
    if k==1
        X[:,k]=x0
    else
        X[:,k]=A_in*X[:,k-1]+w[:,k-1]+B_in*LQGcontrol(X[:,k-1])
    end
    Y[:,k]=C_in*X[:,k]+v[:,k]
    Ya[:,k]=Y[:,k]
    Ya[3,k]=Y[3,k]+a[k]
end


# i=3
# time_axis=[0:MAX_TIME-1].*0.1
# state
# plot(time_axis, X[i,1:time_scale], label = "Oracle State", linecolor = "black", line = (:solid, 1))
# plot(time_axis, X[i,1:MAX_TIME], label = "my State", linecolor = "blue", line = (:solid, 1))
# plot!(time_axis, Xls[i,1:time_scale], label = "my Estimation", linecolor = "blue", line = (:dot, 2))
# plot!(time_axis, Xkm[i,1:time_scale], label = "Kalman State", linecolor = "red", line = (:solid, 1))
# plot!(time_axis, Xkm_hat[i,1:time_scale], label = "Kalman Estimation", linecolor = "red", line = (:dot, 2))

## not transformed kalman
Xkm = zeros(n,time_scale) # real state
Xkm_hat = zeros(n,time_scale) # estimation
Ykm = zeros(m,time_scale) # real state
Ykma = zeros(m,time_scale)
for k=1:time_scale
    # real state
    if k==1
        Xkm[:,k]=x0
    else
        Xkm[:,k]=A_in*Xkm[:,k-1]+w[:,k-1]+B_in*LQGcontrol(Xkm_hat[:,k-1])
    end
    Ykm[:,k]=C_in*Xkm[:,k]+v[:,k]
    Ykma[:,k]=deepcopy(Ykm[:,k])
    Ykma[3,k]=Ykma[3,k]+a[k]
    if k==1
        Xkm_hat[:,k]=x0
    else
        Xkm_hat[:,k]=(A_in-K_km*C_in*A_in)*Xkm_hat[:,k-1]+(I-K_km*C_in)*B_in*LQGcontrol(Xkm_hat[:,k-1])+K_km*Ykma[:,k]
    end

end

## lasso
ζ=complex(zeros(mn,time_scale))
Xls = zeros(n,time_scale)
Xl = zeros(n,time_scale) # real state
Yl = zeros(m,time_scale) # real state
Yla = zeros(m,time_scale)

# F=complex(zeros(n,n,m));
# for i=1:m
#     F[:,:,i]=V*Diagonal(V^(-1)*K[:,i])
# end

for k=1:time_scale
    # real state
    if k==1
        Xl[:,k]=x0
    else
        Xl[:,k]=A_in*Xl[:,k-1]+w[:,k-1]+B_in*LQGcontrol(Xls[:,k-1])
    end
    Yl[:,k]=C_in*Xl[:,k]+v[:,k]
    Yla[:,k]=deepcopy(Yl[:,k])
    Yla[3,k]=Yla[3,k]+a[k]


    println("\n=========================================================")
    println("solving LASSO at k = ", k)
    # zeta estimator
    if k==1
        ζ[:,k]=initialize_ζ( inv(T)*x0 )
    # elseif k<=5
    #     ζ[:,k]=update_ζ( ζ[:,k-1], Yla[:,k], LQGcontrol(Xkm_hat[:,k-1]) )
    else
        ζ[:,k]=update_ζ( ζ[:,k-1], Yla[:,k], LQGcontrol(Xls[:,k-1]) )
    end
    # test_zero=zeros(n,1)
    # for i=1:m
    #     test_zero=test_zero+F[:,:,i]*ζ[(i-1)*n+1:i*n, k]
    # end
    # @show norm(T*test_zero-Xkm_hat[:,k], Inf)

    # solve opt problem
    γ = 10
    x, problem_status= solve_opt(ζ[:,k], γ, 0)
    println("problem_status",problem_status)
    # Xls[:,k] = T*x
    if problem_status || k==1
        Xls[:,k] = T*x
    else
        Xls[:,k] = Xkm_hat[:,k]#A_in*Xls[:,k-1]+B_in*LQGcontrol(X[:,k-1]) # soluation not reliable, use prediction
    end
    println("     real x= ", X[:,k])
    println("estimated x= ", Xls[:,k])
    # end

end

MSEkm=norm(X[:,1:time_scale]-Xkm_hat[:,1:time_scale],2)^2/time_scale
MSEls=norm(X[:,1:time_scale]-Xls[:,1:time_scale],2)^2/time_scale


time_axis=[0:time_scale-1].*Ts
i=1
# compare state
# plot(time_axis, X[i,1:time_scale], label = "Oracle State", linecolor = "black", line = (:solid, 1))
# plot!(time_axis, Xkm[i,1:time_scale], label = "Kalman State", linecolor = "blue", line = (:solid, 1))
# plot!(time_axis, Xls[i,1:time_scale], label = "Kalman Estimation", linecolor = "red", line = (:dot, 2))
# compare error
# plot(time_axis, X[:,1:time_scale]'-Xkm_hat[:,1:time_scale]', label = "kalman est error", linecolor = "black", line = (:solid, 1))
# plot!(time_axis, X[:,1:time_scale]'-Xls[:,1:time_scale]', label = "our est error", linecolor = "blue", line = (:dot, 2))

# kalman control result
plot(time_axis, Xkm[i,1:time_scale], label = "Kalman State", linecolor = "black", line = (:solid, 1))
plot!(time_axis, Xkm_hat[i,1:time_scale], label = "Kalman Estimation", linecolor = "blue", line = (:dot, 2))


# our control result
plot(time_axis, Xl[i,1:time_scale], label = "Our State", linecolor = "black", line = (:solid, 1))
plot!(time_axis, Xls[i,1:time_scale], label = "Our Estimation", linecolor = "blue", line = (:dot, 2))
# plot!(time_axis, Xkm_hat[i,1:time_scale], label = "Kalman Estimation", linecolor = "red", line = (:solid, 1))
plot!(xlabel="time / s", legend=:bottomright, size=(600,400))
# plot!(ylabel="cart position")
# plot!(ylabel="cart velocity")
# plot!(ylabel="pendulum angle")
# plot!(ylabel="pendulum angle velocity")

##
# plot(time_axis, Y[3,:],linecolor=:blue, label="")
# plot!(title="original measurement")
# plot!(xlabel="time / s", ylabel="cart position", size=(600,400))
#
# ##
# plot(time_axis, a, linecolor=:red,label="")
# plot!(title="random signal injection attack")
# plot!(xlabel="time / s", ylabel="cart position",  size=(600,400))
#
# ##
# plot!(time_axis, Ya[3,:],linecolor=:red, label="")
# plot!(title="manipulated measurement")
# plot!(xlabel="time / s", ylabel="cart position", size=(600,400))

# using DelimitedFiles
# writedlm("rand2_Xl.txt", Xl')
# writedlm("rand2_Xkm.txt", Xkm')
writedlm("slope1_Xl.txt", Xl')
writedlm("slope1_Xkm.txt", Xkm')
