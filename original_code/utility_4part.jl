# The used kalman parameter K is transformed from original kalman

using LinearAlgebra, GaussianDistributions, Random
using Convex, SCS

function preprocess(A_in, B_in, C_in, Q_in, R_in, Σ_in, MAX_TIME_IN=1000)
    global n=size(A_in,1)
    global m=size(C_in,1)
    global mn=m*n
    global MAX_TIME = MAX_TIME_IN
    global Λ, B, C, T = pre_transform(A_in, B_in, C_in)
    # global Λ = A_in
    # global C = C_in
    # global T = Matrix(1.0I, n, n)
    # global B_in
    global Σ_in
    global Q_in
    # O1A=obsv(Λ, C[1,:]')*Λ
    global R_in
    global K = get_kalman_K(Q_in ,R_in, Σ_in, B_in)
    # @show K
    # global n_u = n
    if !isempty(findall( <(1), broadcast(abs, diag(Λ)) ))
        global n_u = findall( <(1), broadcast(abs, diag(Λ)) )[1]-1
    else
        global n_u = n
    end
    global n_s=n-n_u
    global Π, S = get_Mtildeinv()
    return Λ, T*K, C, T
end




## pre transformation to form Diagonal state matrix
function pre_transform(A, B, C)
    T = [1.00000000000000	-0.0322072577361264	0.0236972934966003 	-104818.391490053
        0.00000000000000	-0.0894451246483580	-0.0919862406134817 	4764.38143136583
        0.00000000000000	0.360078404080100	-0.257617806082263 	-22.0004197816739
        0.00000000000000	1.00000000000000 	1.00000000000000 	1.00000000000000 ]
    # T=Matrix(1.0I,n,n)
    Λ = T^(-1)*A*T
    C = C*T
    B = T^(-1)*B


    return Λ, B, C, T
    # T^(-1)*A*T==Λ
end


## get Kalman K
function get_kalman_K(Q_in,R_in,Σ_in, B_in)
    # simulate system data ####################################

    w=zeros(n,MAX_TIME)
    v=zeros(m,MAX_TIME)

    X=complex(zeros(n,MAX_TIME))
    Y=complex(zeros(m,MAX_TIME))

    for k in 1:MAX_TIME
        w[:,k]=rand(Gaussian(zeros(n),Q_in)) # (第一个w没用到)
        v[:,k]=rand(Gaussian(zeros(m),R_in))
        if k==1
            X[:,k]=rand(Gaussian(zeros(n),Σ_in))  # 初值随机生成
        else
            X[:,k]=A_in*X[:,k-1]+w[:,k-1]
        end
        Y[:,k]=C_in*X[:,k]+v[:,k]
    end

    Xhat=complex(zeros(n,MAX_TIME)) # Xhat(:,k) means \hat{x} (k)
    Xhatpost=complex(zeros(n,MAX_TIME)) # Xhatpost(:,k) means \hat{x} (k|k-1)
    P=complex(zeros(n,n,MAX_TIME)) # P(:,:,k) means P(k)
    Ppost=complex(zeros(n,n,MAX_TIME) )# Ppost(:,k) means \hat{x} (k|k-1)
    Kmat=complex(zeros(n,m,MAX_TIME))

    for k in 1:MAX_TIME
        if k==1
            Ppost[:,:,k]=Σ_in
            Xhatpost[:,k]=zeros(n,1)
        else
            Ppost[:,:,k]=A_in*P[:,:,k-1]*A_in'+Q_in
            Xhatpost[:,k]=A_in*Xhat[:,k-1]
        end
        Kmat[:,:,k]=Ppost[:,:,k]*C_in'*inv(C_in*Ppost[:,:,k]*C_in'+R_in)
        P[:,:,k]=Ppost[:,:,k]-Kmat[:,:,k]*C_in*Ppost[:,:,k]

        Xhat[:,k]=Xhatpost[:,k]+Kmat[:,:,k]*(Y[:,k]-C_in*Xhatpost[:,k])
    end

    Plim=P[:,:,MAX_TIME]
    Pplus=A_in*Plim*A_in'+Q_in
    K_km=Pplus*C_in'*inv(C_in*Pplus*C_in'+R_in)
    return inv(T)*K_km #

end


## 采用分块SVD，只标准化 unstable states 对应的那些位置
function svd_S(G)
# n_s number of unstable states 后 n_s 个状态是stable的
# println("in svd S, G is ", G)


θ=ones(n_u)
for i=1:n_u
    # @show i, G[:,i]
    if norm(G[:,i], 1)<10^(-8)
        θ[i]=0 # 0列对应0
    end
end
r=Int(sum(θ))
E=Matrix(1.0I, n, n)

for i=1:r # 只对unstable states 的非零列进行操作
    # 第一个零列和第 i个非零列
    zero_index = findall(==(0), θ)
    if !isempty(zero_index) # if there is zero entry
        first_zero_index=zero_index[1]
        first_notzero_index=findall(!=(0), θ)[i]
        # 让此位置之后的第一个非零行换到这个位置
        if first_notzero_index>first_zero_index
            # 交换两列
            temp=E[:,first_zero_index]
            E[:,first_zero_index]=E[:,first_notzero_index]
            E[:,first_notzero_index]=temp
            # 更新theta
            θ[first_notzero_index]=0
            θ[first_zero_index]=1
        end
    end

end

GE=G*E
T1=GE[1:r,1:r]
T2=GE[r+1:end,1:r]
N3=GE[1:r,r+1:end]
N4=GE[r+1:end,r+1:end]
# @show G
# @show GE
# @show T1
if rank(T1)<r || r<1
    println("T1 singular ! SVD can not be performed ! T1 is")
    @show T1
    Ps=Matrix(1.0I,n,n)
    S=G
    # println("in svd S, S directly returned")
else
    Ps=[inv(T1) zeros(r,n-r);  -T2*inv(T1) Matrix(1.0I, n-r, n-r)]
    S=Ps*GE*inv(E)
    S=S.*(broadcast(abs,S).>10^(-8)) # SET LITTEL ERROR TO 0
    # println("in svd S, S is :")
    # @show S
end
return Ps, S #, E
end

## calculate the opt problem cov matrix ###################################################
function get_Mtildeinv() # n_s is the number of stable states
eig_Π = eigen(Λ-K*C*Λ) # TODO: check the eigen values are distinct
# if -10^(-3)<max(eig_Π.values)<0 # small negative values
#     transfer_mat = eigen(Mtinv).vectors
# end
@show eig_Π.values
global Π = Diagonal(eig_Π.values) # TODO 2: image numbers

# @show Π
global V = eig_Π.vectors


# @show K
# @show C
# @show V, rank(V)
# @show V*Π*inv(V)
# @show Λ - K * C * Λ
# @show V*Π*inv(V)-(Λ-K*C*Λ)

# calculate G_i
global G=complex(zeros(n,n,m))
for i in 1:m
    for j in 1:n
        # @show Λ-Π[j,j]*Matrix(1.0I, n, n)
        # @show C[i,:]'*Λ*inv(Λ-Π[j,j]*Matrix(1.0I, n, n))
        G[j,:,i] = (C[i,:]'*Λ*inv(Λ-Π[j,j]*Matrix(1.0I, n, n))) # C[i,:] will be transpose to a column vector in julia
    end
end
# global G_col=zeros(mn,n)
# for i in 1:m
#     G_col[(i-1)*n+1:i*n,:]=G[:,:,i]
# end
# 采用分块SVD，只标准化 unstable states 对应的那些位置
global Pt = complex(zeros(mn, mn)) # Ptilde
S = complex(zeros(mn, n))
for i in 1:m
    Pt[(i-1)*n+1:i*n, (i-1)*n+1:i*n] , S[(i-1)*n+1:i*n, :] = svd_S(G[:,:,i])
end
S=real(S)
# calculate covariance W tilde
global G_C = complex(zeros(mn, n))
for i in 1:m
    G_C[(i-1)*n+1:i*n, :] = G[:,:,i]-ones(n,1)*C[i,:]'
end

global Qt = G_C*Q_in*G_C'+kron(R_in, ones(n,n))
# Qt=((Qt+Qt')/2) #G[j,:,i]
# @show Q
# @show R
# @show Qt
global Πt=kron(Matrix(1.0I,m,m), Π)
# @show Πt
# Wt=dlyap(Πt,Qt)
global Wt = sylvester(-inv(Πt), Matrix(Πt'), inv(Πt)*Qt)
# Wt=(Wt+Wt')/2
# @show Wt
global r=rank(Wt, rtol=10^(-14))
TW=eigen(Wt).vectors # T*Λ*Tinv=Wt
L=sqrt(Diagonal(eigen(Wt).values[mn-r+1:end]))
global D=inv(L)*inv(TW)[mn-r+1:end,:]*inv(Pt)

# global Nc=kron(Matrix(1.0I,m,m),[zeros(n_s,n_u) Matrix(1.0I,n_s,n_s)])
# Mtinv=inv(Pt*Wt*Pt')
# # @show Mtinv-Mtinv'
# global Mtinv=(Mtinv+Mtinv')/2
# global Wc=zeros(m*(n+n_s),m*(n+n_s))
# Wc[1:mn,1:mn]=Mtinv+Nc'*Nc
# Wc[1:mn,mn+1:end]=Nc'
# Wc[mn+1:end,1:mn]=Nc
# Wc[mn+1:end,mn+1:end]=Matrix(1.0I,m*n_s,m*n_s)
# @show eigen(Wc).values
#
# global r=rank(Wc, rtol=10^(-14))
# TW=eigen(Wc).vectors # T*Λ*Tinv=Wt
# global TWL=TW[:,m*(n+n_s)-r+1:end]
# global TWR=inv(TW)[m*(n+n_s)-r+1:end,:]
# global L=sqrt(Diagonal(eigen(Wc).values[m*(n+n_s)-r+1:end]))
# @show TWL*L*L*TWR-Wc
# @show TWL*L*L*TWL'-Wc

# global D=L*TWR
# @show TW*L*L*inv(TW)-Wc
# @show D'*D-inv(Wc)

return Π, S

end ########################################################################################

function initialize_ζ(x0)
    ζ0=complex(zeros(mn,1))
    for i in 1:m
        ζ0[n*(i-1)+1:n*i]=G[:,:,i]*x0
    end
    return ζ0
end

function update_ζ(ζ, y, u)
    return kron(Diagonal(1.0I,m),Π)*ζ + kron(Diagonal(1.0I,m),ones(n,1))*y + G_C*B*u[1]
end

function  solve_opt(ζ, γ, VERBOSE=1) # TODO ignore the unsymmetric of Pt
    global r
    global D
    global Pt
    x = Variable(n, 1)
    μ = ComplexVariable(mn, 1)
    ν = ComplexVariable(mn, 1)
    μ_new = ComplexVariable(r, 1)
    # global r
    # global D
    # x = Variable(n, 1)
    # μ = ComplexVariable(mn, 1)
    # ν = ComplexVariable(mn, 1)
    # big_μ=ComplexVariable(m*(n+n_s), 1)
    # big_μ_l = ComplexVariable(r, 1)
    # big_μ_r = ComplexVariable(r, 1)
    Nc=kron(Matrix(1.0I,m,m),[zeros(n_s,n_u) Matrix(1.0I,n_s,n_s)])
    obj_fun=(sumsquares(real(μ_new))+sumsquares(imag(μ_new)))/2+γ*norm(ν, 1)+real((μ'*Nc'*Nc*S*x)+(μ'*Nc'*Nc*S*x)')+sumsquares(real(Nc*S*x))+sumsquares(imag(Nc*S*x))

    problem = minimize(obj_fun, Pt*ζ==S*x+μ+ν, μ_new==D*μ ) #+γ*norm(ν, 1) +ν
    # problem = minimize((sumsquares(real(big_μ_l))+sumsquares(imag(big_μ_l)))/2, Pt*ζ==S*x+μ, big_μ_l==L'*TWL'*big_μ, big_μ[1:mn]==μ, big_μ[mn+1,end]==Nc*S*x)
    # Solve the problem by calling solve!
    @time solve!(problem, SCS.Optimizer(linear_solver = SCS.DirectSolver, max_iters=100000,verbose=VERBOSE), warmstart=true) #
    # @show big_μ_l.value'*big_μ_l.value-big_μ.value'*Wc*big_μ.value

    # Check the status of the problem
    println("problem status: ", problem.status) # :Optimal, :Infeasible, :unbounded etc.
    println("optimal value: ", problem.optval)
    return x.value, μ.value, ν.value
end
