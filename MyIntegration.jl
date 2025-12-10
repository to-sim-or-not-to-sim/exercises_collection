module MyIntegration
include("MyRootEq.jl")
import .MyRootEq as my
using JSON

function trapezoidal(f::Function,a::Number,b::Number,m::Int64)
    h=(b-a)/m
    sum=(f(a)+f(b))/2
    for j in 1:m-1
        sum+=f(a+j*h)
    end
    return sum*h
end

function simpson(f::Function,a::Number,b::Number,m::Int64)
    h=(b-a)/(2*m)
    sum=f(a)+f(b)
    values=zeros(Float64,m)
    for j in 1:m-1
        sum+=4*f(a+(2*j-1)*h)+2*f(a+2*j*h)
    end
    sum+=4*f(a+(2*m-1)*h)
    return sum*h/3
end

#POLY CLASS-----------------------------------------------------#
                                                                #
struct Poly                                                     #
    coeff::Vector{Float64}                                      #
end                                                             #
                                                                #
function Poly_derivative(p::Poly)                               #
    len=length(p.coeff)                                         #
    der_coeff=zeros(len-1)                                      #
    for i in 1:len-1                                            #
        der_coeff[i]=p.coeff[i+1]*i                             #
    end                                                         #
    der_p=Poly(der_coeff)                                       #
    return der_p                                                #
end                                                             #
                                                                #
function Poly_func(p::Poly)                                     #
    function f(x)                                               #
        len=length(p.coeff)                                     #
        sum=p.coeff[1]                                          #
        for i in 2:len                                          #
            sum+=p.coeff[i]*x^(i-1)                             #
        end                                                     #
        return sum                                              #
    end                                                         #
    return f                                                    #
end                                                             #
                                                                #
import Base: +, -, *                                            #
                                                                #
function +(p::Poly,q::Poly)                                     #
    len_p=length(p.coeff)                                       #
    len_q=length(q.coeff)                                       #
    new_coeff=zeros(maximum([len_p,len_q]))                     #
    if len_p>len_q                                              #
        new_coeff=copy(p.coeff)                                 #
         for i in 1:len_q                                       #
            new_coeff[i]+=q.coeff[i]                            #
        end                                                     #
    else                                                        #
        new_coeff=copy(q.coeff)                                 #
        for i in 1:len_p                                        #
            new_coeff[i]+=p.coeff[i]                            #
        end                                                     #
    end                                                         #
    return Poly(new_coeff)                                      #
end                                                             #
                                                                #
function *(p::Poly,lambda::Number)                              #
    new_coeff=p.coeff.*lambda                                   #
    return Poly(new_coeff)                                      #
end                                                             #
                                                                #
function *(lambda::Number,p::Poly)                              #
    new_coeff=p.coeff.*lambda                                   #
    return Poly(new_coeff)                                      #
end                                                             #
                                                                #
function -(p::Poly,q::Poly)                                     #
    r=-1*q                                                      #
    return p+r                                                  #
end                                                             #
                                                                #
function *(p::Poly,q::Poly)                                     #
    len_p=length(p.coeff)                                       #
    len_q=length(q.coeff)                                       #
    new_coeff=zeros(len_p+len_q-1)                              #
    for i in 1:len_p                                            #
        for j in 1:len_q                                        #
            new_coeff[i+j-1]+=p.coeff[i]*q.coeff[j]             #
        end                                                     #
    end                                                         #
    return Poly(new_coeff)                                      #
end                                                             #
                                                                #
#---------------------------------------------------------------#

function Legendre_poly(k::Int64)
    polys=Vector{Poly}(undef,k+1)
    x=Poly([0,1])
    for i in 0:k
        m=i-1
        if i==0
            polys[i+1]=Poly([1])
        elseif i==1
            polys[i+1]=Poly([0,1])
        else
            polys[i+1]=((2*i-1)/i)*x*polys[i]-((i-1)/i)*polys[i-1]
        end
    end
    return polys 
end

#INITIALIZING LEGENDRE ZEROS AND WEIGHTS----------------------------------------#
function initialize_legendre(;NUM::Int64=100)
    println("Creating file...")
    polys=Legendre_poly(NUM)
    dpolys=Poly_derivative.(polys)
    fs=Poly_func.(polys)
    dfs=Poly_func.(dpolys)
    all_zeros=Vector{Vector{Float64}}(undef,NUM)
    all_w=Vector{Vector{Float64}}(undef,NUM)
    for i in 2:NUM+1
        n=i-1
        println((n/NUM*100)," %")
        f=fs[i]
        df=dfs[i]
        func_zero=zeros(n)
        w=zeros(n)
        for k in 1:n
            phi_k=((4*k-1)/(4*n+2))*pi
            x_k=cos(phi_k)
            func_zero[k]=my.newton(f,df,x_k,0,bracketing=false,set_start="a")[1]
            w[k]=2/((1-func_zero[k]^2)*(df(func_zero[k])^2))
        end
        all_zeros[i-1]=func_zero
        all_w[i-1]=w
    end
    println("Saving...")       
    data=Dict("xk"=>all_zeros,"wk"=>all_w,)
    open("legendre.json","w") do file
        JSON.print(file,data,2)
    end
    println("Done.")
end
#-------------------------------------------------------------------------------#

#initialize_legendre() #<--------gli zeri sono da sistemare!!!
data=JSON.parsefile("legendre.json")
all_zeros=data["xk"]
all_w=data["wk"]

function legendre_gauss(f::Function,a::Number,b::Number,n::Int64)
    g(x)=((b-a)/2)*f(0.5*((b-a)*x+a+b))
    if n>length(all_zeros)
        println("Error, this code works for n up to $NUM.")
        return
    end
    x=all_zeros[n]
    println(x)
    w=all_w[n]
    sum=0
    for i in 1:n
        sum+=w[i]*g(x[i])
    end
    return sum
end

#ADD CLENSHAW-CURTIS

function clenshaw_curtis(f::Function,a::Number,b::Number,n::Int64)
    g(x)=((b-a)/2)*f(0.5*((b-a)*x+a+b))
    k=collect(0:n)
    theta_k=k*pi/n
    x=cos.(theta_k)    
    w=zeros(n+1)
    lim=div(n,2)
    for i in k
        sum=0
        for j in 1:lim
            bj=2
            if j==n/2
                bj=1
            end
            sum+=(bj/(4*j^2-1))*cos(2*j*theta_k[i+1])
        end
        c_k=2
        if i==0 || i==n
            c_k=1
        end
        w[i+1]=(c_k/n)*(1-sum)
    end
    sum=0
    for i in 1:n+1
        sum+=g(x[i])*w[i]
    end
    return sum
end

#USEFUL FUNCTIONS---------------------------------------------------#
                                                                    #
phi_finite(t)=tanh((pi/2)*sinh(t))                                  #
phi_prime_finite(t)=(pi/2)*cosh(t)/(cosh((pi/2)*sinh(t))^2)         #
                                                                    #
phi_semifinite1(t)=exp((pi/2)*sinh(t))                              #
phi_prime_semifinite1(t)=exp((pi/2)*sinh(t))*(pi/2)*cosh(t)         #
                                                                    #
phi_semifinite2(t)=exp(t-exp(-t))                                   #
phi_prime_semifinite2(t)=exp(t-exp(-t))*(1+exp(-t))                 #
                                                                    #
phi_infinite(t)=sinh((pi/2)*sinh(t))                                #
phi_prime_infinite(t)=cosh((pi/2)*sinh(t))*(pi/2)*cosh(t)           #
                                                                    #
#-------------------------------------------------------------------#

plus_inf="pinf"
minus_inf="minf"
#FINITE INTERVAL
#  |  |  |  |  
#  V  V  V  V  
function double_exponential(f::Function,a::Number,b::Number,N::Int64;tol::Float64=my.std_tol)
    g(x)=((b-a)/2)*f(0.5*((b-a)*x+a+b))
    t_M=1
    new_f(t)=g(phi_finite(t))*phi_prime_finite(t)
    while abs(new_f(t_M))>tol || abs(new_f(-t_M))>tol
        t_M+=0.1
    end
    h=t_M/N
    k=collect(-N:N)
    x=phi_finite.(k*h)
    w=h.*phi_prime_finite.(k*h)
    sum=0
    for i in k
        sum+=w[i+N+1]*g(x[i+N+1])
    end
    return sum
end

#SEMIFINITE INTERVAL
#  |  |  |  |  
#  V  V  V  V  
function double_exponential(f::Function,a::Number,b::String,N::Int64;exp_weigth::Bool=false,tol::Float64=my.std_tol)
    g(x)=f(x-a)
    if b!="pinf"
        println("'",b,"' is not a valid input.")
        return
    end
    t_M=1
    if exp_weigth
        new_fw(t)=g(phi_semifinite2(t))*phi_prime_semifinite2(t)
        while abs(new_fw(t_M))>tol || abs(new_fw(-t_M))>tol
            t_M+=0.1
        end    
    else
        new_f(t)=g(phi_semifinite1(t))*phi_prime_semifinite1(t)
        while abs(new_f(t_M))>tol || abs(new_f(-t_M))>tol
            t_M+=0.1
        end
    end
    h=t_M/N
    k=collect(-N:N)
    x=phi_semifinite1.(k*h)
    w=h.*phi_prime_semifinite1.(k*h)
    if exp_weigth
        x=phi_semifinite2.(k*h)
        w=h.*phi_prime_semifinite2.(k*h)
    end
    sum=0
    for i in k
        sum+=w[i+N+1]*g(x[i+N+1])
    end
    return sum
end

function double_exponential(f::Function,a::String,b::Number,N::Int64;exp_weigth::Bool=false,tol::Float64=my.std_tol)
    g(x)=f(-x+b)
    if a!="minf"
        println("'",a,"' is not a valid input.")
        return
    end
    t_M=1
    if exp_weigth
        new_fw(t)=g(phi_semifinite2(t))*phi_prime_semifinite2(t)
        while abs(new_fw(t_M))>tol || abs(new_fw(-t_M))>tol
            t_M+=0.1
        end    
    else
        new_f(t)=g(phi_semifinite1(t))*phi_prime_semifinite1(t)
        while abs(new_f(t_M))>tol || abs(new_f(-t_M))>tol
            t_M+=0.1
        end
    end
    h=t_M/N
    k=collect(-N:N)
    x=phi_semifinite1.(k*h)
    w=h.*phi_prime_semifinite1.(k*h)
    if exp_weigth
        x=phi_semifinite2.(k*h)
        w=h.*phi_prime_semifinite2.(k*h)
    end
    sum=0
    for i in k
        sum+=w[i+N+1]*g(x[i+N+1])
    end
    return sum
end

#INFINITE INTERVAL
#  |  |  |  |  
#  V  V  V  V  
function double_exponential(f::Function,N::Int64;tol::Float64=my.std_tol) 
    t_M=1
    new_f(t)=f(phi_infinite(t))*phi_prime_infinite(t)
    while abs(new_f(t_M))>tol || abs(new_f(-t_M))>tol
        t_M+=0.1
    end
    h=t_M/N
    k=collect(-N:N)
    x=phi_infinite.(k*h)
    w=h.*phi_prime_infinite.(k*h)
    sum=0
    for i in k
        sum+=w[i+N+1]*f(x[i+N+1])
    end
    return sum
end

#---------------------
end