module MyInterpolation

"""Returns a vector filled with 'len' equally spaced numbers from 'a' to 'b'."""
function linspace(a,b,len)
    x=zeros(len)
    step=(b-a)/(len-1)
    for i in 1:len
        x[i]=a+(i-1)*step
    end
    return x
end

"""Calculates the weighted mean of a set of values with a certain set of weights."""
function weight_mean(values::Vector{T},weights::Vector{T}) where T<:Number
    len=length(values)
    if len==length(weights)
        num=0
        den=0
        for i in 1:len
            num+=weights[i]*values[i]
            den+=weights[i]
        end
        return num/den
    else
        println("Error")
    end
end

"""Calculates the lambdas for the interpolation using generally distributed points."""
function lambdas(x::Vector{T}) where T<:Number
    len=length(x)
    omega=ones(Float64,len)
    for i in 1:len
        for j in 1:len
            if i!=j
                omega[i]*=x[i]-x[j]
            end
        end
    end
    lambda=1 ./omega
    return lambda
end

"""Calculates the factorial of a number. If you want you can truncate the factorial using 'start='."""
function fact(n::Int64;start::Int64=1)
    result=1
    for i in start:n
        result*=i
    end
    return result
end

"""Calculates the binomial coefficient (n i)."""
function binomial_coefficient(n::Int64,i::Int64)
    return fact(n,start=i+1)/fact(n-i)
end

"""Calculates the lambdas for the interpolation using equally distributed points. For large n doesn't work."""
function equi_lambdas(n::Int64)
    lambda=zeros(n)
    lambda[1]=n
    for i in 2:n
        lambda[i]=lambda[i-1]*(i-n)/i
    end
    return lambda
end

"""Computes the value of the interpolation for a function described by the points ('xs','ys') in a point 'x', when 'xs' are generally distributed points, using Lagrange Waring interpolation."""
function interpolation(x::Number,xs::Vector{T},ys::Vector{T}) where T<:Number
    len=length(xs)
    if len!=length(ys)
        println("Error")
        return
    end
    while xs[i]<x
        i+=1
    end
    if xs[i]==x
        x+=1e-12
    end
    lambda=lambdas(xs)
    weights=@. lambda/(x-xs)
    return weight_mean(ys,weights)
end

"""Computes the value of the interpolation for a function described by the points ('xs','ys') in a point 'x', when 'xs' are equally distributed points, using Lagrange Waring interpolation."""
function equi_interpolation(x::Number,a::Number,b::Number,n::Int64,f::Function)
    xs=linspace(a,b,n)
    ys=f.(xs)
    i=1
    while xs[i]<x
        i+=1
    end
    if xs[i]==x
        x+=1e-12
    end
    lambda=lambdas(xs)
    weights=@. lambda/(x-xs)
    return weight_mean(ys,weights)
end

function chebyshev_nodes(n::Int64;a::Number=-1,b::Number=1,kind::Int64=2)
    x=zeros(Float64,n)
    lambda=zeros(Float64,n)
    if kind==1
        for i in 1:n
            x[i]=-cos(((2*i-1)/n)*(pi/2))
            lambda[i]=((-1)^i)*sin(((2*i-1)/n)*(pi/2))
        end
    elseif kind==2
        delta=0.5
        for i in 1:n
            x[i]=-cos((i-1)*pi/(n-1))
            if i>1 && i<n
                delta=1
            elseif i==n
                delta=0.5
            end
            lambda[i]=delta*((-1)^i)
        end
    else
        println("Error")
        return
    end
    x=@. a+(b-a)*(x+1)/2
    return x,lambda
end

function chebyshev_interpolation(x::Number,f::Function,n::Int64;a::Number=-1,b::Number=1,kind::Int64=2)
    xs,lambda=chebyshev_nodes(n,a=a,b=b,kind=kind)
    ys=f.(xs)
    i=1
    while xs[i]<x && i<length(xs)
        i+=1
    end
    if xs[i]==x
        x+=1e-12
    end
    weights=@. lambda/(x-xs)
    return weight_mean(ys,weights)
end

function inf_norm(y::Vector{T},y_fit::Vector{T}) where T<:Number
    len=length(y)
    if len!=length(y_fit)
        println("Error")
        return
    end
    max=abs(y[1]-y_fit[1])
    for i in 2:len
        value=abs(y[i]-y_fit[i])
        if value>max
            max=value
        end
    end
    return max
end

#-----------------------------------------------
end