module MyInterpolation

"""Returns a vector filled with 'len' equally spaced numbers from 'a' to 'b'."""
function linspace(a::Number,b::Number,len::Int64)
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

"""Calculates the lambdas for the interpolation using equally distributed points. For large n doesn't work."""
function equi_lambdas(n::Int64)
    lambda=zeros(n)
    lambda[1]=n
    for i in 2:n
        lambda[i]=lambda[i-1]*(i-1-n)/i
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
    lambda=equi_lambdas(length(xs))
    weights=@. lambda/(x-xs)
    return weight_mean(ys,weights)
end

"""Returns chebyshev nodes and the lambdas for the interpolation."""
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

"""Computes the value of the interpolation for a function described by the points ('xs','ys') in a point 'x', when 'xs' are chebyshev nodes (both first and second kind, using "kind=1" or "kind=2") between a and b, using Lagrange Waring interpolation."""
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

"""Returns the inf-norm of two sets of data."""
function inf_norm(y::Vector{T},y_fit::Vector{T}) where T<:Number
    len=length(y)
    if len!=length(y_fit)
        println("Error")
        return
    end
    return maximum(abs.(y.-y_fit))
end

"""This function expand the chebyshev nodes in [-1,1] on the real axis using x'=2x/(1-x^2)."""
function chebyshev2real(x::Vector{T}) where T<:Number
	return @. 2*x/(1-x^2)
end

"""Computes the value of the interpolation for a function described by the points ('xs','ys') in a point 'x' on the real axis using Lagrange Waring interpolation and chebyshev nodes expanded using "chepyshev2real()" method."""
function real_line_interpolation(z::Number,f::Function,n::Int64)
	xs,lambda=chebyshev_nodes(n,kind=1)
	zs=chebyshev2real(xs)
	ys=f.(zs)
	i=1
    while zs[i]<z && i<length(zs)
        i+=1
    end
    if zs[i]==z
        z+=1e-12
    end
    weights=@. lambda/(z-zs)
    return weight_mean(ys,weights)
end

"""Computes trigonometric interpolation between a and b with N points."""
function trigonometric_interpolation(x::Number,f::Function,N::Int64;a::Number=-1,b::Number=1)
    j=linspace(0,N-1,N)
    xs=@. a+j*(b-a)/N
    ys=f.(xs)
    i=1
    while xs[i]<x && i<length(xs)
        i+=1
    end
    if xs[i]==x
        return ys[i]
    end
    if N%2==1
        tau=@. sin(N*pi*(x-xs)/(b-a))/(N*sin(pi*(x-xs)/(b-a)))
    else
        tau=@. sin(N*pi*(x-xs)/(b-a))/(N*tan(pi*(x-xs)/(b-a)))
    end
    sum=0
    for i in 1:N
        sum+=tau[i]*ys[i]
    end
    return sum
end

#-----------------------------------------------
end
