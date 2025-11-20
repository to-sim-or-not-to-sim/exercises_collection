module MyRootEq
include("MyInterpolation.jl")
import .MyInterpolation as my

#with tol=0 I get eps_mach

std_tol=100*eps()
time_lim=30

"""Uses bisection method to find one root of a given function. It has a limit of time (set to 30s) if it doesn't find a root with a certain tollerance (settable using "tol=Float64") and at the end returns the root and all the points analized in the process. By setting "k=Number" it finds the root of the function f(x)-k."""
function bisection(f::Function,a::Number,b::Number;k::Number=0,tol::Float64=std_tol)
    g(x)=f(x)-k
    if g(a)*g(b)>0
        println("Error")
        return 
    end
    m=(b+a)/2
    all_ms=[]
    t_start=time()
    err_max=tol+eps()*maximum([abs(a),abs(b)])
    while abs(f(m))>err_max && time()-t_start<time_lim
        m=(b+a)/2
        push!(all_ms,m)
        if g(a)*g(m)<0
            b=m
        else
            a=m
        end    
    end
    if time()-t_start>=60
        println("Time limit")
    end
    return m,all_ms
end

"""Approximates the derivative of a function."""
function derivative(f::Function;h::Float64=1e-4)
    df(x)=(f(x+h)-f(x-h))/(2*h)
    return df
end

"""Uses newton method to find one root of a given function, using bracketing. It is possible to turn it off by setting "bracketing=false". The starting point can be set as a or b by using set_start="a" or set_start="b". It needs the derivative of the function. It has a time limit (set to 30 s) if it doesn't find a root with a certain tollerance on both the function value and convergence behaviour (settable using "ftol=Float64" and "xtol=Float64") and at the end returns the root and all the points analized in the process. By setting "k=Number" it finds the root of the function f(x)-k."""
function newton(f::Function,df::Function,a::Number,b::Number;k::Number=0,xtol::Float64=std_tol,ftol::Float64=std_tol,bracketing::Bool=true,set_start::String="")
    g(x)=f(x)-k
    if bracketing
        if g(a)*g(b)>0
            bracketing=false
        end
    end
    x_0=b
    x_previous=a
    if set_start=="a"
        x_0=a
        x_previous=b
    elseif set_start=="b"
        x_0=b
        x_previous=a
    else
        if abs(f(a)/df(a))>abs(f(b)/df(b))
            x_0=a
            x_previous=b
        end
    end
    t_start=time()
    all_xs=[]
    while abs(x_previous-x_0)>=xtol || abs(g(x_0))>=ftol || time()-t_start>time_lim
        x_previous=x_0
        x_0-=g(x_0)/df(x_0)
        if bracketing
            if x_0<a || x_0>b
                x_0=(b+a)/2
                if g(a)*g(x_0)<0
                    b=x_0
                else
                    a=x_0
                end
            end
        end
        push!(all_xs,x_0)
    end
    if time()-t_start>=60
        println("Time limit")
    end
    return x_0,all_xs
end

function secant(f::Function,a::Number,b::Number;k::Number=0,xtol::Float64=std_tol,ftol::Float64=std_tol,bracketing::Bool=true)
    g(x)=f(x)-k
    if bracketing
        if g(a)*g(b)>0
            bracketing=false
        end
    end
    all_xs=[a,b]
    t_start=time()
    while abs(all_xs[end]-all_xs[end-1])>=xtol || abs(g(all_xs[end]))>=ftol || time()-t_start>time_lim
        x=all_xs[end]-g(all_xs[end])*(all_xs[end]-all_xs[end-1])/(g(all_xs[end])-g(all_xs[end-1]))
        if bracketing
            if x<a || x>b
                x=(b+a)/2
                if g(a)*g(x)<0
                    b=x
                else
                    a=x
                end
            end
        end
        push!(all_xs,x)
    end
    return all_xs[end],all_xs
end

function inverse_quadratic_interpolation(f::Function,x1::Number,x2::Number,x3::Number;k::Number=0)
    g(x)=f(x)-k
    xs=[x1,x2,x3]
    ys=g.(xs)
    interpol(y)=my.interpolation(y,ys,xs)
    return interpol
end

function interpolation_method(f::Function,x1::Number,x2::Number,x3::Number;k::Number=0,tol::Float64=std_tol)
    g(x)=f(x)-k
    xs=[x1,x2,x3]
    t_start=time()
    while abs(f(xs[end]))>tol && time()-t_start<time_lim
        interpol=inverse_quadratic_interpolation(g,xs[end-2],xs[end-1],xs[end])
        x=interpol(0)
        push!(xs,x)
    end
    return xs[end],xs
end


#------------------------------
end