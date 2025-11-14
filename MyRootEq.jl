module MyRootEq

#with tol=0 I get eps_mach

function bisection(f::Function,a::Number,b::Number;k::Number=0,tol::Number=0)
    g(x)=f(x)-k
    if g(a)*g(b)>0
        println("Error")
        return 
    end
    m=(b+a)/2
    all_ms=[]
    t_start=time()
    err_max=tol+eps()*maximum([abs(a),abs(b)])
    while abs(f(m))>err_max && time()-t_start<120
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

#------------------------------
end