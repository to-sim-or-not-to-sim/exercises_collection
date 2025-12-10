module MyODE

function euler(fs::AbstractVector,u0s::Vector{T},a::Number,b::Number,n::Int64) where T<:Number
    h=(b-a)/n
    len=length(fs)
    if len!=length(u0s)
        println("Size mismatch.")
        return
    end
    uis=zeros(len,n+1)
    uis[:,1]=u0s
    tis=zeros(n+1)
    tis[1]=a
    for i in 1:n
        tis[i+1]=a+i*h
        for k in 1:len
            uis[k,i+1]=uis[k,i]+h*fs[k](tis[i+1],uis[:,i]...)
        end
    end
    return uis,tis
end

#--------------
end