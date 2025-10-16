module MyLinearAlgebra

function diagm(diagonal::Vector{T} where T<:Number)
    n=length(diagonal)
    A=zeros(Float64,n,n)
    for i in 1:n
        A[i,i]=diagonal[i]
    end
    return A
end

function lower_triangular_matrix_solver(A::AbstractMatrix{T} where T<:Number,b::Vector{T} where T<:Number)
    n=length(b)
    x=zeros(Float64,n)
    for i in 1:n
        x_i=b[i]
        for j in 1:i-1
            x_i-=A[i,j]*x[j]
        end
        x[i]=x_i/A[i,i]
    end
    return x
end

function upper_triangular_matrix_solver(A::AbstractMatrix{T} where T<:Number,b::Vector{T} where T<:Number)
    n=length(b)
    x=zeros(Float64,n)
    for i in n:-1:1
        x_i=b[i]
        for j in n:-1:i+1
            x_i-=A[i,j]*x[j]
        end
        x[i]=x_i/A[i,i]
    end
    return x
end   

function lower_triangular_matrix_solver(A::AbstractMatrix{T} where T<:Number,b::AbstractMatrix{T} where T<:Number)
    p,n=size(b)
    x=zeros(Float64,n,p)
    for k in 1:p
        for i in 1:n
            x_ik=b[i,k]
            for j in 1:i-1
                x_ik-=A[i,j]*x[j,k]
            end
            x[i,k]=x_ik/A[i,i]
        end
    end
    return x
end

function upper_triangular_matrix_solver(A::AbstractMatrix{T} where T<:Number,b::AbstractMatrix{T} where T<:Number)
    p,n=size(b)
    x=zeros(Float64,n,p)
    for k in 1:p
        for i in n:-1:1
            x_ik=b[i,k]
            for j in n:-1:i+1
                x_ik-=A[i,j]*x[j,k]
            end
            x[i,k]=x_ik/A[i,i]
        end
    end
    return x
end  

function LU_decomposition(A::AbstractMatrix{T} where T<:Number)
    n,m=size(A)
    L=diagm(ones(n))
    U=zeros(Float64,n,n)
    for i in 1:n-1
        U[i,:]=A[i,:]
        L[:,i]=A[:,i]/U[i,i]
        A=A-L[:,i]*U[i,:]'
    end
    U[n,n]=A[n,n]
    return L, U
end 

function display_LU(A::AbstractMatrix{T} where T<:Number)
    println("Matrix A:")
    display(A)
    L,U=LU_decomposition(A)
    println("Matrix L:")
    display(L)
    println("Matrix U:")
    display(U)
    println("Matrix A-L*U:")
    display(A-L*U)
end

function matrix_solver_LU(A::AbstractMatrix{T} where T<:Number,b::Vector{T} where T<:Number)
    L,U=LU_decomposition(A)
    z=lower_triangular_matrix_solver(L,b)
    x=upper_triangular_matrix_solver(U,z)
    return x
end

function matrix_solver_LU(A::AbstractMatrix{T} where T<:Number,B::AbstractMatrix{T} where T<:Number)
    L,U=LU_decomposition(A)
    n,m=size(B)
    X=zeros(Float64,n,m)
    for i in 1:m
        b=B[:,i]
        z=lower_triangular_matrix_solver(L,b)
        x=upper_triangular_matrix_solver(U,z)
        X[:,i]=x
    end
    return X
end

function display_LU_solver(A::AbstractMatrix{T} where T<:Number,b::Vector{T} where T<:Number)
    x=matrix_solver_LU(A,b)
    println("Matrix A:")
    display(A)
    println("Vector b:")
    display(b)
    println("Results:")
    display(x)
end

function translation_matrix(x::Number,y::Number)
    T=diagm(ones(3))
    T[1,3]=x
    T[2,3]=y
    return T
end

function rotation_matrix(theta::Number)
    c=cos(theta)
    s=sin(theta)
    R=[ c s 0
       -s c 0
        0 0 1]
    return R
end


function det_LU(A::AbstractMatrix{T} where T<:Number)
    L,U=LU_decomposition(A)
    n,m=size(U)
    det=1
    for i in 1:n
        det*=U[i,i]
    end
    return det
end

function inverse(A::AbstractMatrix{T} where T<:Number)
    n,m=size(A)
    identity=diagm(ones(n))
    inverted=matrix_solver_LU(A,identity)
    return inverted
end

end