module MyLinearAlgebra

"""Returns the matrix with the elements of a vector on its diagonal: diagm(diagonal)."""
function diagm(diagonal::Vector{T} where T<:Number)
    n=length(diagonal)
    A=zeros(Float64,n,n)
    for i in 1:n
        A[i,i]=diagonal[i]
    end
    return A
end

"""Returns x where Ax=b, where A is a lower triangular matrix and b is a vector."""
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

"""Returns x where Ax=b, where A is a upper triangular matrix and b is a vector."""
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

"""Returns x where Ax=b, where A is a lower triangular matrix and b is a matrix."""
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

"""Returns x where Ax=b, where A is a upper triangular matrix and b is a matrix."""
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

"""Returns L and U where LU=A with L lower triangular matrix and U upper triangular matrix."""
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

"""Returns x where Ax=b using LU=A, where A is a matrix and b a vector."""
function matrix_solver_LU(A::AbstractMatrix{T} where T<:Number,b::Vector{T} where T<:Number)
    L,U=LU_decomposition(A)
    z=lower_triangular_matrix_solver(L,b)
    x=upper_triangular_matrix_solver(U,z)
    return x
end

"""Returns x where Ax=b using LU=A, where A is a matrix and b a matrix."""
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

"""Returns a 3x3 translation matrix T:\n
\t\t  1 0 x \n
T(x,y)= 0 1 y \n
\t\t  0 0 1 \n
translation_matrix(x,y).
"""
function translation_matrix(x::Number,y::Number)
    T=diagm(ones(3))
    T[1,3]=x
    T[2,3]=y
    return T
end

"""Returns a 3x3 rotation matrix R:\n
\t\t cos(a) sin(a) 0 \n
R(a)=  -sin(a) cos(a) 0 \n
\t\t      0      0 1 \n
rotation_matrix(theta).
"""
function rotation_matrix(theta::Number)
    c=cos(theta)
    s=sin(theta)
    R=[ c s 0
       -s c 0
        0 0 1]
    return R
end

"""Returns the determinant of a matrix A applying LU decomposition."""
function det_LU(A::AbstractMatrix{T} where T<:Number)
    L,U=LU_decomposition(A)
    n,m=size(U)
    det=1
    for i in 1:n
        det*=U[i,i]
    end
    return det
end

"""Returns the inverse of a matrix A applying LU decomposition and solving the system: AX=I."""
function inverse_LU(A::AbstractMatrix{T} where T<:Number)
    n,m=size(A)
    identity=diagm(ones(n))
    inverted=matrix_solver_LU(A,identity)
    return inverted
end

"""Returns P, L and U where LU=PA with L lower triangular matrix, U upper triangular matrix and P permutation matrix. Returns also the permutation number."""
function PLU_decomposition(A::AbstractMatrix{T} where T<:Number)
    n,m=size(A)
    L=diagm(ones(n))
    P=diagm(ones(n))
    U=zeros(Float64,n,n)
    s=0
    indexes=zeros(Int64,n)
    exp_indexes=collect(1:n)
    for i in 1:n
        #---Permutations count---------------------#
        indexes[i]=argmax(abs.(A[:,i]))            #
        if exp_indexes[i]!=indexes[i]              #
            s+=1                                   #
            exp_index=exp_indexes[i]               #
            exp_indexes[i]=indexes[i]              #
            k=i+1                                  #
            while exp_indexes[k]!=indexes[i]       #
                k+=1                               #
            end                                    #
            exp_indexes[k]=exp_index               #
        end                                        #
        #------------------------------------------#
        U[i,:]=A[indexes[i],:]
        L[:,i]=A[:,i]/U[i,i]
        A=A-L[:,i]*U[i,:]'
    end
    P=P[indexes,:]
    L=L[indexes,:]
    return P, L, U, s
end 

function display_PLU(A::AbstractMatrix{T} where T<:Number)
    println("Matrix A:")
    display(A)
    P,L,U,s=PLU_decomposition(A)
    println("Matrix P:")
    display(P)
    println("Matrix L:")
    display(L)
    println("Matrix U:")
    display(U)
    println("Matrix P*A-L*U:")
    display(P*A-L*U)
end

"""Returns x where Ax=b using LU=PA, where A is a matrix and b a vector."""
function matrix_solver_PLU(A::AbstractMatrix{T} where T<:Number,b::Vector{T} where T<:Number)
    P,L,U,s=PLU_decomposition(A)
    z=lower_triangular_matrix_solver(L,P*b)
    x=upper_triangular_matrix_solver(U,z)
    return x
end

"""Returns x where Ax=b using LU=PA, where A is a matrix and b a matrix."""
function matrix_solver_PLU(A::AbstractMatrix{T} where T<:Number,B::AbstractMatrix{T} where T<:Number)
    P,L,U,s=PLU_decomposition(A)
    n,m=size(B)
    X=zeros(Float64,n,m)
    for i in 1:m
        b=B[:,i]
        z=lower_triangular_matrix_solver(L,P*b)
        x=upper_triangular_matrix_solver(U,z)
        X[:,i]=x
    end
    return X
end

"""Returns the determinant of a matrix A applying LU decomposition and the permutation matrix P."""
function det_PLU(A::AbstractMatrix{T} where T<:Number)
    P,L,U,s=PLU_decomposition(A)
    n,m=size(U)
    det=1
    for i in 1:n
        det*=U[i,i]
    end
    return det*((-1)^s)
end

"""Returns the inverse of a matrix A solving the system: AX=I, applying LU decomposition and the permutation matrix P."""
function inverse_PLU(A::AbstractMatrix{T} where T<:Number)
    n,m=size(A)
    identity=diagm(ones(n))
    inverted=matrix_solver_PLU(A,identity)
    return inverted
end

"""Evaluates a vector's maximum value."""
function max(v::Vector{T} where T<:Number)
    maximum=v[1]
    len=length(v)
    for i in 2:len
        if v[i]>maximum
            maximum=v[i]
        end
    end
    return maximum
end

"""Returns the 1-norm of a vector v."""
function one_norm(v::Vector{T} where T<:Number)
    norm=0
    len=length(v)
    for i in 1:len
        norm+=abs(v[i])
    end
    return norm
end

"""Returns the 1-norm of a matrix A."""
function one_norm(A::AbstractMatrix{T} where T<:Number)
    n,m=size(A)
    norm=zeros(m)
    for j in 1:m
        for i in 1:n
            norm[j]+=abs(A[i,j])
        end
    end
    return max(norm)
end

"""Evaluates the condition number of a square matrix A using 1-norm."""
function one_norm_matrix_condition_number(A::AbstractMatrix{T} where T<:Number)
    invA=inverse_PLU(A)
    norm_A=one_norm(A)
    norm_invA=one_norm(invA)
    return norm_A*norm_invA
end

"""Returns the R matrix for cholensky decomposition defined as A=R^T R."""
function cholesky_decomposition(A::AbstractMatrix{T} where T<:Number)
    n,m=size(A)
    R=zeros(Float64,n,m)
    for i in 1:m
        for j in 1:i
            element=0
            sum=0
            if i==j
                for k in 1:i-1
                    sum+=R[k,i]^2
                end
                if A[j,i]-sum>=0
                    element=sqrt(A[j,i]-sum)
                end
            else
                for k in 1:i-1
                    sum+=R[k,i]*R[k,j]
                end
                element=(A[j,i]-sum)/R[j,j]
            end
            R[j,i]=element
        end
    end
    return R
end

"""Solves the linear system caused by a positive definite matrix A: Ax=b."""
function matrix_solver_cholesky(A::AbstractMatrix{T} where T<:Number,b::Vector{T} where T<:Number)
    R=cholesky_decomposition(A)
    u=lower_triangular_matrix_solver(R',b)
    x=upper_triangular_matrix_solver(R,u)
    return x
end

"""Says if a matrix is positive definite."""
function is_positive_definite(A::AbstractMatrix{T} where T<:Number,threshold=1e-15)
    n=size(A,1)
    b_test=ones(n)
    x_PLU=matrix_solver_PLU(A,b_test)
    x_chole=matrix_solver_cholesky(A,b_test)
    if all(abs.(x_chole-x_PLU).<=threshold)
        return true
    end
    return false
end

"""Solves the linear system Nx=z with N=A^tA and z=A^tb."""
function linear_leastsquares(A::AbstractMatrix{T} where T<:Number,b::Vector{T} where T<:Number)
    At=A'
    N=At*A
    z=At*b
    x=matrix_solver_cholesky(N,z)
    return x
end

"""Returns a matrix filled with the standard polynomials of a set of data."""
function poly_matrix(set::Vector{T} where T<:Number,n::Int64)
    m=length(set)
    A=zeros(m,n)
    for i in 1:m
        for j in 1:n
            A[i,j]=set[i]^(j-1)
        end
    end
    return A
end

"""Returns a matrix filled with sin and cos of a set of data."""
function fourier_matrix(set::Vector{T} where T<:Number,n::Int64)
    m=length(set)
    p=1+2*n
    A=zeros(m,p)
    for i in 1:m
        for j in 1:p
            value=1
            if j==1
                value=1
            elseif j%2==1
                k=(j-1)/2
                value=sin(k*set[i])
            else
                k=j/2
                value=cos(k*set[i])
            end
            A[i,j]=value
        end
    end
    return A
end

end