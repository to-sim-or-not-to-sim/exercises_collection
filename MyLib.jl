module MyLib

    #--- LINEAR SYSTEMS -----------------------------------------------------------------------------------------------

    module MyLinearAlgebra
    export diagm, lower_triangular_matrix_solver, upper_triangular_matrix_solver, LU_decomposition,
           matrix_solver_LU, translation_matrix, rotation_matrix, det_LU, inverse_LU, PLU_decomposition,
           matrix_solver_PLU, det_PLU, inverse_PLU, one_norm, one_norm_matrix_condition_number,
           cholesky_decomposition, matrix_solver_cholesky, is_positive_definite, linear_leastsquares,
           Vandermonde_matrix, fourier_matrix, poly_linear_fit, fourier_linear_fit, norm, QR_decomposition,
           matrix_solver_QR, matrix_solver_Qless, eigenvalues

        using PyPlot
        std_tol=100*eps()

        """Returns the matrix with the elements of a vector on its diagonal: diagm(diagonal)."""
        function diagm(diagonal::Vector{T}) where T<:Number
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
        function LU_decomposition(A::AbstractMatrix{T}) where T<:Number
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

        function display_LU(A::AbstractMatrix{T}) where T<:Number
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
        function det_LU(A::AbstractMatrix{T}) where T<:Number
            L,U=LU_decomposition(A)
            n,m=size(U)
            det=1
            for i in 1:n
                det*=U[i,i]
            end
            return det
        end

        """Returns the inverse of a matrix A applying LU decomposition and solving the system: AX=I."""
        function inverse_LU(A::AbstractMatrix{T}) where T<:Number
            n,m=size(A)
            identity=diagm(ones(n))
            inverted=matrix_solver_LU(A,identity)
            return inverted
        end

        """Returns P, L and U where LU=PA with L lower triangular matrix, U upper triangular matrix and P permutation matrix. Returns also the permutation number."""
        function PLU_decomposition(A::AbstractMatrix{T}) where T<:Number
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

        function display_PLU(A::AbstractMatrix{T}) where T<:Number
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
        function det_PLU(A::AbstractMatrix{T}) where T<:Number
            P,L,U,s=PLU_decomposition(A)
            n,m=size(U)
            det=1
            for i in 1:n
                det*=U[i,i]
            end
            return det*((-1)^s)
        end

        """Returns the inverse of a matrix A solving the system: AX=I, applying LU decomposition and the permutation matrix P."""
        function inverse_PLU(A::AbstractMatrix{T}) where T<:Number
            n,m=size(A)
            identity=diagm(ones(n))
            inverted=matrix_solver_PLU(A,identity)
            return inverted
        end

        """Evaluates a vector's maximum value."""
        function max(v::Vector{T}) where T<:Number
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
        function one_norm(v::Vector{T}) where T<:Number
            norm=0
            len=length(v)
            for i in 1:len
                norm+=abs(v[i])
            end
            return norm
        end

        """Returns the 1-norm of a matrix A."""
        function one_norm(A::AbstractMatrix{T}) where T<:Number
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
        function one_norm_matrix_condition_number(A::AbstractMatrix{T}) where T<:Number
            invA=inverse_PLU(A)
            norm_A=one_norm(A)
            norm_invA=one_norm(invA)
            return norm_A*norm_invA
        end

        """Returns the R matrix for cholensky decomposition defined as A=R^T R."""
        function cholesky_decomposition(A::AbstractMatrix{T}) where T<:Number
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
        function is_positive_definite(A::AbstractMatrix{T};tol::Float64=std_tol) where T<:Number
            R=cholesky_decomposition(A)
            if all(R[i,i]>tol for i in 1:size(R)[1])
                return true
            end
            return false
        end

        """Solves the linear system Nx=z with N=A^tA and z=A^tb."""
        function linear_leastsquares(A::AbstractMatrix{T} where T<:Number,b::Vector{T} where T<:Number;cholensky::Bool=false)
            if cholensky
                At=A'
                N=At*A
                z=At*b            
                x=matrix_solver_cholesky(N,z)
                return x
            end
            x=matrix_solver_QR(A,b)
            return x
        end

        """Returns a matrix filled with the standard polynomials of a set of data."""
        function Vandermonde_matrix(set::Vector{T},n::Int64) where T<:Number
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
        function fourier_matrix(set::Vector{T},n::Int64) where T<:Number
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

        function plot_matrix(x::Vector{T} where T<:Number,A::AbstractMatrix{T} where T<:Number,z::Vector{T} where T<:Number,color::String,label::String) 
            plot(x,A*z,color=color,label=label)
        end

        """Evaluates a linear fit with polynomials for a function f given a certain number of point x, y=f(x). Eventually can plot it."""
        function poly_linear_fit(x::Vector{T} where T<:Number,y::Vector{T} where T<:Number,n::Int64;plot::Bool=false,color::String="blue",label::String="Polynomial fit",cholensky::Bool=false) 
            A=Vandermonde_matrix(x,n)
            z=linear_leastsquares(A,y,cholensky=cholensky)
            if plot
                plot_matrix(x,A,z,color,label)
            end
            return z
        end

        """Evaluates a linear fit with sin(kx) and cos(kx) for a function f given a certain number of point x, y=f(x). Eventually can plot it."""
        function fourier_linear_fit(x::Vector{T} where T<:Number,y::Vector{T} where T<:Number,n::Int64;plot::Bool=false,color::String="blue",label::String="Fourier fit",cholensky::Bool=false) 
            A=fourier_matrix(x,n)
            z=linear_leastsquares(A,y,cholensky=cholensky)
            if plot
                plot_matrix(x,A,z,color,label)
            end
            return z
        end

        """Calculates the norm of a vector."""
        function norm(x::Vector{T}) where T<:Number
            n=length(x)
            sum=0
            for i in 1:n
                sum+=x[i]^2
            end
            return sqrt(sum)
        end

        """Applys the QR decomposition to a matrix A."""
        function QR_decomposition(A::AbstractMatrix{T}) where T<:Number
            m,n=size(A)
            Q=zeros(Float64,m,n)
            R=zeros(Float64,n,n)
            for j in 1:n
                v_j=A[:,j]
                for i in 1:j-1
                    R[i,j]=Q[:,i]'*v_j
                    v_j=v_j-Q[:,i]*R[i,j] 
                end
                R[j,j]=norm(v_j)
                if R[j,j]!=0
                    Q[:,j]=v_j/R[j,j]
                else
                    Q[:,j]=zeros(length(v_j))
                end
            end
            return Q,R
        end

        """Solves a linear system using QR decomposition."""
        function matrix_solver_QR(A::AbstractMatrix{T} where T<:Number,b::Vector{T} where T<:Number)
            Q,R=QR_decomposition(A)
            z=Q'*b
            x=upper_triangular_matrix_solver(R,z)
            return x
        end

        """Solves a linear system using QR decomposition and Qless algorithm."""
        function matrix_solver_Qless(A::AbstractMatrix{T} where T<:Number,b::Vector{T} where T<:Number)
            n,m=size(A)
            A2=zeros(n,m+1)
            A2[:,1:m]=A[:,:]
            A2[:,m+1]=b
            Q1,R1=QR_decomposition(A2)
            z=R1[1:m,m+1]
            R=R1[1:m,1:m]
            x=upper_triangular_matrix_solver(R,z)
            return x
        end

        """Returns a try to diagonalize the matrix and the matrix with eigenvectors."""
        function eigenvalues(A::AbstractMatrix{T};k::Int64=100) where T<:Number
            n,m=size(A)
            newQ=diagm(ones(n))
            newA=copy(A)
            for i in 1:k
                Q,R=QR_decomposition(newA)
                newA=R*Q
                newQ*=Q
            end
            return newA,newQ
        end

    end

    #------------------------------------------------------------------------------------------------------------------

    #--- INTERPOLATION ------------------------------------------------------------------------------------------------

    module MyInterpolation
    export linspace, interpolation, inf_norm, lambdas, equi_lambdas, equi_interpolation,
           newton_polynomials, newton_interpolation, chebyshev_nodes, chebyshev_interpolation,
           real_line_interpolation, trigonometric_interpolation

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
                lambda[i]=lambda[i-1]*(i-1-n)/(i-1)
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
            i=1
            while xs[i]<x && i<length(xs)
                i+=1
            end
            if xs[i]==x
                return f(x)
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
            while xs[i]<x && i<length(xs)
                i+=1
            end
            if xs[i]==x
                return f(x)
            end
            lambda=equi_lambdas(length(xs))
            weights=@. lambda/(x-xs)
            return weight_mean(ys,weights)
        end

        """Calculates the divided differences for a set of data xs and a given function f."""
        function divided_differences(xs::Vector{T},f::Function) where T<:Number
            len=length(xs)
            if len==1
                return f(xs[1])
            elseif len==2
                return (f(xs[2])-f(xs[1]))/(xs[2]-xs[1])
            else
                x_1=zeros(len-1)
                x_2=zeros(len-1)
                x_1[:]=xs[1:len-1]
                x_2[:]=xs[2:len]
                return (divided_differences(x_1,f)-divided_differences(x_2,f))/(xs[1]-xs[end])
            end
            return
        end

        """Evaluates a Newton polynomial in x."""
        function newton_polynomials(x::Number,xs::Vector{T}) where T<:Number
            len=length(xs)
            poly=1
            for i in 2:len
                poly*=(x-xs[i-1])
            end
            return poly
        end

        """Evaluates Newton interpolation in x."""
        function newton_interpolation(x::Number,xs::Vector{T},f::Function) where T<:Number
            p=0
            len=length(xs)
            for k in 1:len
                xs_k=xs[1:k]
                c_k=divided_differences(xs_k,f)
                phi_k=newton_polynomials(x,xs_k)
                p+=c_k*phi_k
            end
            return p 
        end

        """Returns the chebyshev nodes in [a,b] and their weights. It is possible to select the kind of nodes."""
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
                    elseif i==n || i==1
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

        """Evaluate the Lagrange-Waring interpolation of a function in x, done using chebyshev nodes. It is possible to select the kind of the nodes."""
        function chebyshev_interpolation(x::Number,f::Function,n::Int64;a::Number=-1,b::Number=1,kind::Int64=2)
            xs,lambda=chebyshev_nodes(n,a=a,b=b,kind=kind)
            ys=f.(xs)
            i=1
            while xs[i]<x && i<length(xs)
                i+=1
            end
            if xs[i]==x
                return f(x)
            end
            weights=@. lambda/(x-xs)
            return weight_mean(ys,weights)
        end

        """Returns inf-norm of the difference between two vectors."""
        function inf_norm(y::Vector{T},y_fit::Vector{T}) where T<:Number
            len=length(y)
            if len!=length(y_fit)
                println("Error")
                return
            end
            return maximum(abs.(y.-y_fit))
        end

        """Useful to transform chebyshev points over the real axis."""
        function chebyshev2real(x::Vector{T}) where T<:Number
            return @. 2*x/(1-x^2)
        end

        """Evaluates the interpolation over a real axis in z."""
        function real_line_interpolation(z::Number,f::Function,n::Int64)
            xs,lambda=chebyshev_nodes(n,kind=1)
            zs=chebyshev2real(xs)
            ys=f.(zs)
            i=1
            while zs[i]<z && i<length(zs)
                i+=1
            end
            if zs[i]==z
                return f(z)
            end
            weights=@. lambda/(z-zs)
            return weight_mean(ys,weights)
        end

        """Evaluate trigonometric interpolation of a function in x. N can be both odd and even."""
        function trigonometric_interpolation(x::Number,f::Function,N::Int64;a::Number=-1,b::Number=1)
            k=linspace(0,N,N+1)
            xk=@. a+k*(b-a)/N
            yk=f.(xk)
            len=length(xk)
            for i in 1:len
                if xk[i]==x
                    return yk[i]
                end
            end
            tauk=zeros(len)
            if N%2==1
                tauk=@. sin(N*pi*(x-xk)/(b-a))/(N*sin(pi*(x-xk)/(b-a)))
            else
                tauk=@. sin(N*pi*(x-xk)/(b-a))/(N*tan(pi*(x-xk)/(b-a)))
            end
            p=0
            for i in 1:N
                p+=tauk[i]*yk[i]
            end
            return p
        end
        
    end
    #------------------------------------------------------------------------------------------------------------------

    #--- ROOT FINDING -------------------------------------------------------------------------------------------------

    module MyRootFinding
    export bisection, newton, secant, interpolation_method, multiple_newton

        using ..MyInterpolation: interpolation

        std_tol=100*eps()
        std_rep_lim=10000

        """Uses bisection method to find one root of a given function. It has a limit of repetitions (settable using "rep_lim=Int64") if it doesn't find a root with a certain tollerance (settable using "tol=Float64") and at the end returns the root and all the points analized in the process. By setting "k=Number" it finds the root of the function f(x)-k."""
        function bisection(f::Function,a::Number,b::Number;k::Number=0,tol::Float64=std_tol,rep_lim::Int64=std_rep_lim)
            g(x)=f(x)-k
            if g(a)*g(b)>0
                println("Error")
                return 
            end
            m=(b+a)/2
            all_ms=[]
            err_max=tol+eps()*maximum([abs(a),abs(b)])
            rep=0
            while abs(f(m))>err_max && rep<rep_lim
                rep+=1
                m=(b+a)/2
                push!(all_ms,m)
                if g(a)*g(m)<0
                    b=m
                else
                    a=m
                end    
            end
            if rep>=rep_lim
                println("Repetition limit.")
            end
            return m,all_ms
        end

        #=
        """Approximates the derivative of a function."""
        function derivative(f::Function;h::Float64=1e-4)
            df(x)=(f(x+h)-f(x-h))/(2*h)
            return df
        end
        =#

        """Uses newton method to find one root of a given function, using bracketing. It is possible to turn it off by setting "bracketing=false". The starting point can be set as a or b by using set_start="a" or set_start="b". It needs the derivative of the function. It has a limit of repetitions (settable using "rep_lim=Int64"); if it doesn't find a root with a certain tollerance on both the function value and convergence behaviour (settable using "ftol=Float64" and "xtol=Float64") at the end returns the last guess and all the points analized in the process. By setting "k=Number" it finds the root of the function f(x)-k."""
        function newton(f::Function,df::Function,a::Number,b::Number;k::Number=0,xtol::Float64=std_tol,ftol::Float64=std_tol,bracketing::Bool=true,set_start::String="",rep_lim::Int64=std_rep_lim)
            g(x)=f(x)-k
            if bracketing
                if g(a)*g(b)>0
                    println("Bracketing is not possible.")
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
                if abs(g(a)/df(a))>abs(g(b)/df(b))
                    x_0=a
                    x_previous=b
                end
            end
            all_xs=[Float64(x_0)]
            rep=0
            while (abs(x_previous-x_0)>=xtol && abs(g(x_0))>=ftol) && rep<rep_lim
                rep+=1
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
            if rep>=rep_lim
                println("Repetition limit.")
            end
            return x_0,all_xs
        end

        """Uses newton method to find one root of a given function with an unknown multiplicity, using bracketing. It is possible to turn it off by setting "bracketing=false". The starting point can be set as a or b by using set_start="a" or set_start="b". It needs the derivative and the second derivative of the function. It has a limit of repetitions (settable using "rep_lim=Int64"); if it doesn't find a root with a certain tollerance on both the function value and convergence behaviour (settable using "ftol=Float64" and "xtol=Float64") at the end returns the last guess and all the points analized in the process. By setting "k=Number" it finds the root of the function f(x)-k."""
        function multiple_newton(f::Function,df::Function,d2f::Function,a::Number,b::Number;k::Number=0,xtol::Float64=std_tol,ftol::Float64=std_tol,bracketing::Bool=true,set_start::String="",rep_lim::Int64=std_rep_lim)
            g(x)=f(x)-k
            if bracketing
                if g(a)*g(b)>0
                    println("Bracketing is not possible.")
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
            all_xs=[Float64(x_0)]
            rep=0
            while (abs(x_previous-x_0)>=xtol && abs(g(x_0))>=ftol) && rep<rep_lim
                rep+=1
                x_previous=x_0
                fx=g(x_0)
                dfx=df(x_0)
                d2fx=d2f(x_0)
                x_0-=fx/(dfx-fx*d2fx/dfx)
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
            if rep>=rep_lim
                println("Repetition limit.")
            end
            return x_0,all_xs
        end

        """Uses secant method to find one root of a given function, using bracketing. It is possible to turn it off by setting "bracketing=false". It has a limit of repetitions (settable using "rep_lim=Int64"); if it doesn't find a root with a certain tollerance on both the function value and convergence behaviour (settable using "ftol=Float64" and "xtol=Float64") at the end returns the last guess and all the points analized in the process. By setting "k=Number" it finds the root of the function f(x)-k."""
        function secant(f::Function,a::Number,b::Number;k::Number=0,xtol::Float64=std_tol,ftol::Float64=std_tol,bracketing::Bool=true,rep_lim::Int64=std_rep_lim)
            g(x)=f(x)-k
            if bracketing
                if g(a)*g(b)>0
                    bracketing=false
                end
            end
            all_xs=[a,b]
            rep=0
            while (abs(all_xs[end]-all_xs[end-1])>=xtol && abs(g(all_xs[end]))>=ftol) && rep<rep_lim
                rep+=1
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

        """Interpolates three points with a parabola: x=ay^2+by+c."""
        function inverse_quadratic_interpolation(f::Function,x1::Number,x2::Number,x3::Number;k::Number=0)
            g(x)=f(x)-k
            xs=[x1,x2,x3]
            ys=g.(xs)
            interpol(y)=interpolation(y,ys,xs)
            return interpol
        end

        """Uses inverse quadratic interpolation method to find one root of a given function. It has a limit of repetitions (settable using "rep_lim=Int64"); if it doesn't find a root with a certain tollerance on both the function value and convergence behaviour (settable using "ftol=Float64" and "xtol=Float64") at the end returns the last guess and all the points analized in the process. By setting "k=Number" it finds the root of the function f(x)-k."""
        function interpolation_method(f::Function,x1::Number,x2::Number,x3::Number;k::Number=0,tol::Float64=std_tol,rep_lim::Int64=std_rep_lim)
            g(x)=f(x)-k
            xs=[x1,x2,x3]
            rep=0
            while abs(g(xs[end]))>tol && rep<rep_lim
                rep+=1
                interpol=inverse_quadratic_interpolation(g,xs[end-2],xs[end-1],xs[end])
                x=interpol(0)
                push!(xs,x)
            end
            return xs[end],xs
        end

    end

    #------------------------------------------------------------------------------------------------------------------

    #--- INTEGRATION --------------------------------------------------------------------------------------------------

    module MyIntegration
    export trapezoidal, simpson, legendre_gauss, clenshaw_curtis, double_exponential 

        using JSON

        std_tol=100*eps()

        """Approximate the integral of a function f over the interval [a,b] using trapezoidal method with m nodes."""
        function trapezoidal(f::Function,a::Number,b::Number,m::Int64)
            h=(b-a)/m
            sum=(f(a)+f(b))/2
            for j in 1:m-1
                sum+=f(a+j*h)
            end
            return sum*h
        end

        """Approximate the integral of a function f over the interval [a,b] using simpson rule with m nodes."""
        function simpson(f::Function,a::Number,b::Number,m::Int64)
            h=(b-a)/(2*m)
            sum=f(a)+f(b)
            for j in 1:m-1
                sum+=4*f(a+(2*j-1)*h)+2*f(a+2*j*h)
            end
            sum+=4*f(a+(2*m-1)*h)
            return sum*h/3
        end

        """Evaluates n-th Legendre polynomial in x."""
        function legendre(x::Number,n::Int64)
            P_0=1
            P_1=x
            if n==0
                return 1,0
            end
            for m in 1:n-1
                P_new=((2*m+1)/(m+1))*x*P_1-(m/(m+1))*P_0
                P_0=P_1
                P_1=P_new
            end
            P_prime=n*(x*P_1-P_0)/(x^2-1)
            return P_1,P_prime
        end

        """Creates a file containing the zeros of the first Legendre Polynomials and the relative weights."""
        function initialize_legendre(;NUM::Int64=100)
            println("Creating file...")
            all_zeros=Vector{Vector{Float64}}(undef,NUM)
            all_w=Vector{Vector{Float64}}(undef,NUM)
            for i in 2:NUM+1
                n=i-1
                println((n/NUM*100)," %")
                func_zero=zeros(n)
                w=zeros(n)
                for k in 1:n
                    phi_k=((4*k-1)/(4*n+2))*pi
                    x_k=cos(phi_k)
                    P_1,P_prime=legendre(x_k,n)
                    count=0
                    while abs(P_1)>100*eps() && count<10000
                        count+=1
                        x_k-=P_1/P_prime
                        P_1,P_prime=legendre(x_k,n)
                    end
                    func_zero[k]=x_k
                    w[k]=2/((1-func_zero[k]^2)*(P_prime^2))
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

        #initialize_legendre() #USE THIS THE FIRST TIME YOU OPEN THE FILE

        # --- Load the legendre polynomials zeros and weights
        data=JSON.parsefile("legendre.json")
        all_zeros=data["xk"]
        all_w=data["wk"]
        # ---------------------------------------------------

        """Approximate the integral of a function f over the interval [a,b] using gauus-legendre method with n-th order Legendre polynomial."""
        function legendre_gauss(f::Function,a::Number,b::Number,n::Int64)
            g(x)=((b-a)/2)*f(0.5*((b-a)*x+a+b))
            if n>length(all_zeros)
                println("Error, this code works for n up to $NUM.")
                return
            end
            x=all_zeros[n]
            w=all_w[n]
            sum=0
            for i in 1:n
                sum+=w[i]*g(x[i])
            end
            return sum
        end

        """Approximate the integral of a function f over the interval [a,b] using clenshaw-curtis rule with n nodes."""
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
        """Approximate the integral of a function f over the interval [a,b] using double exponential quadrature method with N nodes."""
        function double_exponential(f::Function,a::Number,b::Number,N::Int64;tol::Float64=eps())
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
            result=sum(g.(x).*w)
            if isnan(result)
                mask=(abs.(x) .< 1.0)
                x=x[mask]
                w=w[mask]
                result=sum(g.(x).*w)
            end
            return result
        end

        #SEMIFINITE INTERVAL
        #  |  |  |  |  
        #  V  V  V  V  
        """Approximate the integral of a function f over the interval [a,+inf] using double exponential quadrature method with N nodes."""
        function double_exponential(f::Function,a::Number,b::String,N::Int64;exp_weigth::Bool=false,tol::Float64=std_tol)
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

        """Approximate the integral of a function f over the interval [-inf,b] using double exponential quadrature method with N nodes."""
        function double_exponential(f::Function,a::String,b::Number,N::Int64;exp_weigth::Bool=false,tol::Float64=std_tol)
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
        """Approximate the integral of a function f over the real line using double exponential quadrature method with N nodes."""
        function double_exponential(f::Function,N::Int64;tol::Float64=std_tol) 
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

    end

    #------------------------------------------------------------------------------------------------------------------

    #--- ODE ----------------------------------------------------------------------------------------------------------

    module MyODE
    export euler, backward_euler, IE2, RK4, BS23, RKDP
    
        using ..MyRootFinding: secant
        using ..MyInterpolation: inf_norm

        """Uses Euler's method to approximate the solution of an initial values problem over [a,b]."""
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
                    uis[k,i+1]=uis[k,i]+h*fs[k](tis[i],uis[:,i]...)
                end
            end
            return uis,tis
        end

        """Uses backward Euler's method to approximate the solution of an initial values problem over [a,b]."""
        function backward_euler(f::Function,u0::Number,a::Number,b::Number,n::Int64)
            h=(b-a)/n
            ui=zeros(n+1)
            ui[1]=u0
            ti=zeros(n+1)
            ti[1]=a
            for i in 1:n
                ti[i+1]=a+i*h
                g(u)=ui[i]+h*f(ti[i+1],u)-u
                initial_guess=ui[i]+h*f(ti[i+1],ui[i])
                ui[i+1]=secant(g,ui[i],initial_guess,bracketing=false)[1]
            end
            return ui,ti
        end

        """Uses IE2 method to approximate the solution of an initial values problem over [a,b]."""
        function IE2(fs::AbstractVector,u0s::Vector{T},a::Number,b::Number,n::Int64) where T<:Number
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
                k1=zeros(len)
                v=zeros(len)
                k2=zeros(len)
                for k in 1:len
                    k1[k]=h*fs[k](tis[i],uis[:,i]...)
                end
                v=@. uis[:,i]+0.5*k1
                for k in 1:len
                    k2[k]=h*fs[k](tis[i]+0.5*h,v...)
                end
                uis[:,i+1]=@. uis[:,i]+k2      
            end
            return uis,tis
        end

        """Uses RK4 method to approximate the solution of an initial values problem over [a,b]."""
        function RK4(fs::AbstractVector,u0s::Vector{T},a::Number,b::Number,n::Int64) where T<:Number
            h=(b-a)/n
            len=length(fs)
            if len!=length(u0s)
                println("Size mismatch.")
                return
            end
            C=[0,0.5,0.5,1]
            A=[0,0.5,0.5,1]
            B=[1/6,1/3,1/3,1/6]
            uis=zeros(len,n+1)
            uis[:,1]=u0s
            tis=zeros(n+1)
            tis[1]=a
            for i in 1:n
                tis[i+1]=a+i*h
                ks=zeros(5,len)
                for j in 1:4
                    for k in 1:len
                        ks[j+1,k]=h*fs[k](tis[i]+C[j]*h,uis[:,i].+A[j]*ks[j,k]...)
                    end
                end
                for j in 1:len
                    uis[j,i+1]=uis[j,i]+sum(B.*ks[2:5,j])
                end
            end
            return uis,tis
        end

        """Uses B2S3 method to approximate the solution of an initial values problem over [a,b]."""
        function BS23(fs::AbstractVector,u0s::Vector{T},a::Number,b::Number;tol::Float64=1e-12) where T<:Number
            h_max=(b-a)/1000
            len=length(fs)
            if len!=length(u0s)
                println("Size mismatch.")
                return
            end
            h=tol^(1/3)/2
            C=[0,0.25,0.75,1]
            A=[   0    0   0   0
                0.5    0   0   0
                  0 0.75   0   0
                2/9  1/3 4/9   0]
            B_3=[2/9,1/3,4/9,0]
            B_2=[7/24,0.25,1/3,0.125]
            uis=[]
            for i in 1:len
                push!(uis,Float64[u0s[i]])
            end
            ti=Float64[]
            push!(ti,a)
            while ti[end]<b
                ks=zeros(4,len) 
                for i in 1:4
                    t_step=ti[end]+C[i]*h
                    u_step=zeros(len)
                    for l in 1:len
                        u_step[l]=sum(A[i,:].*ks[:,l])+uis[l][end]                    
                    end
                    for l in 1:len
                        f=fs[l]
                        ks[i,l]=h*f(t_step,u_step...)
                    end
                end
                u_2=zeros(len)
                u_3=zeros(len)
                for l in 1:len
                    u_2[l]=sum(B_2.*ks[:,l])+uis[l][end]
                    u_3[l]=sum(B_3.*ks[:,l])+uis[l][end]
                end
                errs=abs.(u_2-u_3)
                err=maximum(errs)
                if err<tol
                    for l in 1:len
                        push!(uis[l],u_3[l])
                    end
                    push!(ti,ti[end]+h) 
                end
                h=0.8*h*((tol/err)^(1/3))
                h=minimum([h,h_max])
                if ti[end]+h==ti[end]
                    println("Error during execution: step zero.")
                    return uis,ti
                end 
                if ti[end]+h>b
                    h=b-ti[end]
                end
            end
            return uis,ti
        end

        """Uses RKDP method to approximate the solution of an initial values problem over [a,b]."""
        function RKDP(fs::AbstractVector,u0s::Vector{T},a::Number,b::Number;tol::Float64=1e-12) where T<:Number
            h_max=(b-a)/1000
            len=length(fs)
            if len!=length(u0s)
                println("Size mismatch.")
                return
            end
            h=tol^(1/5)/2
            C=[0,0.2,0.3,0.8,8/9,1,1]
            A=[    0 0 0 0 0 0 0
                 0.2 0 0 0 0 0 0
                3/40 9/40 0 0 0 0 0
               44/45  -56/15 32/9 0 0 0 0
               19372/6561 -25360/2187 64448/6561 -212/729 0 0 0
               9017/3168 -355/33 46732/5247 49/176 -5103/18656 0 0
                35/384 0 500/1113 125/192 -2187/6784 11/84 0]
            B_5=[35/384 0 500/1113 125/192 -2187/6784 11/84 0]
            B_4=[5179/57600 0 7571/16695 393/640 -92097/339200 187/2100 1/40]
            uis=[]
            for i in 1:len
                push!(uis,Float64[u0s[i]])
            end
            ti=Float64[]
            push!(ti,a)
            while ti[end]<b
                ks=zeros(7,len) 
                for i in 1:7
                    t_step=ti[end]+C[i]*h
                    u_step=zeros(len)
                    for l in 1:len
                        u_step[l]=sum(A[i,:].*ks[:,l]')+uis[l][end]                    
                    end
                    for l in 1:len
                        f=fs[l]
                        ks[i,l]=h*f(t_step,u_step...)
                    end
                end
                
                u_4=zeros(len)
                u_5=zeros(len)
                for l in 1:len
                    u_4[l]=sum(B_4.*ks[:,l]')+uis[l][end]
                    u_5[l]=sum(B_5.*ks[:,l]')+uis[l][end]
                end
                errs=abs.(u_4-u_5)
                err=maximum(errs)
                if err<tol
                    for l in 1:len
                        push!(uis[l],u_5[l])
                    end
                    push!(ti,ti[end]+h) 
                end
                
                h=0.8*h*((tol/err)^(1/5))
                h=minimum([h,h_max])
                if ti[end]+h==ti[end]
                    println("Error during execution: step zero.")
                    return uis,ti
                end 
                if ti[end]+h>b
                    h=b-ti[end]
                end
            end
            return uis,ti
        end

    end

    #------------------------------------------------------------------------------------------------------------------

end