include("../MyLib.jl")
using LinearAlgebra
using PyPlot
import .MyLib.MyLinearAlgebra as my

#=
 THE CODE CONTAINS ALL THE EXERCISE RELATED TO
 LINEAR SYSTEMS SECTION (BOTH MANDATORY AND 
 OPTIONAL). IT EXECUTES ALL THE EXERCISES. IF 
 YOU WANT TO MODIFY THIS YOU CAN GO TO LINE 561.
=#

# USEFUL FUNCTIONS --------------------------------------------------

function linspace(a::Number,b::Number,len::Int;type::Type=Float64)
        x=zeros(type,len) 
        step=(b-a)/(len-1)
        for i in 1:len
            x[i]=a+(i-1)*step
        end
        return x
    end

#--------------------------------------------------------------------

# EX 2.1 ------------------------------------------------------------

function ex1()
    A_1=[-2 0 0
        1 -1 0
        3 2 1]

    b_1=[-4,2,1]

    println("Problem 1:")
    println(my.lower_triangular_matrix_solver(A_1,b_1))
    println("--------------------")

    A_2=[4 0 0 0
        1 -2 0 0
        -1 4 4 0
        2 -5 5 1]

    b_2=[-4,1,-3,5]

    println("Problem 2:")
    println(my.lower_triangular_matrix_solver(A_2,b_2))
    println("--------------------")

    A_3=[3 1 0
        0 -1 -2
        0 0 3]

    b_3=[1,1,6]

    println("Problem 3:")
    println(my.upper_triangular_matrix_solver(A_3,b_3))
    println("--------------------")

    A_4=[3 1 0 6
        0 -1 -2 7
        0 0 3 4
        0 0 0 5]

    b_4=[4,1,1,5]

    println("Problem 4:")
    println(my.upper_triangular_matrix_solver(A_4,b_4))
    println("--------------------")

    b_5=[-4 1 0
        2 1 0
        1 6 0]

    println("Problem 5:")
    display(my.lower_triangular_matrix_solver(A_1,b_5))
    println("--------------------")

    println("Problem 6:")
    display(my.upper_triangular_matrix_solver(A_3,b_5))
    println("--------------------")
end

#--------------------------------------------------------------------

# EX 2.2 ------------------------------------------------------------

function ex2()

    #FIRST PART
    
    A_1=[2 3  4
        4 5 10
        4 8  2]

    A_2=[  6 -2  -4  4
        3 -3  -6  1
        -12  8  21 -8
        -6  0 -10  7]

    A_3=[ 1  4  5 -5
        -1  0 -1 -5
        1  3 -1  2
        1 -1  5 -1]

    println("EX a")
    my.display_LU(A_1)
    println("EX b")
    my.display_LU(A_2)
    println("EX c")
    my.display_LU(A_3)

    #SECOND PART

    T1=my.translation_matrix(3,-1)
    R=my.rotation_matrix(pi/5)
    T2=my.translation_matrix(-3,1)
    A=T1*R*T2
    z=[2
       2
       1]
    b=A*z
    display(b)
    L,U=my.LU_decomposition(A)
    x=my.matrix_solver_LU(A,b)
    display(x)

    #THIRD PART

    A_1=[2 3  4
         4 5 10
         4 8  2]

    A_2=[  6 -2  -4  4
           3 -3  -6  1
         -12  8  21 -8
          -6  0 -10  7]

    A_3=[ 1  4  5 -5
         -1  0 -1 -5
          1  3 -1  2
          1 -1  5 -1]

    println("My det(A_1)=",my.det_LU(A_1),"\nLinearAlgebra det(A_1)=",det(A_1))
    println("My det(A_2)=",my.det_LU(A_2),"\nLinearAlgebra det(A_2)=",det(A_2))
    println("My det(A_3)=",my.det_LU(A_3),"\nLinearAlgebra det(A_3)=",det(A_3))
end

#--------------------------------------------------------------------

# EX 2.3 ------------------------------------------------------------

function ex3()
    A_1=[2 3  4
         4 5 10
         4 8  2]

    A_2=[  6 -2  -4  4
           3 -3  -6  1
         -12  8  21 -8
          -6  0 -10  7]

    A_3=[ 1  4  5 -5
         -1  0 -1 -5
          1  3 -1  2
          1 -1  5 -1]

    println("EX a")
    my.display_PLU(A_1)
    println("EX b")
    my.display_PLU(A_2)
    println("EX c")
    my.display_PLU(A_3)

    T1=my.translation_matrix(3,-1)
    R=my.rotation_matrix(pi/5)
    T2=my.translation_matrix(-3,1)
    A=T1*R*T2
    z=[2
    2
    1]
    b=A*z
    display(b)
    x=my.matrix_solver_PLU(A,b)
    display(x)


    println("My det(A_1)=",my.det_PLU(A_1),"\nLinearAlgebra det(A_1)=",det(A_1))
    println("My det(A_2)=",my.det_PLU(A_2),"\nLinearAlgebra det(A_2)=",det(A_2))
    println("My det(A_3)=",my.det_PLU(A_3),"\nLinearAlgebra det(A_3)=",det(A_3))

end

#--------------------------------------------------------------------

# EX 2.4 ------------------------------------------------------------

function ex4()
    function exercise(A::AbstractMatrix{T} where T<:Number)
        x=[1,1]
        b=A*x
        L,U=my.LU_decomposition(A)
        z=my.lower_triangular_matrix_solver(L,b)
        x_tilde=my.upper_triangular_matrix_solver(U,z)
        println("Delta x")
        display(abs.(x-x_tilde))
        println("A-LU")
        display(A-L*U)
        println("Mine")
        println("k_A=",my.one_norm_matrix_condition_number(A))
        println("k_L=",my.one_norm_matrix_condition_number(L))
        println("k_U=",my.one_norm_matrix_condition_number(U))
        println("LinearAlgebra")
        println("k_A=",cond(A,1))
        println("k_L=",cond(L,1))
        println("k_U=",cond(U,1))
    end

    EPS=[1e-12,1e-20]
    println("-"^30)
    for eps in EPS
        A=[-eps  1
            1 -1]
        exercise(A)
        println("-"^30)
    end
    for eps in EPS
        A=[   1  1
        -eps -1]
        exercise(A)
        println("-"^30)
    end
end

#--------------------------------------------------------------------

# EX 2.5 ------------------------------------------------------------

function ex5()
    #--------------------------------
    println("--- PART A ---")
    function exA()
        A= [2 -1
            0  1
        -2  2]
        b=[1,-5,6]
        x=my.linear_leastsquares(A,b,cholensky=true)
        display(x)
    end
    exA()
    #--------------------------------
    println("--- PART B ---")
    function exB()
        R=[57.59,108.11,149.57,227.84,778.14,1427,2870.3,4499.9]
        tau=[87.99,224.7,365.26,686.98,4332.4,10759,30684,60188]
        b=log.(tau)
        A=zeros(length(R),2)
        A[:,1]=ones(length(R))
        A[:,2]=log.(R)
        x=my.linear_leastsquares(A,b,cholensky=true)
        c=exp(x[1])
        alpha=x[2]
        println(alpha)
        xs=linspace(57,4600,200)
        plot(xs,c.*xs.^(alpha),color="darkslategrey",label="Linear fit")
        scatter(R,tau,color="crimson",label="Data")
        xlabel("R",fontsize=15)
        ylabel("\$\\tau\$",fontsize=15)
        xscale("log")
        yscale("log")
        grid()
        legend(fontsize=15)
        savefig("r-tau.png",dpi=300)
        show()
    end
    exB()
    #--------------------------------
    println("--- PART C ---")
    function g(t)
        return exp(sin(t-1))
    end
    function exC()
        t_i=linspace(0,2*pi,61)
        b=g.(t_i)
        grid()
        x_1=my.poly_linear_fit(t_i,b,7,plot=true,color="crimson",cholensky=true)
        x_2=my.fourier_linear_fit(t_i,b,2,plot=true,color="indigo",cholensky=true)
        display(x_1)
        display(x_2)
        scatter(t_i,b,s=20,color="slategrey",label="Data")
        xlabel("t",fontsize=15)
        ylabel("g(t)",fontsize=15)
        legend(fontsize=15)
        savefig("fit.png",dpi=300)
        show()
    end
    exC()
    #--------------------------------
end

#--------------------------------------------------------------------

# EX 2.6 ------------------------------------------------------------

function ex6()
    A=[ 2 -1
        0  1
       -2  2]

    B=[2 3  4
       4 5 10
       4 8  2]

    C=[2 -1  4
       0  3  5
       7  2 -6
       1 -4  0]

    Q,R=my.QR_decomposition(A)
    println("Q=")
    display(Q)
    println("I=")
	display(Q'*Q)
    println("R=")
    display(R)
    println("0=")
    display(A-Q*R)

    Q,R=my.QR_decomposition(B)
    println("\nQ=")
    display(Q)
    println("I=")
	display(Q'*Q)
    println("R=")
    display(R)
    println("0=")
    display(B-Q*R)

    Q,R=my.QR_decomposition(C)
    println("\nQ=")
    display(Q)
    println("I=")
	display(Q'*Q)
    println("R=")
    display(R)
    println("0=")
    display(C-Q*R)
end

#--------------------------------------------------------------------

# EX 2.7 ------------------------------------------------------------

function ex7()
    function h(t)
        return @. exp(sin(4*t))
    end
    i=collect(0:99)
    t=i./99
    y=h(t)./2006.787453080206
    A=my.Vandermonde_matrix(t,15)
    At=A'
    N=At*A
    z=At*y

	z_LU=my.matrix_solver_PLU(N,z)
	println("LU: ",z_LU)
	
    z_CHOLE=my.matrix_solver_cholesky(N,z)
	println("CHOLESKY: ",z_CHOLE)

    z_QR=my.matrix_solver_QR(A,y)
    println("QR: ",z_QR)

    z_QLESS=my.matrix_solver_Qless(A,y)
    println("Q-LESS: ",z_QLESS)
	
	function f(x,z)
		len=length(x)
		y=zeros(len)
		for i in 1:len
			for j in 1:length(z)
				y[i]+=z[j]*(x[i]^(j-1))
			end
		end
		return y
	end
	fig,ax=subplots(ncols=2,figsize=(10,5))
    ax[1].scatter(t,y,color="darkslategrey",s=20,label="Data")
	ax[1].plot(t,f(t,z_LU),color="darkorange",label="LU")
	ax[1].plot(t,f(t,z_CHOLE),color="hotpink",label="Cholesky")
	ax[1].plot(t,f(t,z_QR),color="rebeccapurple",label="QR")
	ax[1].plot(t,f(t,z_QLESS),color="mediumseagreen",label="Q-less")
	ax[1].grid()
    ax[1].set_xlabel("t",fontsize=15)
    ax[1].set_ylabel("y(t)",fontsize=15)
	ax[1].legend(fontsize=15)
    y_QR=abs.(y.-f(t,z_QR))
    y_QLESS=abs.(y.-f(t,z_QLESS))
    y_LU=abs.(y.-f(t,z_LU))
    ax[2].plot(t,y_LU,color="darkorange",label="LU")
	ax[2].plot(t,y_QR,color="rebeccapurple",label="QR")
	ax[2].plot(t,y_QLESS,color="mediumseagreen",label="Q-less")
	ax[2].grid()
    ax[2].set_xlabel("t",fontsize=15)
    ax[2].set_ylabel("|y(t)-fit(t)|",fontsize=15)
	ax[2].legend(fontsize=15)
    tight_layout()
    savefig("c15.png",dpi=300)
	show()
end

#--------------------------------------------------------------------

# EX 2.8 ------------------------------------------------------------

function ex8()
    A=[5 1 0 0
       1 4 1 0
       0 1 3 1
       0 0 1 2]
    display(A)
    diagA,Q=my.eigenvalues(A)
    display(diagA)
    display(Q)
    for i in 1:4
        display((A*Q[:,i])-diagA[i,i]*Q[:,i])
    end
end

#--------------------------------------------------------------------

# OPTIONAL 1 --------------------------------------------------------

function ex1_optional()

    A_1=[2 3  4
         4 5 10
         4 8  2]

    A_2=[  6 -2  -4  4
           3 -3  -6  1
         -12  8  21 -8
          -6  0 -10  7]

    A_3=[ 1  4  5 -5
         -1  0 -1 -5
          1  3 -1  2
          1 -1  5 -1]

    display(my.inverse_LU(A_1))
    display(inv(A_1))
    display(my.inverse_LU(A_2))
    display(inv(A_2))
    display(my.inverse_LU(A_3))
    display(inv(A_3))

end

#--------------------------------------------------------------------

# OPTIONAL 2 --------------------------------------------------------

function ex2_optional()
    A=[[ 1 0 -1
         0 4  5
        -1 5 10],
        [1 0  1
         0 4  5
         1 5 10],
        [1 0 1
         0 4 5
         1 5 1], 
        [6 2 1 0
         2 6 2 1
         1 2 5 2
         0 1 2 4]]
    len=length(A)
    for i in 1:len
        println(my.is_positive_definite(A[i]))
    end
end

#--------------------------------------------------------------------

# OPTIONAL 3 --------------------------------------------------------

function ex3_optional()
    t=linspace(0,10,21)
    y=tanh.(t)
    A=my.Vandermonde_matrix(t,4)
    x=my.linear_leastsquares(A,y,cholensky=true)
    a=@. (t^2/(1+t^2))
    B=my.Vandermonde_matrix(a,4)
    z=my.linear_leastsquares(B,y,cholensky=true)
    display(x)
    display(z)
    scatter(t,y,color="indigo",label="Data")
    t_fit=collect(0:0.1:10)
    z_fit=@. (t_fit^2/(1+t_fit^2))
    y_t=zeros(length(t_fit))
    y_z=zeros(length(t_fit))
    for i in 1:4
        y_t=@. y_t+x[i]*t_fit^(i-1)
        y_z=@. y_z+z[i]*z_fit^(i-1)
    end
    plot(t_fit,y_t,color="violet",label="Polynomial fit with t")
    plot(t_fit,y_z,color="crimson",label="Polynomial fit with z")
    grid()
    xlabel("t",fontsize=15)
    legend(fontsize=15)
    savefig("poly_fit_tvsz.png",dpi=300)
    show()
end

#--------------------------------------------------------------------

# EXECUTION ---------------------------------------------------------

function main()
    println("--- MAIN ---")
    println("Exercise 2.1") #LINES 28-82
    ex1()
    println("-"^30)
    println("\nExercise 2.2") #LINES 88-147
    ex2()
    println("-"^30)
    println("\nExercise 2.3") #LINES 153-192
    ex3()
    println("-"^30)
    println("\nExercise 2.4") #LINES 198-233
    ex4()
    println("-"^30)
    println("\nExercise 2.5") #LINES 239-299
    ex5()
    println("-"^30)
    println("\nExercise 2.6") #LINES 305-348
    ex6()
    println("-"^30)
    println("\nExercise 2.7") #LINES 354-411
    ex7()
    println("-"^30)
    println("\nExercise 2.8") #LINES 417-429
    ex8()
    println("-"^30)
end

function optional()
    println("--- OPTIONAL ---")
    println("Exercise B.1") #LINES 435-458
    ex1_optional()
    println("-"^30)
    println("\nExercise B.2") #LINES 464-482
    ex2_optional()
    println("-"^30)
    println("\nExercise B.3") #LINES 488-514
    ex3_optional()
    println("-"^30)
end

#=
 IT IS POSSIBLE TO ANALYZE SINGLE EXERCISES BY
 REPLACING THE NEXT FEW LINE WITH THE FUNCTION
 RELATED TO THE EXERCISE YOU WANT TO ANALYZE.
=#

main() 
optional()