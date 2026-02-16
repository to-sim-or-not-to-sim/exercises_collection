include("../MyLib.jl")
using PyPlot,LinearAlgebra
import .MyLib.MyInterpolation as my
import .MyLib.MyLinearAlgebra as myLA

# USEFUL FUNCTIONS --------------------------------------------------

function mean(x::Vector{T}) where T<:Number
    sum=0
    len=length(x)
    for i in 1:len
        sum+=x[i]
    end
    return sum/len
end

function cheby_poly(x::Number,n::Int64)
    if n==0
        return 1
    elseif n==1
        return x
    elseif n>=2
        return 2*x*cheby_poly(x,n-1)-cheby_poly(x,n-2)
    else
        println("Error, n must be a positive integer.")
    end
    return
end

#--------------------------------------------------------------------

# EX 3.1 ------------------------------------------------------------

function ex1()
    colors=["darkolivegreen",
            "green",
            "springgreen",
            "aquamarine",
            "aqua",
            "deepskyblue",
            "dodgerblue",
            "blue",
            "navy",
            "darkmagenta"]
    for n in 4:4:40
        h=2/n
        i=collect(1:n+1)
        x=@. -1+((i-1)*h)
        V=myLA.Vandermonde_matrix(x,n+1)
        println("cond(V,1)=",cond(V,1),"\nmy cond(V)=",myLA.one_norm_matrix_condition_number(V))
        println("(probably too ill conditioned to work with LU-decomposition)\n")
        println("cond(V)=",cond(V))
        title("n=$n",fontsize=20)
        plot([0,0],[-1,1],linestyle="--",color="gainsboro")
        println("$n")
        grid()
        value_of_interest=zeros(10)
        for k in 10:10:100
            y=@. cos(k*x)
            z=myLA.matrix_solver_Qless(V,y)
            function g(x)
                sum=0
                len=length(z)
                for i in 1:len
                    sum+=z[i]*x^(i-1)        
                end
                return sum
            end
            value_of_interest[Int(k/10)]=g(0)
            println("k=$k, f(0)=",g(0))
            x_plot=collect(-1:0.01:1)
            plot(x_plot,g.(x_plot),color=colors[Int(k/10)],label="k=$k")
        end
        println("min=",minimum(value_of_interest),"\nmax=",maximum(value_of_interest),"\nmean=",mean(value_of_interest))
        println("\n")
        xlim(-1,2)
        ylim(-1.25,1.25)
        legend(loc="upper right",fontsize=15)
        savefig("n$n.png",dpi=300)
        show()
    end
end

#--------------------------------------------------------------------

# EX 3.2 ------------------------------------------------------------

function ex2()
    function ex(f,a,b,n_a,n_b,label,name;step=1)
        x_0=my.linspace(a,b,100)
        plot(x_0,f.(x_0),color="indigo",linewidth=3,label=label)
        colors=["violet","magenta","crimson"]
        k=1
        for n in n_a:step:n_b
            x=my.linspace(a,b,n)
            y=f.(x)
            y_0=zeros(100)
            for j in 1:100
                y_0[j]=my.equi_interpolation(x_0[j],a,b,n,f)
            end
            plot(x_0,y_0,color=colors[k],label="n = $n")
            k+=1
        end
        grid()
        legend(fontsize=15)
        xlabel("x",fontsize=15)
        ylabel("y",fontsize=15)
        savefig("$name.png",dpi=200)
        show()
    end

    f(x)=log(x)
    g(x)=tanh(x)
    h(x)=cosh(x)
    i(x)=abs(x)
    labels=["\$\\log(x)\$","\$\\tanh(x)\$","\$\\cosh(x)\$","\$|x|\$"]
    names=["log","tanh","cosh","abs"]
    ex(f,1,10,2,4,labels[1],names[1])
    ex(g,-3,2,2,4,labels[2],names[2])
    ex(h,-1,3,2,4,labels[3],names[3])
    ex(i,-2,1,3,7,labels[4],names[4],step=2)
end

#--------------------------------------------------------------------

# EX 3.3.1 ----------------------------------------------------------

function ex3(points::Int64)
    a(x)=1/(25*x^2+1)
    b(x)=tanh(5*x+2)
    c(x)=cosh(sin(x))
    d(x)=sin(cosh(x))
    functions=[a,b,c,d]
    x=my.linspace(-1,1,points)
    len=length(x)
    y_fit=zeros(Float64,len)
    y_fit2=zeros(Float64,len)
    N=collect(4:4:60)
    err=zeros(Float64,length(N))
    err_eq=zeros(Float64,length(N))
    for j in 1:4
        f=functions[j]
        y=f.(x)
        k=1
        fig,ax=subplots(ncols=2,figsize=(11,5))
        for n in N
            for i in 1:len
                y_fit[i]=my.chebyshev_interpolation(x[i],f,n)
                y_fit2[i]=my.equi_interpolation(x[i],-1,1,n,f)
            end
            err[k]=my.inf_norm(y,y_fit)
            err_eq[k]=my.inf_norm(y,y_fit2)
            k+=1
        end	

		if j<=2
			A=myLA.Vandermonde_matrix(N,2)                
			z=myLA.linear_leastsquares(A,log.(err))       
			println(z,"  C=",exp(z[1]),"  K=",exp(-z[2]))  
			ax[2].plot(N,exp(z[1])*exp(-z[2]).^(-N),color="black",label="Error estimation")  						
		else
			A=myLA.Vandermonde_matrix(N[1:5],2)                
			z=myLA.linear_leastsquares(A,log.(err[1:5]))       
			println(z,"  C=",exp(z[1]),"  K=",exp(-z[2])) 
			ax[2].plot(N[1:5],exp(z[1])*exp(-z[2]).^(-N[1:5]),color="black",label="Error estimation")  			
		end
        ax[1].grid()
        ax[2].grid()
        ax[1].set_title("Interpolation",fontsize=13.5)
        ax[1].plot(x,y,color="aqua",label="f")
        ax[1].plot(x,y_fit,color="navy",linestyle="--",label="Chebyshev n=60")
        ax[1].plot(x,y_fit2,color="forestgreen",linestyle=":",label="Equally distributed nodes, n=60")
        ax[1].set_xlabel("x",fontsize=13.5)
        ax[1].set_ylabel("f(x)",fontsize=13.5)
        ax[1].legend(fontsize=13.5)
        if maximum(abs.(y_fit2))>2
            ax[1].set_ylim(-2,2)
        end
        ax[2].set_title("Inf-norm",fontsize=13.5)
        ax[2].scatter(N,err_eq,color="forestgreen",label="Equally distributed nodes")
        ax[2].scatter(N,err,color="crimson",label="Chebyshev nodes")
        ax[2].set_yscale("log")
        ax[2].set_xlabel("n",fontsize=13.5)
        ax[2].set_ylabel("\$||f-p||_\\infty\$",fontsize=13.5)
        ax[2].legend(fontsize=13.5)
        savefig("chebyshev$j.png",dpi=500)
		show()
    end
end

#--------------------------------------------------------------------

# EX 3.3.2 ----------------------------------------------------------

function ex4()
    f(x)=cosh(sin(x))
    n=40
    x=my.linspace(0,2*pi,100)
    len=length(x)
    y=zeros(Float64,len)
    for i in 1:len
        y[i]=my.chebyshev_interpolation(x[i],f,n,a=0,b=2*pi)
    end
    plot(x,f.(x),color="indigo",linewidth=3,label="\$f(x)\$")
    plot(x,y,linestyle="--",color="violet",label="\$p(x)\$")
    legend(fontsize=15,bbox_to_anchor=(1,1))
    xlabel("x",fontsize=15)
    ylabel("y",fontsize=15)
    tight_layout()
    grid()
    savefig("chebyshev_cosh.png",dpi=500)
    show()
end

#--------------------------------------------------------------------

# EX 3.3.3 ----------------------------------------------------------

function ex5()
    f(x)=(x^2-2*x+2)^(-1)
    n=30
	x=my.linspace(-6,6,5000)
    len=length(x)
    y=zeros(Float64,len)
    for i in 1:len
        y[i]=my.chebyshev_interpolation(x[i],f,n,a=-6,b=6,kind=1)
    end
    fig,ax=subplots(ncols=2,figsize=(9,5))
    ax[1].set_title("Interpolation",fontsize=15)
    ax[1].plot(x,f.(x),color="black",label="f(x)")
	ax[1].plot(x,y,color="navy",label="p(x) over \$[-6,6]\$")
    ax[1].plot(x,my.real_line_interpolation.(x,f,n),color="crimson",label="p(x) over \$\\mathbb{R}\$")
    ax[1].set_xlabel("x",fontsize=15)
    ax[1].grid()
    ax[1].set_ylabel("y",fontsize=15)
	ax[2].set_title("Zoom",fontsize=15)
    ax[2].plot(x,f.(x),color="black")
	ax[2].plot(x,y,color="navy")
    ax[2].plot(x,my.real_line_interpolation.(x,f,n),color="crimson")
    ax[2].set_xlabel("x",fontsize=15)
    ax[2].set_ylabel("y",fontsize=15)
    ax[2].grid()
	ax[2].set_xlim(6,5.5)
    ax[2].set_ylim(0.03,0.05)
    fig.legend(fontsize=15,loc="upper center")
    tight_layout()
    savefig("real_line_interpolation.png",dpi=500)
    show()
    err_real=abs.(my.real_line_interpolation.(x,f,n).-f.(x))
    err_cheby=abs.(y.-f.(x))
    plot(x,err_real,color="deepskyblue",label="\$ |p_{\\mathbb{R}}(x)-f(x)| \$")
    plot(x,err_cheby,color="navy",label="\$ |p_{[-6,6]}(x)-f(x)| \$")
    xlabel("x",fontsize=15)
    ylabel("error",fontsize=15)
    legend(fontsize=15,loc=(0.75,0.75))
    grid()
    tight_layout()
    savefig("real_line_vs_cheby_finite.png")
    show()
end

#--------------------------------------------------------------------

# EX 3.3.4 ----------------------------------------------------------

function ex6()
    x=my.linspace(-1,1,4000)
    len=length(x)
    y=zeros(Float64,len)
    N=collect(10:10:100)
    err=zeros(Float64,length(N))
	for m in 1:2:11
        fig,ax=subplots()
		f(x)=(abs(x))^m
		y_exp=f.(x)    
		plot(x,y_exp,label="\$|x|^{$m}\$")
		k=1
		for n in 10:10:100
			for i in 1:len
				y[i]=my.chebyshev_interpolation(x[i],f,n)
			end
			err[k]=my.inf_norm(y_exp,y)
			k+=1
			plot(x,y,label="$n")
		end
        grid()
		fig.legend(fontsize=15,loc="right")
        tight_layout()
		show()
        fig,ax=subplots()
        if m<9
            A=myLA.Vandermonde_matrix(log.(N),2)                
            z=myLA.linear_leastsquares(A,log.(err))       
            println("\nC=",exp(z[1]),"  m=",-z[2]) 
            plot(N,exp(z[1])*N.^(z[2]),color="black",label="\$Cn^{-m}\$")
        elseif m==9
            A=myLA.Vandermonde_matrix(log.(N[1:5]),2)                
            z=myLA.linear_leastsquares(A,log.(err[1:5]))       
            println("C=",exp(z[1]),"  m=",-z[2]) 
            plot(N[1:5],exp(z[1])*N[1:5].^(z[2]),color="black",label="\$Cn^{-m}\$")
        else
            A=myLA.Vandermonde_matrix(log.(N[1:4]),2)                
            z=myLA.linear_leastsquares(A,log.(err[1:4]))       
            println("C=",exp(z[1]),"  m=",-z[2]) 
            plot(N[1:4],exp(z[1])*N[1:4].^(z[2]),color="black",label="\$Cn^{-m}\$")
        end
        println("\nAnother way for m=",@. -log(err[2:end]/err[1:end-1])/log(N[2:end]/N[1:end-1]))
		scatter(N,err,label="\$||f-p||_\\infty \$")
        plot(N,1 ./N.^m,color="orangered",linestyle="--",label="\$o(n^{-m})\$")
        fig.legend(fontsize=15,loc="right")
		xscale("log")
		yscale("log")
        grid()
        tight_layout()
        savefig("fit$m.png")
		show()
	end
end

#--------------------------------------------------------------------

# OPTIONAL 1 --------------------------------------------------------

function ex1_optional()
    cond_V_std=zeros(17)
    cond_V_cos=zeros(17)
    cond_P_cheby=zeros(17)
    cond_P_cheby_cos=zeros(17)
    ns=collect(4:20)
    for n in ns
        step=2/n
        x=collect(-1:step:1)
        i=collect(1:n+1)
        x_cos=@. -cos(((2*i-1)/(n+1))*pi/2)
        
        V_std=myLA.Vandermonde_matrix(x,n+1)
        V_cos=myLA.Vandermonde_matrix(x_cos,n+1)
        P_cheby=zeros(n+1,n+1)
        P_cheby_cos=zeros(n+1,n+1)
        for j in 1:n+1
            P_cheby[:,j]=cheby_poly.(x,j-1)
            P_cheby_cos[:,j]=cheby_poly.(x_cos,j-1)
        end
        cond_V_std[n-3]=myLA.one_norm_matrix_condition_number(V_std)
        cond_V_cos[n-3]=myLA.one_norm_matrix_condition_number(V_cos)
        cond_P_cheby[n-3]=myLA.one_norm_matrix_condition_number(P_cheby)
        cond_P_cheby_cos[n-3]=myLA.one_norm_matrix_condition_number(P_cheby_cos)

        println(n," & ",
                round(cond_V_std[n-3],sigdigits=3)," & ",
                round(cond_V_cos[n-3],sigdigits=3)," & ",
                round(cond_P_cheby[n-3],sigdigits=3)," & ",
                round(cond_P_cheby_cos[n-3],sigdigits=3)," \\\\")
    end
    fig,ax=subplots(ncols=2,figsize=(10,5))
    ax[1].plot(ns,cond_V_std,color="royalblue",label="\$V_{std}\$")
    ax[1].plot(ns,cond_V_cos,color="lime",label="\$V_{cos}\$")
    ax[1].plot(ns,cond_P_cheby,color="darkorange",label="\$P_{std}\$")
    ax[1].plot(ns,cond_P_cheby_cos,color="blueviolet",label="\$P_{cos}\$")
    ax[1].grid()
    ax[1].set_xlabel("n",fontsize=15)
    ax[1].set_ylabel("\$\\kappa\$",fontsize=15)
    ax[1].set_yscale("log")
    ax[1].legend(fontsize=15)
    ax[2].plot(ns,cond_V_std,color="royalblue",label="\$V_{std}\$")
    ax[2].plot(ns,cond_V_cos,color="lime",label="\$V_{cos}\$")
    ax[2].plot(ns,cond_P_cheby,color="darkorange",label="\$P_{std}\$")
    ax[2].plot(ns,cond_P_cheby_cos,color="blueviolet",label="\$P_{cos}\$")
    ax[2].set_xlabel("n",fontsize=15)
    ax[2].set_ylabel("\$\\kappa\$",fontsize=15)
    ax[2].set_ylim(0,60)
    ax[2].grid()
    ax[2].legend(fontsize=15)
    savefig("cond.png")
    show()
end

#--------------------------------------------------------------------

# OPTIONAL 2 --------------------------------------------------------

function ex2_optional()
    yscale("log")
    colors=["indianred","slateblue","mediumspringgreen"]
    k=1
    for n in 30:30:90
        xs=my.linspace(-1,1,n)
        lambda=my.lambdas(xs)
        plot(xs,abs.(lambda),color=colors[k],label="n=$n")
        k+=1
    end
    legend(fontsize=15,loc=(0.7,0.7))
    grid()
    xlabel("\$x_i\$",fontsize=15)
    ylabel("\$|\\lambda_i|\$",fontsize=15)
    savefig("weights.png")
    show()
end

#--------------------------------------------------------------------

# OPTIONAL 3 --------------------------------------------------------

function ex3_optional()

    function ex(f,a,b,n_a,n_b,label,name;step=1)
        colors=["violet","magenta","crimson"]
        x_0=my.linspace(a,b,100)
        k=1
        for n in n_a:step:n_b
            x=my.linspace(a,b,n)
            y_0_n=zeros(100)
            y_0_l=zeros(100)
            for j in 1:100
                y_0_n[j]=my.newton_interpolation(x_0[j],x,f)
                y_0_l[j]=my.equi_interpolation(x_0[j],a,b,n,f)
            end
            plot(x_0,abs.(y_0_n.-y_0_l),color=colors[k],label="n = $n")
            k+=1
        end
        grid()
        legend(fontsize=15)
        xlabel("x",fontsize=15)
        ylabel("\$ |p_n(x)-p_l(x)| \$",fontsize=15)
        savefig("newton_$name.png",dpi=200)
        show()
    end

    f(x)=log(x)
    g(x)=tanh(x)
    h(x)=cosh(x)
    i(x)=abs(x)
    labels=["\$\\log(x)\$","\$\\tanh(x)\$","\$\\cosh(x)\$","\$|x|\$"]
    names=["log","tanh","cosh","abs"]
    ex(f,1,10,2,4,labels[1],names[1])
    ex(g,-3,2,2,4,labels[2],names[2])
    ex(h,-1,3,2,4,labels[3],names[3])
    ex(i,-2,1,3,7,labels[4],names[4],step=2)
end

#--------------------------------------------------------------------

# OPTIONAL 4 --------------------------------------------------------

function ex4_optional()
    a(x)=exp(sin(2*pi*x))
    b(x)=log(2+sin(3*pi*x))
    c(x)=(cos(pi*(x-0.2)))^12
    fs=[a,b,c]
    x=my.linspace(0,2,1000)
    Ns=collect(2:30)
    for k in 1:3
        f=fs[k]
        y=f.(x)
        fig,ax=subplots(ncols=2,figsize=(12,5))
        err_even=zeros(Float64,div(Ns[end]-Ns[1],2)+1)
        err_odd=zeros(Float64,div(Ns[end]-Ns[1],2))
        for N in Ns
            v=my.inf_norm(y,my.trigonometric_interpolation.(x,f,N,a=0,b=2))
            if N%2==0
                err_even[Int(N/2)]=v
            else
                err_odd[Int((N+1)/2)-1]=v
            end
        end
        ax[1].set_title("Trigonometric interpolation",fontsize=15)
        ax[2].set_title("Inf-norm over N",fontsize=15)
        ax[1].plot(x,y,color="indigo",label="f(x)")
        ax[1].plot(x,my.trigonometric_interpolation.(x,f,3,a=0,b=2),color="violet",label="Fit with N=3")
        ax[1].plot(x,my.trigonometric_interpolation.(x,f,6,a=0,b=2),color="magenta",label="Fit with N=6")
        ax[1].plot(x,my.trigonometric_interpolation.(x,f,9,a=0,b=2),color="crimson",label="Fit with N=9")
        ax[2].plot(collect(Ns[1]:2:Ns[end]),err_even,color="navy",marker="o",label="N even")
        ax[2].plot(collect(Ns[1]+1:2:Ns[end]),err_odd,color="darkorange",marker="o",label="N odd")
        ax[1].legend(fontsize=15,loc=(0.9,0.65))
        if k==3
            ax[2].legend(fontsize=15,loc=(0.25,0.25))
        else
            ax[2].legend(fontsize=15,loc=(0.5,0.75))
        end
        ax[1].set_xlabel("x",fontsize=15)
        ax[1].set_ylabel("y",fontsize=15)
        ax[2].set_xlabel("N",fontsize=15)
        ax[2].set_ylabel("inf-norm",fontsize=15)
        ax[1].grid()
        ax[2].grid()
        ax[2].set_yscale("log")
        tight_layout()
        savefig("trigon$k.png",dpi=300)
        show()
    end
end

#--------------------------------------------------------------------

# OPTIONAL 5 --------------------------------------------------------

function ex5_optional()

    function width(y_interp::Vector{Float64})
        len=length(y_interp)
        y=y_interp.-1
        k_start=0
        k_end=0
        is_negative=true
        if y[1]>0
            is_negative=false
        end
        is_negative_previous=is_negative
        for i in 1:len
            if y[i]>0
                is_negative=false
            else
                is_negative=true
            end

            if is_negative_previous!=is_negative
                if k_start==0
                    k_start=i
                else
                    k_end=i
                    return k_end,k_start
                end
            end

            is_negative_previous=is_negative
        end
        return k_end,k_start
    end

    f(x)=sign(x+eps())
    a=-0.05
    b=0.15
    len=20000
    step=(b-a)/len
    x=(a:step:b)
    #n=[30,80,180]
    n=collect(10:10:300)
    Ns=2n.+1
    colors=["violet","magenta","crimson"]
    plot(x,f.(x),color="indigo",label="f(x)")
    k=0
    w=zeros(length(Ns))
    h=zeros(length(Ns))
    for j in 1:length(Ns)
        N=Ns[j]
        y=my.trigonometric_interpolation.(x,f,N,a=a,b=b)
        if N==61 || N==161 || N==361
            k+=1
            plot(x,y,color=colors[k],label="N = $N")
        end
        k_end,k_start=width(y)
        w[j]=(k_end-k_start)*step
        h[j]=maximum(y[k_start:k_end])
    end
    grid()
    xlabel("x",fontsize=15)
    ylabel("y",fontsize=15)
    ylim(-1.35,1.35)
    legend(fontsize=15)
    tight_layout()
    savefig("sign1.png")
    show()
    
    fig,ax=subplots(ncols=2,figsize=(10,5))
    ax[1].plot(Ns,w,marker="o",color="navy",label="width")
    ax[2].plot(Ns,h.-1,marker="o",color="navy",label="height")
    ax[1].grid()
    ax[1].set_xlabel("N",fontsize=15)
    ax[1].set_ylabel("width",fontsize=15)
    ax[2].grid()
    ax[2].set_xlabel("N",fontsize=15)
    ax[2].set_ylabel("height",fontsize=15)
    ax[1].legend(fontsize=15)
    ax[2].legend(fontsize=15)
    tight_layout()
    savefig("sign2.png")
    show()

end

#--------------------------------------------------------------------

# EXECUTION ---------------------------------------------------------

function main()
    println("--- MAIN ---")
    println("Exercise 3.1") #LINES 34-82
    ex1()
    println("-"^30)
    println("\nExercise 3.2") #LINES 88-122
    ex2()
    println("-"^30)
    println("\nExercise 3.3.1") #LINES 128-189
    ex3(4000)
    println("-"^30)
    println("\nExercise 3.3.2") #LINES 195-213
    ex4()
    println("-"^30)
    println("\nExercise 3.3.3") #LINES 219-260
    ex5()
    println("-"^30)
    println("\nExercise 3.3.4") #LINES 266-318
    ex6()
    println("-"^30)
end

function optional()
    println("--- OPTIONAL ---")
    println("Exercise C.1") #LINES 324-376
    ex1_optional()
    println("-"^30)
    println("\nExercise C.2") #LINES 382-398
    ex2_optional()
    println("-"^30)
    println("\nExercise C.3") #LINES 404-439
    ex3_optional()
    println("-"^30)
    println("\nExercise C.4.1") #LINES 445-491
    ex4_optional()
    println("-"^30)
    println("\nExercise C.4.2") #LINES 497-579
    ex5_optional()
    println("-"^30)
end

#=
 IT IS POSSIBLE TO ANALYZE SINGLE EXERCISES BY
 REPLACING THE NEXT FEW LINE WITH THE FUNCTION
 RELATED TO THE EXERCISE YOU WANT TO ANALYZE.
=#

main() 
optional()