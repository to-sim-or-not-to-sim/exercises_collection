include("../MyLib.jl")
using PyPlot,PrettyTables
import .MyLib.MyRootFinding as my
import .MyLib.MyLinearAlgebra as LA

# EX 4.1 ------------------------------------------------------------

function ex1()
    f(x)=x^3-7*x^2+14*x-6

    m,all_ms=my.bisection(f,0,1,tol=eps())
    rounded_m=round(m,digits=8)
    println(m,all_ms)
    x=collect(0:0.01:1)
    fig,ax=subplots(ncols=2,figsize=(11,5))
    ax[1].plot(x,zeros(length(x)),linestyle="--",color="purple")
    ax[1].plot(m*ones(2),[f(0),f(1)],linestyle="--",color="purple")
    ax[1].plot(x,f.(x),color="violet",label="Function")
    ax[1].scatter(all_ms,f.(all_ms),color="indigo",s=10,label="Bisection points")
    ax[1].scatter(m,f(m),color="purple",s=50,label="\$x\\sim\$$rounded_m")
    ax[1].grid()
    ax[1].legend(fontsize=15)
    n=collect(1:length(all_ms))
    y=abs.(m.-all_ms)
    ax[2].plot(n,y,color="darkorange",label="Convergence")
    ax[2].scatter(n,y,s=20,color="darkorange")
    ax[2].plot(n,2.0.^(-n),color="navy",linestyle="--",label="Theoretical convergence")
    ax[2].set_yscale("log")
    ax[2].grid()
    ax[2].legend(fontsize=15)
    savefig("bisection1.png")
    show()

    start=2
    y_fit=zeros(length(all_ms)-1-start)
    x_fit=zeros(length(all_ms)-1-start)
    for i in start:length(all_ms)-2
        y_fit[i+1-start]=-log10(abs(all_ms[i+1]-m))
        x_fit[i+1-start]=-log10(abs(all_ms[i]-m))
    end
    scatter(x_fit,y_fit,label="\$(a_n,a_{n+1})\$")
    grid()
    z=LA.poly_linear_fit(x_fit,y_fit,2,plot=true,label="Interpolation")
    legend(fontsize=15)
    println("C= ",10^(-z[1]),"\nq= ",z[2])
    est_C=abs.(all_ms[3:end].-all_ms[2:end-1])./abs.(all_ms[2:end-1].-all_ms[1:end-2])
    println("another way: q=1 => C= ", est_C)
    savefig("bisection2.png")
    show()
    new_ms=zeros(12)
    for i in 1:6
        a=i*0.05
        new_m,new_all_m=my.bisection(f,a,1)
        new_ms[i]=new_m
    end
    for i in 1:6
        b=1-i*0.05
        new_m,new_all_m=my.bisection(f,0,b)
        new_ms[i+6]=new_m
    end
    println(new_ms)
end

#--------------------------------------------------------------------

# EX 4.2.1 ----------------------------------------------------------

function ex2()
    function ex(f::Function,df::Function,exact_q::Number,a::Number,b::Number,name::String,start::String)
        step=(b-a)/99
        x=collect(a:step:b)
        fig,ax=subplots(ncols=2,figsize=(11,5))
        ax[1].plot(x,f.(x),color="violet",label="f(x)")
        x_0,xs=my.newton(f,df,a,b,bracketing=false,set_start=start)
        println("x$name= ",x_0)
        ax[1].scatter(xs,f.(xs),s=10,color="indigo",label="x during Newton")
        ax[1].scatter(x_0,f(x_0),color="purple",label="\$\\bar x: f(\\bar x)=0\$")
        ax[1].grid()
        ax[1].legend(fontsize=15)
        len=length(xs)-1
        fracs=zeros(Float64,len-2)
        x_fit=zeros(Float64,len-2)
        y_fit=zeros(Float64,len-2)
        n=collect(1:len-2)
        an=-log10.(abs.(xs.-x_0))
        for i in n
            fracs[i]=an[i+1]/an[i]
            y_fit[i]=an[i+1]
            x_fit[i]=an[i]
        end
        ax[2].plot(n,fracs,color="navy",label="\$\\dfrac{a_{n+1}}{a_n}\$")
        ax[2].scatter(n,fracs,color="navy")
        ax[2].grid()
        ax[2].legend(fontsize=15)
        savefig("newton2_$name$start.png",dpi=500)
        show()
        q=y_fit[end]/x_fit[end]
        println("q= ",q)
        println("C= ",@. 10^(exact_q*x_fit-y_fit))
    end

    f1(x)=x^2-exp(-x)
    df1(x)=2*x+exp(-x)
    f2(x)=2*x-tan(x)
    df2(x)=2-(1/(cos(x))^2)
    f3(x)=exp(x+1)-2-x
    df3(x)=exp(x+1)-1
    as=[-2,-0.2,-2]
    bs=[2,1.4,2]
    fs=[f1,f2,f3]
    dfs=[df1,df2,df3]
    exact_q=[2,3,1]
    names=["a","b","c"]

    for n in 1:3
        a=as[n] 
        b=bs[n]
        f=fs[n]
        df=dfs[n]
        q=exact_q[n]
        name=names[n]
        if n==2
            ex(f,df,3,a,b,name,"a")
            ex(f,df,2,a,b,name,"b")
        else
            ex(f,df,q,a,b,name,"a")
            ex(f,df,q,a,b,name,"b")
        end 
    end
end

#--------------------------------------------------------------------

# EX 4.2.2 ----------------------------------------------------------

function ex3(;bracketing::Bool=false)
    f(x)=1/(x^2)-sin(x)
    df(x)=-2/(x^3)-cos(x)
    x=collect(0.5:0.025:10)
    x_start=collect(1:7)
    x0s=zeros(length(x_start))
    plot(x,f.(x),color="indigo")
    axhline(y=0,xmin=0,xmax=1,color="violet",linestyle="--")
    vlines(x=x_start,ymin=-1,ymax=3.4,color="plum",linestyle="--")
    scatter(x_start,f.(x_start),color="indigo")
    xlabel("x",fontsize=15)
    ylabel("y",fontsize=15)
    grid()
    savefig("newton3.png")
    show()
    for a in x_start
        x0,xs=my.newton(f,df,a,10,bracketing=bracketing,set_start="a")
        x0s[a]=x0
        println(x0)
    end
    data=hcat(x_start,x0s)
    column_labels=[["Starting x","x: f(x)=0"]]
    style=TextTableStyle(first_line_column_label=crayon"green bold")
    table_format=TextTableFormat(borders=text_table_borders__unicode_rounded)
    highlighter=TextHighlighter(
           (data, i, j) -> (j == 2),
           crayon"yellow bold"
       )
    pretty_table(
        data;
        column_labels=column_labels,
        style=style,
        highlighters=[highlighter],
        table_format=table_format
    )
end

#--------------------------------------------------------------------

# EX 4.2.3 ----------------------------------------------------------

function ex4()
    function equation(epsilon,t,tau)
        f(psi)=psi-epsilon*sin(psi)-(2*pi*t)/tau
        df(psi)=1-epsilon*cos(psi)
        return f,df
    end
    tau=1.7610
    epsilon=0.2230
    step=tau/99
    t=collect(0:step:tau)
    function theta(t)
        f,df=equation(epsilon,t,tau)
        psi=my.newton(f,df,0,tau)[1]
        return mod(2*atan(sqrt((1+epsilon)/(1-epsilon))*tan(psi/2)),2*pi)
    end
    title("True anomaly of Eros",fontsize=15)
    plot(t,theta.(t),color="navy")
    xlabel("t [years]",fontsize=15)
    ylabel("\$\\theta\$ [rad]",fontsize=15)
    grid()
    savefig("eros.png")
    show()
end

#--------------------------------------------------------------------

# EX 4.3 ------------------------------------------------------------

function ex5()

    function ex(f::Function,a::Number,b::Number,q_exp::Number,name::String)
        step=(b-a)/99
        x=collect(a:step:b)
        fig,ax=subplots(ncols=2,figsize=(11,5))
        ax[1].plot(x,f.(x),color="violet",label="\$f_$name(x)\$")
        x_0,xs=my.secant(f,a,b,bracketing=false)
        ax[1].scatter(xs,f.(xs),s=10,color="indigo",label="x during secant")
        ax[1].scatter(x_0,f(x_0),color="purple",label="x: f(x)=0")
        ax[1].axhline(y=0,xmin=0,xmax=1,color="purple",linestyle="--")
        ax[1].grid()
        ax[1].legend(fontsize=15)
        len=length(xs)-1
        fracs=zeros(Float64,len-2)
        x_fit=zeros(Float64,len-2)
        y_fit=zeros(Float64,len-2)
        n=collect(1:len-2)
        an=log10.(abs.(xs.-x_0))
        println("x= ",x_0)
        for i in n
            fracs[i]=an[i+1]/an[i]
            y_fit[i]=an[i+1]
            x_fit[i]=an[i]
        end
        ax[1].set_xlim(a,b)
        ax[1].set_ylim(minimum(f.(x))-0.5,maximum(f.(x))+0.5)
        ax[2].plot(n,fracs,marker="o",color="darkorange",label="\$\\dfrac{a_{n+1}}{a_n}\$")
        ax[2].set_ylim(0,5)
        q_rounded=round(q_exp,digits=3)
        ax[2].plot(n,q_exp*ones(length(n)),color="navy",linestyle="--",label="y=$q_rounded")
        ax[2].grid()
        ax[2].legend(fontsize=15)
        savefig("secant_$name.png")
        show()

        println("q= ",fracs)
        println("q_exp= ",q_exp)
        println("c= ",@. abs(y_fit)/(abs(x_fit)^q_exp))
    end
        
    f1(x)=x^2-exp(-x)
    f2(x)=2*x-tan(x)
    f3(x)=exp(x+1)-2-x
    as=[[-2],[-0.2,0.8],[-2]]
    bs=[[2],[0.8,1.4],[2]]
    names=[["a"],["b^{(1)}","b^{(2)}"],["c"]]
    fs=[f1,f2,f3]
    r=(1+sqrt(5))/2
    q_exp=[r,r,1]

    for n in 1:3
        len=length(as[n])
        for i in 1:len
            f=fs[n]
            a=as[n][i] 
            b=bs[n][i]
            q=q_exp[n]
            name=names[n][i]
            ex(f,a,b,q,name)
        end
    end
end

#--------------------------------------------------------------------

# EX 4.4 ------------------------------------------------------------

function ex6()

    function ex(f::Function,x1::Number,x2::Number,x3::Number,name::String)
        a=x1
        b=x3
        step=(b-a)/99
        x=collect(a:step:b)
        fig,ax=subplots(ncols=2,figsize=(11,5))
        ax[1].plot(x,f.(x),color="violet",label="\$f_$name(x)\$")
        x_0,xs=my.interpolation_method(f,x1,x2,x3)
        ax[1].scatter(xs,f.(xs),s=10,color="indigo",label="x during interpolation")
        ax[1].scatter(x_0,f(x_0),color="purple",label="x: f(x)=0")
        ax[1].axhline(y=0,xmin=0,xmax=1,color="purple",linestyle="--")
        ax[1].grid()
        ax[1].legend(fontsize=15)
        len=length(xs)-1
        fracs=zeros(Float64,len-2)
        x_fit=zeros(Float64,len-2)
        y_fit=zeros(Float64,len-2)
        n=collect(1:len-2)
        an=-log10.(abs.(xs.-x_0))
        println("x= ",x_0)
        for i in n
            fracs[i]=an[i+1]/an[i]
            y_fit[i]=an[i+1]
            x_fit[i]=an[i]
        end
        ax[1].set_xlim(a,b)
        ax[1].set_ylim(minimum(f.(x))-0.5,maximum(f.(x))+0.5)
        ax[2].plot(n,fracs,marker="o",color="darkorange",label="\$\\dfrac{a_{n+1}}{a_n}\$")
        ax[2].set_ylim(0,5)
        ax[2].grid()
        ax[2].legend(fontsize=15)
        savefig("inverse_quad_interpol$name.png")
        show()

        println("q= ",fracs[end])
    end

    f1(x)=x^2-exp(-x)
    f2(x)=2*x-tan(x)
    f3(x)=exp(x+1)-2-x
    x1s=[[-2],[-0.2,0.8],[-2]]
    x2s=[[-0.5],[0.2,1],[0]]
    x3s=[[2],[0.8,1.4],[2]]
    fs=[f1,f2,f3]
    names=[["a"],["b^{(1)}","b^{(2)}"],["c"]]

    for n in 1:3
        len=length(x1s[n])
        for i in 1:len
            f=fs[n]
            x1=x1s[n][i] 
            x2=x2s[n][i] 
            x3=x3s[n][i]
            name=names[n][i]
            ex(f,x1,x2,x3,name)
        end
    end
end

#--------------------------------------------------------------------

# OPTIONAL 1 --------------------------------------------------------

function ex1_optional()
    f(x)=x*exp(x)
    df(x)=(x-1)*exp(x)
    n=200
    
    step_std=2/(n-1)
    x_std=collect(0:step_std:2)
    plot(x_std,f.(x_std),color="cyan",label="\$ f(x) \$")
    y_0=f(0)
    y_1=f(2)
    step=(y_1-y_0)/(n-1)
    ys=collect(y_0:step:y_1)
    xs=zeros(length(ys))
    i=1
    for y in ys
        xs[i]=my.newton(f,df,0,2,k=y)[1]
        i+=1
    end
    plot([0,maximum(ys)],[0,maximum(ys)],linestyle="--",color="gray")
    plot(ys,xs,color="magenta",label="\$ f^{-1}(x) \$")
    legend(fontsize=15)
    xlabel("x",fontsize=15)
    ylabel("y",fontsize=15)
    grid()
    savefig("inverse.png")
    show()
end

#--------------------------------------------------------------------

# OPTIONAL 2 --------------------------------------------------------

function ex2_optional()
    function ex(f::Function,df::Function,d2f::Function,exact_q::Number,a::Number,b::Number,name::String,start::String)
        step=(b-a)/99
        x=collect(a:step:b)
        fig,ax=subplots(ncols=2,figsize=(11,5))
        ax[1].plot(x,f.(x),color="violet",label="f(x)")
        x_0,xs=my.multiple_newton(f,df,d2f,a,b,bracketing=false,set_start=start)
        println("x$name= ",x_0)
        ax[1].scatter(xs,f.(xs),s=10,color="indigo",label="x during Newton")
        ax[1].scatter(x_0,f(x_0),color="purple",label="\$\\bar x: f(\\bar x)=0\$")
        ax[1].set_xlim(a,b)
        ax[1].set_ylim(minimum(f.(x))-0.1,maximum(f.(x))+0.1)
        ax[1].grid()
        ax[1].legend(fontsize=15)
        len=length(xs)-1
        fracs=zeros(Float64,len-2)
        x_fit=zeros(Float64,len-2)
        y_fit=zeros(Float64,len-2)
        n=collect(1:len-2)
        an=-log10.(abs.(xs.-x_0))
        for i in n
            fracs[i]=an[i+1]/an[i]
            y_fit[i]=an[i+1]
            x_fit[i]=an[i]
        end
        ax[2].plot(n,fracs,color="navy",label="\$\\dfrac{a_{n+1}}{a_n}\$")
        ax[2].scatter(n,fracs,color="navy")
        ax[2].set_ylim(0,5)
        ax[2].grid()
        ax[2].legend(fontsize=15)
        savefig("newton2_$name$start.png",dpi=500)
        show()
    end

    f1(x)=x^2-exp(-x)
    df1(x)=2*x+exp(-x)
    d2f1(x)=2-exp(-x)
    f2(x)=2*x-tan(x)
    df2(x)=2-(1/(cos(x))^2)
    d2f2(x)=-2*tan(x)/((cos(x))^2)
    f3(x)=exp(x+1)-2-x
    df3(x)=exp(x+1)-1
    d2f3(x)=exp(x+1)
    as=[-2,-0.2,-2]
    bs=[2,1.3,2]
    fs=[f1,f2,f3]
    dfs=[df1,df2,df3]
    d2fs=[d2f1,d2f2,d2f3]
    exact_q=[2,3,1]
    names=["a","b","c"]

    for n in 1:3
        a=as[n] 
        b=bs[n]
        f=fs[n]
        df=dfs[n]
        d2f=d2fs[n]
        q=exact_q[n]
        name=names[n]
        if n==2
            ex(f,df,d2f,3,a,b,name,"a")
            ex(f,df,d2f,2,a,b,name,"b")
        else
            ex(f,df,d2f,q,a,b,name,"a")
            ex(f,df,d2f,q,a,b,name,"b")
        end 
    end
end

#--------------------------------------------------------------------

# OPTIONAL 3 --------------------------------------------------------

function ex3_optional()
    f(x)=sin(x)
    df(x)=cos(x)
    a=-pi/2
    b=pi/2
    x0=1.16556118521
    x,xs=my.newton(f,df,a,x0,set_start="b",bracketing=false)
    println(x)
    println(xs)
    println(x/pi)
end

#--------------------------------------------------------------------

# EXECUTION ---------------------------------------------------------

function main()
    println("--- MAIN ---")
    println("Exercise 4.1") #LINES 8-62
    ex1()
    println("-"^30)
    println("\nExercise 4.2.1") #LINES 68-130
    ex2()
    println("-"^30)
    println("\nExercise 4.2.2") #LINES 136-171
    ex3()
    println("-"^30)
    println("\nExercise 4.2.3") #LINES 177-199
    ex4()
    println("-"^30)
    println("\nExercise 4.3") #LINES 205-267
    ex5()
    println("-"^30)
    println("\nExercise 4.4") #LINES 273-332
    ex6()
    println("-"^30)
end

function optional()
    println("--- OPTIONAL ---")
    println("Exercise D.1.1") #LINES 338-364
    ex1_optional()
    println("-"^30)
    println("\nExercise D.1.2") #LINES 370-437
    ex2_optional()
    println("-"^30)
    println("\nExercise D.1.3") #LINES 443-453
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