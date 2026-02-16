include("../MyLib.jl")
using PyPlot,PrettyTables,SpecialFunctions
import .MyLib.MyIntegration as my
import .MyLib.MyLinearAlgebra as LA
import .MyLib.MyInterpolation as interpol 

# EX 5.1 ------------------------------------------------------------

function ex1()
    k=collect(1:10)
    ms=@. 10*2^k
    fa(x)=x*log(1+x)
    fb(x)=(x^2)*atan(x)
    fc(x)=exp(x)*cos(x)
    fd(x)=atan(sqrt(2+x^2))/((1+x^2)*sqrt(2+x^2))
    fe(x)=sqrt(x)*log(x)
    ff(x)=sqrt(1-x^2)
    fs=[fa,fb,fc,fd,fe,ff]
    expected_values=[0.25,
                     (pi-2+2*log(2))/12,
                     (exp(pi/2)-1)/2,
                     5*pi^2/96,
                     -4/9,
                     pi/4]
    as=[0,0,0,0,eps(),0]
    bs=[1,1,pi/2,1,1,1]
    titles=["\$x\\log(1+x)\$",
            "\$x^2\\arctan(x)\$",
            "\$e^x\\cos(x)\$",
            "\$\\frac{\\arctan\\left(\\sqrt{2+x^2}\\right)}{(1+x^2)\\sqrt{x+x^2}}\$",
            "\$\\sqrt{x}\\log(x)\$",
            "\$\\sqrt{1-x^2}\$"]
    names=["a","b","c","d","e","f"]
    value_trap=0
    A=zeros(10)
    z_trap=0
    for i in 1:6
        t=titles[i]
        name=names[i]
        title(t,fontsize=15)
        f=fs[i]
        a=as[i]
        b=bs[i]
        exp_value=expected_values[i]
        err_trap=zeros(Float64,10)
        for j in 1:10
            m=ms[j]
            value_trap=my.trapezoidal(f,a,b,m)
            err_trap[j]=abs(value_trap-exp_value)
        end
        x=(b-a)./ms
        scatter(x,err_trap,color="darkorange",label="\$|I[f]-T[h]|\$")
        xlabel("h",fontsize=15)
        grid()
        A[:]=x.^2
        At=A'
        N=At*A
        z=At*err_trap
        z_trap=z/N 
        println("trap= ",z_trap)
        x_fit=collect(x[end]:0.0001:x[1])
        plot(x_fit,z_trap[1].*(x_fit.^2),color="navy",label="\$Ch^2\$")
        legend(fontsize=15)
        xscale("log")
        yscale("log")
        savefig("trapezoidal_$name.png")
        show()
    end
end

#--------------------------------------------------------------------

# EX 5.2 ------------------------------------------------------------

function ex2()
    k=collect(1:10)
    ms=@. 10*2^k
    fa(x)=x*log(1+x)
    fb(x)=(x^2)*atan(x)
    fc(x)=exp(x)*cos(x)
    fd(x)=atan(sqrt(2+x^2))/((1+x^2)*sqrt(2+x^2))
    fe(x)=sqrt(x)*log(x)
    ff(x)=sqrt(1-x^2)
    fs=[fa,fb,fc,fd,fe,ff]
    expected_values=[0.25,
                     (pi-2+2*log(2))/12,
                     (exp(pi/2)-1)/2,
                     5*pi^2/96,
                     -4/9,
                     pi/4]
    as=[0,0,0,0,eps(),0]
    bs=[1,1,pi/2,1,1,1]
    titles=["\$x\\log(1+x)\$",
            "\$x^2\\arctan(x)\$",
            "\$e^x\\cos(x)\$",
            "\$\\frac{\\arctan\\left(\\sqrt{2+x^2}\\right)}{(1+x^2)\\sqrt{x+x^2}}\$",
            "\$\\sqrt{x}\\log(x)\$",
            "\$\\sqrt{1-x^2}\$"]
    names=["a","b","c","d","e","f"]
    value_simp=0
    A=zeros(10)
    z_simp=0
    for i in 1:6
        t=titles[i]
        name=names[i]
        title(t,fontsize=15)
        f=fs[i]
        a=as[i]
        b=bs[i]
        exp_value=expected_values[i]
        err_simp=zeros(Float64,10)
        for j in 1:10
            m=ms[j]
            value_simp=my.simpson(f,a,b,m)
            err_simp[j]=abs(value_simp-exp_value)
        end
        x=(b-a)./ms
        scatter(x,err_simp,color="darkorange",label="\$|I[f]-T[h]|\$")
        xlabel("h",fontsize=15)
        grid()
        A[:]=x.^4
        At=A'
        N=At*A
        z=At*err_simp
        z_simp=z/N 
        println("simp= ",z_simp)
        x_fit=collect(x[end]:0.0001:x[1])
        plot(x_fit,z_simp[1].*(x_fit.^4),color="navy",label="\$Ch^2\$")
        legend(fontsize=15)
        xscale("log")
        yscale("log")
        savefig("simpson_$name.png")
        show()
    end 
end

#--------------------------------------------------------------------

# EX 5.3.1 ----------------------------------------------------------

function ex3()
    all_zeros=my.all_zeros
    all_w=my.all_w

    n=collect(2:5)
    x=Vector{Vector{Float64}}(undef,4)
    w=Vector{Vector{Float64}}(undef,4)

    for i in n
        x[i-1]=all_zeros[i]
        w[i-1]=all_w[i]
    end

    data=hcat(n,x,w)
    column_labels=[["n","x","w"]]
    style=TextTableStyle(first_line_column_label=crayon"yellow bold")
    table_format=TextTableFormat(borders=text_table_borders__unicode_rounded)
    highlighter_x=TextHighlighter(
            (data, i, j) -> (j == 2),
            crayon"magenta bold"
        )
    highlighter_w=TextHighlighter(
            (data, i, j) -> (j == 3),
            crayon"cyan bold"
        )
    pretty_table(
        data;
        column_labels=column_labels,
        style=style,
        highlighters=[highlighter_x,highlighter_w],
        table_format=table_format
    )
end

#--------------------------------------------------------------------

# EX 5.3.2 ----------------------------------------------------------

function ex4()
    fa(x)=exp(-4*x)
    fb(x)=exp(-9x^2)
    fc(x)=sech(x)
    fd(x)=1/(1+9*x^2)
    fe(x)=x^2*sin(8*x)
    fs=[fa,fb,fc,fd,fe]
    as=[-1,-1,-1,-1,pi/2]
    bs=[1,1,1,1,pi]
    exp_values=[sinh(4)/2,sqrt(pi)*erf(3)/3,2*atan(sinh(1)),(2/3)*atan(3),-3*pi^2/32]
    labels=["a","b","c","d","e"]
    colors=["royalblue","crimson","forestgreen","goldenrod","indigo"]
    A=zeros(37,1)
    function_d=zeros(37)
    ns=collect(4:40)
    error=zeros(37)
    for i in 1:5
        f=fs[i]
        a=as[i]
        b=bs[i]
        I=exp_values[i]
        for k in 1:37
            n=ns[k]
            I_n=my.legendre_gauss(f,a,b,n)
            error[k]=abs(I-I_n)
        end
        plot(ns,error,marker="o",color=colors[i],label=labels[i])
        if i==4
            function_d[:]=error[:]
            A[:]=((1/3)+sqrt(1+(1/9))).^(-2 .*ns)
            At=A'
            N=At*A
            b=At*error
            global z=b/N
        end
    end
    xlabel("n",fontsize=15)
    ylabel("\$|I_n-I[f]|\$",fontsize=15)        
    legend(fontsize=15)
    grid()
    yscale("log")
    tight_layout()
    savefig("gauss_legendre.png")
    show()
    scatter(ns,function_d,color="darkorange",marker="o",label="\$|I^{(d)}_n-I[f_d]|\$")
    plot(ns,z.*A,color="navy",label="\$C\\rho^{-2n}\$")
    xlabel("n",fontsize=15)
    ylabel("\$|I_n-I[f]|\$",fontsize=15)
    yscale("log")
    legend(fontsize=15)
    grid()
    tight_layout()
    savefig("bernstein_ellipse_plot.png")
    show()
end

#--------------------------------------------------------------------

# EX 5.4.1 ----------------------------------------------------------

function ex5()
    labels=["a","b","c","d","e"]
    colors=["royalblue","crimson","forestgreen","goldenrod","indigo"]    
    fa(x)=1/(1+x^2+x^4)
    fb(x)=exp(-x^2)*cos(x)
    fc(x)=(1+x^2)^(-2/3)
    fd(x)=1/(1+x^2)
    fe(x)=exp(-x)/sqrt(x)
    fs=[fa,fb,fc,fd,fe]
    as=[my.minus_inf,my.minus_inf,my.minus_inf,0,0]
    exp_values=[pi/sqrt(3),sqrt(pi)*exp(-0.25),sqrt(pi)*gamma(1/6)/gamma(2/3),pi/2,sqrt(pi)]
    ns=collect(4:2:60)
    for i in 1:5
        f=fs[i]
        a=as[i]
        exp_v=exp_values[i]
        l=labels[i]
        c=colors[i]
        len=length(ns)
        error=zeros(len)
        last_v=0
        for j in 1:len
            N=div(ns[j],2)
            v=0
            if a==my.minus_inf
                v=my.double_exponential(f,N)
            else
                if i==5 
                    v=my.double_exponential(f,a,my.plus_inf,N,exp_weigth=true)
                else
                    v=my.double_exponential(f,a,my.plus_inf,N)
                end
            end
            error[j]=abs(v-exp_v)
            last_v=v
        end
        println(l,"--> V=",last_v,"; V_exp=",exp_v)
        plot(ns,error,marker="o",color=c,label=l)
    end
    grid()
    xlabel("n",fontsize=15)
    ylabel("\$|I_n-I[f]|\$",fontsize=15)
    legend(fontsize=15)
    yscale("log")
    savefig("double_exp1.png")
    tight_layout()
    show()
end

#--------------------------------------------------------------------

# EX 5.4.2 ----------------------------------------------------------

function ex6()
    labels=["a","b","c","d","e"]
    colors=["royalblue","crimson","forestgreen","goldenrod","indigo"]    
    fa(x)=sqrt(x)*log(x)
    fb(x)=sqrt(1-x^2)
    fc(x)=(log(x))^2
    fd(x)=log(cos(x))
    fe(x)=sqrt(tan(x))
    fs=[fa,fb,fc,fd,fe]
    a=0
    bs=[1,1,1,pi/2,pi/2]
    exp_values=[-4/9,pi/4,2,-(pi/2)*log(2),pi/sqrt(2)]
    ns=collect(4:2:60)
    fig,ax=subplots(ncols=2,figsize=(10,5))
    for i in 1:5
        f=fs[i]
        b=bs[i]
        exp_v=exp_values[i]
        l=labels[i]
        c=colors[i]
        len=length(ns)
        error_de=zeros(len)
        error_gl=zeros(len)
        last_v_de=0
        last_v_gl=0
        println(l)
        for j in 1:len
            N=div(ns[j],2)
            v_de=my.double_exponential(f,a,b,N)
            v_gl=my.legendre_gauss(f,a,b,ns[j])
            error_de[j]=abs(v_de-exp_v)
            error_gl[j]=abs(v_gl-exp_v)
            last_v_de=v_de
            last_v_gl=v_gl
        end
        println(l,"--> V_double=",last_v_de,"; V_gauss=",last_v_gl,"; V_exp=",exp_v)
        ax[1].plot(ns,error_de,marker="o",color=c,label=l)
        ax[2].plot(ns,error_gl,marker="o",color=c,label=l)
    end
    ax[1].set_title("Double-exponential")
    ax[2].set_title("Gauss-Legendre")
    ax[1].grid()
    ax[2].grid()
    ax[1].legend()
    ax[2].legend()
    ax[1].set_yscale("log")
    ax[2].set_yscale("log")
    savefig("DE_vs_gauss.png")
    show()

end

#--------------------------------------------------------------------

# OPTIONAL 1 --------------------------------------------------------

function ex1_optional()
    points=4000
    all_zeros=my.all_zeros
    a(x)=1/(25*x^2+1)
    b(x)=tanh(5*x+2)
    c(x)=cosh(sin(x))
    d(x)=sin(cosh(x))
    functions=[a,b,c,d]
    x=interpol.linspace(-1,1,points)
    len=length(x)
    y_fit=zeros(Float64,len)
    y_fit2=zeros(Float64,len)
    N=collect(4:4:60)
    err=zeros(Float64,length(N))
    err_leg=zeros(Float64,length(N))
    for j in 1:4
        f=functions[j]
        y=f.(x)
        k=1
        for n in N
            xs_legendre=Float64.(all_zeros[n])
            for i in 1:len
                y_fit[i]=interpol.chebyshev_interpolation(x[i],f,n)
                y_fit2[i]=interpol.interpolation(x[i],xs_legendre,f.(xs_legendre))
            end
            err[k]=interpol.inf_norm(y,y_fit)
            err_leg[k]=interpol.inf_norm(y,y_fit2)
            k+=1
        end	
        grid()
        plot(N,err_leg,marker="o",color="forestgreen",label="Legendre nodes")
        plot(N,err,marker="o",color="crimson",label="Chebyshev nodes")
        yscale("log")
        xlabel("n",fontsize=15)
        ylabel("\$||f-p||_\\infty\$",fontsize=15)
        legend(fontsize=15)
        savefig("chebyshev_vs_legendre$j.png",dpi=500)
        tight_layout()
		show()
    end

end

#--------------------------------------------------------------------

# OPTIONAL 2 --------------------------------------------------------

function ex2_optional()
    fa(x)=exp(-4*x)
    fb(x)=exp(-9x^2)
    fc(x)=sech(x)
    fd(x)=1/(1+9*x^2)
    fe(x)=x^2*sin(8*x)
    fs=[fa,fb,fc,fd,fe]
    as=[-1,-1,-1,-1,pi/2]
    bs=[1,1,1,1,pi]
    exp_values=[sinh(4)/2,sqrt(pi)*erf(3)/3,2*atan(sinh(1)),(2/3)*atan(3),-3*pi^2/32]
    labels=["a","b","c","d","e"]
    colors=["royalblue","crimson","forestgreen","goldenrod","indigo"]
    A=zeros(37,1)
    function_d=zeros(37)
    ns=collect(4:40)
    error=zeros(37)
    for i in 1:5
        f=fs[i]
        a=as[i]
        b=bs[i]
        I=exp_values[i]
        for k in 1:37
            n=ns[k]
            I_n=my.clenshaw_curtis(f,a,b,n)
            error[k]=abs(I-I_n)
        end
        plot(ns,error,marker="o",color=colors[i],label=labels[i])
        if i==4
            function_d[:]=error[:]
            A[:]=((1/3)+sqrt(1+(1/9))).^(-2 .*ns)
            At=A'
            N=At*A
            b=At*error
            global z=b/N
        end
    end
    xlabel("n",fontsize=15)
    ylabel("\$|I_n-I[f]|\$",fontsize=15)        
    legend(fontsize=15)
    grid()
    yscale("log")
    tight_layout()
    savefig("clenshaw_curtis.png")
    show()
end

#--------------------------------------------------------------------

# OPTIONAL 3 --------------------------------------------------------

function ex3_optional()
    fa(x)=exp(-4*x)
    fb(x)=exp(-9x^2)
    fc(x)=sech(x)
    fd(x)=1/(1+9*x^2)
    fe(x)=x^2*sin(8*x)
    fs=[fa,fb,fc,fd,fe]
    as=[-1,-1,-1,-1,pi/2]
    bs=[1,1,1,1,pi]
    exp_values=[sinh(4)/2,sqrt(pi)*erf(3)/3,2*atan(sinh(1)),(2/3)*atan(3),-3*pi^2/32]
    labels=["a","b","c","d","e"]
    colors=["royalblue","crimson","forestgreen","goldenrod","indigo"]
    Ns=collect(2:20)
    error=zeros(19)
    for i in 1:5
        f=fs[i]
        a=as[i]
        b=bs[i]
        I=exp_values[i]
        for k in 1:19
            N=Ns[k]
            I_n=my.double_exponential(f,a,b,N)
            error[k]=abs(I-I_n)
        end
        plot(Ns,error,marker="o",color=colors[i],label=labels[i])
    end
    xlabel("N",fontsize=15)
    ylabel("\$|I_n-I[f]|\$",fontsize=15)        
    legend(fontsize=15)
    grid()
    yscale("log")
    tight_layout()
    savefig("DE_old.png")
    show()
end

#--------------------------------------------------------------------

# OPTIONAL 4 --------------------------------------------------------

function ex4_optional()
    f(x)=1/(sin(exp(sin(x))))
    a=0
    b=2*pi
    N_max=100
    N_min=10
    values=zeros(N_max-N_min+1)
    Ns=collect(N_min:N_max)
    for N in Ns
        values[N-N_min+1]=my.trapezoidal(f,a,b,N) 
    end
    fig,ax=plt.subplots(ncols=2,figsize=(10,5))
    ax[1].set_title("\$ N \\in [10,100]\$")
    ax[1].plot(Ns,values,color="navy")
    ax[1].set_xlabel("N",fontsize=15)
    ax[1].set_ylabel("\$I_N[f]\$",fontsize=15)
    ax[1].grid()
    ax[2].set_title("\$ N \\in [90,100]\$")
    ax[2].plot(Ns[end-10:end],values[end-10:end],color="navy")
    ax[2].set_ylim(minimum(values[end-10:end])-1e-14,maximum(values[end-10:end])+1e-14)
    ax[2].set_xlabel("N",fontsize=15)
    ax[2].set_ylabel("\$I_N[f]\$",fontsize=15)
    ax[2].grid()
    tight_layout()
    savefig("trapezoidal_periodic.png")
    show()
    println(values[end-10:end])
end

#--------------------------------------------------------------------

# EXECUTION ---------------------------------------------------------

function main()
    println("--- MAIN ---")
    println("Exercise 5.1") #LINES 9-69
    ex1()
    println("-"^30)
    println("\nExercise 5.2") #LINES 75-135
    ex2()
    println("-"^30)
    println("\nExercise 5.3.1") #LINES 141-173
    ex3()
    println("-"^30)
    println("\nExercise 5.3.2") #LINES 179-233
    ex4()
    println("-"^30)
    println("\nExercise 5.4.1") #LINES 239-286
    ex5()
    println("-"^30)
    println("\nExercise 5.4.2") #LINES 292-342
    ex6()
    println("-"^30)
end

function optional()
    println("--- OPTIONAL ---")
    println("Exercise E.1") #LINES 348-389
    ex1_optional()
    println("-"^30)
    println("\nExercise E.2") #LINES 395-439
    ex2_optional()
    println("-"^30)
    println("\nExercise E.3") #LINES 445-479
    ex3_optional()
    println("-"^30)
    println("\nExercise E.4") #LINES 485-512
    ex4_optional()
    println("-"^30)
end

#=
 IT IS POSSIBLE TO ANALYZE SINGLE EXERCISES BY
 REPLACING THE NEXT FEW LINE WITH THE FUNCTION
 RELATED TO THE EXERCISE YOU WANT TO ANALYZE.
=#

main() 
optional()
