include("../MyLib.jl")
using PyPlot,PrettyTables,SpecialFunctions
import .MyLib.MyODE as my
using .MyLib.MyInterpolation: inf_norm

# EX 6.1.1 ----------------------------------------------------------

function ex1()
    fa(t,u)=-2*t*u
    sol_a(t)=2*(exp(-t^2))
    fb(t,u)=u+t
    sol_b(t)=-1-t+3*exp(t)
    fc(t,u)=(t^2/((1+t^3)*u))
    sol_c(t)=sqrt(1+(2/3)*log(1+t^3))
    fs=[fa,fb,fc]
    u0s=[2,2,1]
    a=0
    bs=[2,1,3]
    sol_fs=[sol_a,sol_b,sol_c]
    titles=["a","b","c"]

    fig,ax=subplots(ncols=3,figsize=(12,4))
    for i in 1:3
        ax[i].set_title(titles[i],fontsize=15)
        f=fs[i]
        u0=u0s[i]
        b=bs[i]
        sol=sol_fs[i]
        uis,tis=my.euler([f],[u0],a,b,320)
        ax[i].plot(tis,sol.(tis),color="navy",label="\$\\hat u(t)\$")
        ax[i].scatter(tis,uis[1,:],s=10,color="darkorange",label="\$u_i\$")
        ax[i].grid()
        ax[i].set_xlabel("t",fontsize=15)
        ax[i].legend(fontsize=15)
    end
    tight_layout()
    savefig("euler1.png")
    show()

    k=collect(2:10)
    ns=@. 10*(2^k)

    fig,ax=subplots(ncols=3,figsize=(12,4))
    for i in 1:3
        f=fs[i]
        u0=u0s[i]
        b=bs[i]
        sol=sol_fs[i]
        error=zeros(length(ns))
        final_error=zeros(length(ns))
        for k in 1:length(ns)
            n=ns[k]
            uis,tis=my.euler([f],[u0],a,b,n)
            error[k]=inf_norm(uis[1,:],sol.(tis))
            final_error[k]=@. abs(uis[1,end]-sol(tis[end]))
        end
        ax[i].set_title(titles[i])
        ax[i].plot((b-a)./ns,error,marker="o",label="\$||u-\\hat{u}||_\\infty\$")
        ax[i].plot((b-a)./ns,final_error,marker="o",label="\$|u_n-\\hat{u}(t_n)|\$")
        ax[i].plot((b-a)./ns,(b-a)./ns,linestyle="--",label="\$o(h)\$")
        ax[i].set_xscale("log")
        ax[i].set_yscale("log")
        ax[i].set_xlabel("h",fontsize=15)
        ax[i].set_ylabel("\$||u-\\hat u||_\\infty\$",fontsize=15)
        ax[i].grid()
        ax[i].legend(fontsize=15)
    end
    tight_layout()
    savefig("euler2.png")
    show()
end

#--------------------------------------------------------------------

# EX 6.1.2 ----------------------------------------------------------

function ex2()
    titles=["\$y''+9y=\\sin(2t)\$","\$y''-4y=4t\$","\$y''+4y'+4y=t\$"]
    fa(t,u1,u2)=sin(2*t)-9*u1
    sol_a(t)=0.2*sin(3*t)+2*cos(3*t)+0.2*sin(2*t)
    sold_a(t)=0.6*cos(3*t)-6*sin(3*t)+0.4*cos(2*t)
    fb(t,u1,u2)=4*t+4*u1
    sol_b(t)=exp(2*t)+exp(-2*t)-t
    sold_b(t)=2*exp(2*t)-2*exp(-2*t)-1
    fc(t,u1,u2)=t-4*u1-4*u2
    sol_c(t)=(3*t+1.25)*exp(-2*t)+0.25*(t-1)
    sold_c(t)=3*exp(-2*t)-2*(3*t+1.25)*exp(-2*t)+0.25
    prime(t,u1,u2)=u2
    fs=[fa,fb,fc]
    sol_fds=[sold_a,sold_b,sold_c]
    u0s=[[2,1],[2,-1],[1,0.75]]
    a=0
    bs=[2*pi,1.5,4]
    sol_fs=[sol_a,sol_b,sol_c]
    fig,ax=subplots(ncols=3,nrows=2,figsize=(12,6))
    for i in 1:3
        f=fs[i]
        u0=u0s[i]
        b=bs[i]
        sol=sol_fs[i]
        sold=sol_fds[i]
        n=1000
        uis,tis=my.euler([prime,f],u0,a,b,n)
        ax[1,i].set_title(titles[i],fontsize=15)
        ax[1,i].plot(tis,sol.(tis),color="navy",linestyle="--",label="\$\\hat{y}(t)\$")
        ax[1,i].plot(tis,uis[1,:],color="navy",label="y(t)")
        ax[1,i].plot(tis,sold.(tis),color="orangered",linestyle="--",label="\$\\hat{y}'(t)\$")
        ax[1,i].plot(tis,uis[2,:],color="orangered",label="y'(t)")
        ax[1,i].grid()
        ax[2,i].plot(tis,abs.(uis[1,:].-sol.(tis)),color="navy",label="\$|y(x)-\\hat y(x)|\$")
        ax[2,i].plot(tis,abs.(uis[2,:].-sold.(tis)),color="orangered",label="\$|y'(x)-\\hat y'(x)|\$")
        ax[2,i].set_xlabel("t",fontsize=15)
        ax[1,i].set_xlabel("t",fontsize=15)
        ax[2,i].grid()
        if i==3
            ax[1,i].legend(fontsize=15,loc=(1,0.25))
            ax[2,i].legend(fontsize=15,loc=(0.75,0.5))
        end
    end
    tight_layout()
    savefig("euler3.png")
    show()    
end

#--------------------------------------------------------------------

# EX 6.2.1 ----------------------------------------------------------

function ex3()
    f(t,u)=-2*t*u
    sol(t)=2*(exp(-t^2))
    u0=2
    a=0
    b=2
    ns=collect(30:30:300)
    error_IE2=zeros(length(ns))
    error_RK4=zeros(length(ns))
    for k in 1:length(ns)
        n=ns[k]
        uis,tis=my.IE2([f],[u0],a,b,n)
        error_IE2[k]=my.inf_norm(vec(uis),sol.(tis))
        uis,tis=my.RK4([f],[u0],a,b,n)
        error_RK4[k]=my.inf_norm(vec(uis),sol.(tis))
    end 
    fig,ax=subplots(ncols=2,figsize=(10,5))
    ax[1].plot(2 .*ns,error_IE2,color="indigo",marker="o",label="\$||y-\\hat{y}||_\\infty\$")
    ax[1].plot(2 .*ns,1 ./(ns.^2),color="silver",linestyle="--",label="\$o(h^2)\$")
    ax[1].set_xlabel("Number of function evaluations",fontsize=15)
    ax[1].set_title("IE2",fontsize=15)
    ax[1].legend(fontsize=15)
    ax[1].set_xscale("log")
    ax[1].set_yscale("log")
    ax[1].grid()
    ax[2].plot(4 .*ns,error_RK4,color="indigo",marker="o",label="\$||y-\\hat{y}||_\\infty\$")
    ax[2].plot(4 .*ns,1 ./(ns.^4),color="silver",linestyle="--",label="\$o(h^4)\$")
    ax[2].set_xlabel("Number of function evaluations",fontsize=15)
    ax[2].set_title("RK4",fontsize=15)
    ax[2].legend(fontsize=15)
    ax[2].set_xscale("log")
    ax[2].set_yscale("log")
    ax[2].grid()
    savefig("IE2vsRK4.png")
    show()
    number_of_evaluation=collect(1:100:10000)
    plot(number_of_evaluation,1 ./(number_of_evaluation.^2),color="navy",linestyle="--",label="\$o(h^2)\$")
    plot(number_of_evaluation,1 ./(number_of_evaluation.^4),color="darkorange",linestyle="--",label="\$o(h^4)\$")
    xlabel("Number of function evaluations",fontsize=15)
    legend(fontsize=15)
    grid()
    xscale("log")
    yscale("log")
    savefig("h4vsh2.png")
    show()
end

#--------------------------------------------------------------------

# EX 6.2.2 ----------------------------------------------------------

function ex4()

    function exp1()
        fx(t,x,y)=-4*y+x*(1-x^2-y^2)
        fy(t,x,y)=4*x+y*(1-x^2-y^2)
        x0s=[0.1,0]
        y0s=[0,1.9]
        theta_0=[0,pi/2]
        r_0=[0.1,1.9]
        a=0
        b=10
        fig,ax=subplots(ncols=2,nrows=2,figsize=(10,10))
        n=100000
        for i in 1:2
            x0=x0s[i]
            y0=y0s[i]
            uis,tis=my.RK4([fx,fy],[x0,y0],a,b,n)
            ax[1,i].set_title("\$(x_0,y_0)=\$($x0,$y0)",fontsize=15)
            ax[1,i].plot(uis[1,:],uis[2,:],color="deepskyblue")
            ax[1,i].set_xlabel("x",fontsize=15)
            ax[1,i].set_ylabel("y",fontsize=15)
            ax[1,i].grid()

            t=(10/n)*collect(0:n)
            x=@. sqrt(1/(((1-r_0[i]^2)/r_0[i]^2)*exp(-2*t)+1))*cos(4*t+theta_0[i])
            y=@. sqrt(1/(((1-r_0[i]^2)/r_0[i]^2)*exp(-2*t)+1))*sin(4*t+theta_0[i])
            ax[2,i].plot(t,abs.(uis[1,:].-x),color="navy",label="x")
            ax[2,i].plot(t,abs.(uis[2,:].-y),color="darkorange",label="y")
            ax[2,i].grid()
            ax[2,i].set_ylabel("error",fontsize=15)
            ax[2,i].set_xlabel("t",fontsize=15)
            ax[2,i].legend(fontsize=15)
        end
        tight_layout()
        savefig("RK4_1.png")
        show()
    end

    function exp2()
        fx(t,x,y)=-4*y-0.25*x*(1-x^2-y^2)*(4-x^2-y^2)
        fy(t,x,y)=4*x-0.25*y*(1-x^2-y^2)*(4-x^2-y^2)
        x0s=[0.95,0,-2.5]
        y0s=[0,1.05,0]
        a=0
        b=10
        fig,ax=subplots(ncols=3,nrows=2,figsize=(9,6))
        n=100000
        for i in 1:3
            x0=x0s[i]
            y0=y0s[i]
            uis,tis=my.RK4([fx,fy],[x0,y0],a,b,n)
            ax[1,i].set_title("\$(x_0,y_0)=\$($x0,$y0)",fontsize=15)
            ax[1,i].plot(uis[1,:],uis[2,:],color="deepskyblue")
            ax[1,i].set_xlabel("x",fontsize=15)
            ax[1,i].set_ylabel("y",fontsize=15)
            ax[1,i].grid()

            t=(10/n)*collect(0:n)
            ax[2,i].plot(t,sqrt.(uis[1,:].^2 .+ uis[2,:].^2),color="navy")
            ax[2,i].set_xlabel("t",fontsize=15)
            ax[2,i].set_ylabel("\$ \\sqrt{x^2(t)+y^2(t)} \$",fontsize=15)
            ax[2,i].grid()
        end
        tight_layout()
        savefig("RK4_2.png")
        show()
    end

    function test()
        fx(t,x,y)=-4*y-0.25*x*(1-x^2-y^2)*(4-x^2-y^2)
        fy(t,x,y)=4*x-0.25*y*(1-x^2-y^2)*(4-x^2-y^2)
        x0=1
        y0=0
        a=0
        b=10
        n=100000
        fig,ax=subplots(ncols=2,figsize=(6,3))
        uis,tis=my.RK4([fx,fy],[x0,y0],a,b,n)
        ax[1].plot(uis[1,:],uis[2,:],color="deepskyblue")
        ax[1].set_xlabel("x",fontsize=15)
        ax[1].set_ylabel("y",fontsize=15)
        ax[1].grid()
        t=(10/n)*collect(0:n)
        ax[2].plot(t,sqrt.(uis[1,:].^2 .+ uis[2,:].^2),color="navy")
        ax[2].set_xlabel("t",fontsize=15)
        ax[2].set_ylabel("\$ \\sqrt{x^2(t)+y^2(t)} \$",fontsize=15)
        ax[2].grid()
        tight_layout()
        savefig("RK4_2_problems.png")
        show()
    end

    exp1()
    exp2()
    test()
end

#--------------------------------------------------------------------

# EX 6.2.3 ----------------------------------------------------------

function ex5(;alpha::Number=0.1,beta::Number=0.25)
    fy(t,y,z)=y*(1-alpha*y)-y*z/(1+beta*y)
    fz(t,y,z)=-z+y*z/(1+beta*y)
    a=0
    b=60
    y0=1
    z0=0.01
    n=10000
    uis,tis=my.RK4([fy,fz],[y0,z0],a,b,n)
    fig,ax=subplots(ncols=2,figsize=(10,5))
    ax[1].plot(uis[1,:],uis[2,:],color="deepskyblue")
    ax[1].set_xlabel("y",fontsize=15)
    ax[1].set_ylabel("z",fontsize=15)
    ax[1].grid()
    ax[2].plot(tis,uis[1,:],color="forestgreen",label="y(t)")
    ax[2].plot(tis,uis[2,:],color="crimson",label="z(t)")
    ax[2].set_xlabel("t",fontsize=15)
    ax[2].grid()
    ax[2].legend(fontsize=15,loc=(1,0.5))
    tight_layout()
    savefig("predator_prey.png")
    show()
end

#--------------------------------------------------------------------

# OPTIONAL 1 --------------------------------------------------------

function ex1_optional()
    fa(t,u)=-2*t*u
    sol_a(t)=2*(exp(-t^2))
    fb(t,u)=u+t
    sol_b(t)=-1-t+3*exp(t)
    fc(t,u)=(t^2/((1+t^3)*u))
    sol_c(t)=sqrt(1+(2/3)*log(1+t^3))
    fs=[fa,fb,fc]
    u0s=[2,2,1]
    a=0
    bs=[2,1,3]
    sol_fs=[sol_a,sol_b,sol_c]
    titles=["a","b","c"]

    fig,ax=subplots(ncols=3,figsize=(12,4))
    for i in 1:3
        ax[i].set_title(titles[i],fontsize=15)
        f=fs[i]
        u0=u0s[i]
        b=bs[i]
        sol=sol_fs[i]
        ui,ti=my.backward_euler(f,u0,a,b,320)
        ax[i].plot(ti,sol.(ti),color="navy",label="\$\\hat u(t)\$")
        ax[i].scatter(ti,ui,s=10,color="darkorange",label="\$u_i\$")
        ax[i].grid()
        ax[i].set_xlabel("t",fontsize=15)
        ax[i].legend(fontsize=15)
    end
    tight_layout()
    savefig("backward_euler1.png")
    show()

    k=collect(2:10)
    ns=@. 10*(2^k)

    error_back=zeros(3,length(ns))
    final_error=zeros(length(ns))
    error_euler=zeros(3,length(ns))
    error_average=zeros(3,length(ns))

    fig,ax=subplots(ncols=3,figsize=(12,4))
    for i in 1:3
        f=fs[i]
        u0=u0s[i]
        b=bs[i]
        sol=sol_fs[i]
        for k in 1:length(ns)
            n=ns[k]
            ui,ti=my.backward_euler(f,u0,a,b,n)    
            uis_euler,tis=my.euler([f],[u0],a,b,n)
            ui_euler=uis_euler[1,:]
            error_back[i,k]=inf_norm(ui,sol.(ti))
            final_error[k]=@. abs(ui[end]-sol(ti[end]))
            error_euler[i,k]=inf_norm(ui_euler,sol.(ti))
            error_average[i,k]=inf_norm((ui+ui_euler)/2,sol.(ti))
        end
        ax[i].set_title(titles[i])
        ax[i].plot((b-a)./ns,error_back[i,:],marker="o",label="\$||u-\\hat{u}||_\\infty\$")
        ax[i].plot((b-a)./ns,final_error,marker="o",label="\$|u_n-\\hat{u}(t_n)|\$")
        ax[i].plot((b-a)./ns,(b-a)./ns,linestyle="--",label="\$o(h)\$")
        ax[i].set_xscale("log")
        ax[i].set_yscale("log")
        ax[i].set_xlabel("h",fontsize=15)
        ax[i].set_ylabel("\$||u-\\hat u||_\\infty\$",fontsize=15)
        ax[i].grid()
        ax[i].legend(fontsize=15)
    end
    tight_layout()
    savefig("backward_euler2.png")
    show()
    for i in 1:3
        title(titles[i],fontsize=15)
        plot(ns,error_back[i,:],marker="o",color="darkorange",label="Backward Euler")
        plot(ns,error_euler[i,:],marker="o",color="lime",label="Euler")
        plot(ns,error_average[i,:],marker="o",color="rebeccapurple",label="Average")
        plot(ns,1 ./ns.^2,color="navy",linestyle="--",label="\$o(n^{-2})\$")
        plot(ns,1 ./ns,color="deepskyblue",linestyle="--",label="\$o(n^{-1})\$")
        xscale("log")
        yscale("log")
        xlabel("n",fontsize=15)
        ylabel("\$ ||u-\\hat u||_\\infty \$",fontsize=15)
        legend(fontsize=15,loc=(0,-0.25))
        grid()
        tight_layout()
        savefig("backward_vs_euler$i.png")
        show()
    end
end

#--------------------------------------------------------------------

# OPTIONAL 2 --------------------------------------------------------

function ex2_optional()
    fv(t,v,w)=0.2*(1-v)-3*v*w
    fw(t,v,w)=(3*v-1)*w
    v0=0.95
    w0=0.05
    a=0
    b=60
    n=600000
    uis,tis=my.RK4([fv,fw],[v0,w0],a,b,n)
    fig,ax=subplots(ncols=2,figsize=(10,5))
    ax[1].plot(uis[1,:],uis[2,:],color="mediumvioletred")
    ax[1].set_xlabel("v",fontsize=15)
    ax[1].set_ylabel("w",fontsize=15)
    ax[1].grid()
    ax[2].plot(tis,uis[1,:],color="royalblue",label="v(t)")
    ax[2].plot(tis,uis[2,:],color="darkolivegreen",label="w(t)")
    ax[2].set_xlabel("t",fontsize=15)
    ax[2].grid()
    ax[2].legend(fontsize=15)
    savefig("SIR.png")
    show()
    println("v:",uis[1,end],"\nw:",uis[2,end])
end

#--------------------------------------------------------------------

# OPTIONAL 3 --------------------------------------------------------

function ex3_optional()
    f(t,u)=-2*t*u
    a=0
    b=2
    u0=2
    sol(t)=2*exp(-t^2)
    uis,ti=my.BS23([f],[u0],a,b,tol=1e-8)
    plot(ti,uis[1],color="darkorange",label="\$u_i\$")
    plot(ti,sol.(ti),linestyle="--",color="navy",label="\$\\hat u\$")
    legend(fontsize=15)
    xlabel("t",fontsize=15)
    ylabel("f(t)",fontsize=15)
    grid()
    savefig("BS23_1.png")
    show()
    plot(ti,abs.(uis[1]-sol.(ti)))
    yscale("log")
    grid()
    show()
    hi=ti[3:end-1]-ti[2:end-2]
    plot(ti[3:end-1],hi,color="indigo")
    xlabel("t",fontsize=15)
    ylabel("h",fontsize=15)
    grid()
    tight_layout()
    savefig("h_t.png")
    show()
    delta=[1e-6,1e-7,1e-8,1e-9,1e-10,1e-11,1e-12]
    error_last=zeros(length(delta))
    error=zeros(length(delta))
    k=1
    for tol in delta
        uis,ti=my.BS23([f],[u0],a,b,tol=tol)
        error_last[k]=abs(uis[1][end]-sol(b))
        error[k]=inf_norm(uis[1],sol.(ti))
        k+=1
    end
    plot(delta,error_last,color="darkorange",label="\$||\\hat u(b)-u_n||\$")
    plot(delta,error,color="navy",label="\$||\\hat u(t)-u_i||\$")
    plot(delta,delta.^(1/2),color="gray",linestyle="--",label="\$o(\\delta^{1/2})\$")
    xlabel("\$\\delta\$",fontsize=15)
    ylabel("\$||\\cdot||_\\infty\$",fontsize=15)
    yscale("log")
    xscale("log")
    legend(fontsize=15)
    grid()
    tight_layout()
    savefig("BS23_error.png")
    show()
end

#--------------------------------------------------------------------

# OPTIONAL 4 --------------------------------------------------------

function ex4_optional()
    u_prime(t,u,v)=v
    v_prime(t,u,v)=-u*(1+v)^3
    u0s=[0.1,0.5,0.75,0.95]
    v0=0
    a=0
    b=4*pi
    fs=[u_prime,v_prime]
    for k in 1:4
        fig,ax=subplots(ncols=3,figsize=(12,4))
        u0=u0s[k]
        uis,ti=my.BS23(fs,[u0,v0],a,b,tol=1e-8)
        yis=uis[1]
        yis_prime=uis[2]
        fig.suptitle("\$y_0=$u0\$   \$y'_0=$v0\$",fontsize=15)
        ax[1].set_title("\$y_i\$",fontsize=15)
        ax[1].plot(ti,yis,color="indigo")
        ax[1].set_xlabel("x",fontsize=15)
        ax[1].set_ylabel("y",fontsize=15)
        ax[1].grid()
        ax[2].set_title("\$y'_i\$",fontsize=15)
        ax[2].plot(ti,yis_prime,color="indigo")
        ax[2].set_xlabel("x",fontsize=15)
        ax[2].set_ylabel("y",fontsize=15)
        ax[2].grid()
        hi=ti[3:end-1]-ti[2:end-2]
        ax[3].set_title("\$h(t)\$",fontsize=15)
        ax[3].plot(ti[3:end-1],hi,color="indigo")
        ax[3].set_xlabel("t",fontsize=15)
        ax[3].set_ylabel("h",fontsize=15)
        ax[3].grid()
        tight_layout()
        savefig("BS23_2_$k")
        show()
        println("min(h)= ",minimum(hi))
        println("mean(h)= ",sum(hi)/length(hi))
    end
end

#--------------------------------------------------------------------

# OPTIONAL 5 --------------------------------------------------------

function ex5_optional()
    f(t,u)=100*u^2-u^3
    u0=0.0002
    a=0
    b=100
    fig,ax=subplots(ncols=2,figsize=(8,4))
    uis,ti=my.BS23([f],[u0],a,b,tol=1e-8)
    fig.suptitle("\$y_0=$u0\$",fontsize=15)
    ax[1].set_title("\$y_i\$",fontsize=15)
    ax[1].plot(ti,uis[1],color="indigo")
    ax[1].set_xlabel("u",fontsize=15)
    ax[1].set_ylabel("f(t)",fontsize=15)
    ax[1].grid()
    hi=ti[3:end-1]-ti[2:end-2]
    ax[2].set_title("\$h(t)\$",fontsize=15)
    ax[2].plot(ti[3:end-1],hi,color="indigo")
    println(hi[end-20:end])
    ax[2].set_xlabel("t",fontsize=15)
    ax[2].set_ylabel("h",fontsize=15)
    ax[2].grid()
    tight_layout()
    savefig("BS23_const.png")
    show()
end

#--------------------------------------------------------------------

# OPTIONAL 6 --------------------------------------------------------

function ex6_optional()
    f(t,u)=-2*t*u
    a=0
    b=2
    u0=2
    sol(t)=2*exp(-t^2)
    uis,ti=my.RKDP([f],[u0],a,b,tol=1e-8)
    plot(ti,uis[1],color="darkorange",label="\$u_i\$")
    plot(ti,sol.(ti),linestyle="--",color="navy",label="\$\\hat u\$")
    legend(fontsize=15)
    xlabel("t",fontsize=15)
    ylabel("f(t)",fontsize=15)
    grid()
    savefig("RKDP_1.png")
    show()
    plot(ti,abs.(uis[1]-sol.(ti)))
    yscale("log")
    grid()
    show()
    hi=ti[3:end-1]-ti[2:end-2]
    plot(ti[3:end-1],hi,color="indigo")
    xlabel("t",fontsize=15)
    ylabel("h",fontsize=15)
    grid()
    tight_layout()
    savefig("h_t_RKDP.png")
    show()
    delta=[1e-6,1e-7,1e-8,1e-9,1e-10,1e-11,1e-12]
    error_last=zeros(length(delta))
    error=zeros(length(delta))
    k=1
    for tol in delta
        uis,ti=my.RKDP([f],[u0],a,b,tol=tol)
        error_last[k]=abs(uis[1][end]-sol(b))
        error[k]=inf_norm(uis[1],sol.(ti))
        k+=1
    end
    plot(delta,error_last,color="darkorange",label="\$||\\hat u(b)-u_n||\$")
    plot(delta,error,color="navy",label="\$||\\hat u(t)-u_i||\$")
    plot(delta,delta.^(1/2),color="gray",linestyle="--",label="\$o(\\delta^{1/2})\$")
    xlabel("\$\\delta\$",fontsize=15)
    ylabel("\$||\\cdot||_\\infty\$",fontsize=15)
    yscale("log")
    xscale("log")
    legend(fontsize=15)
    grid()
    tight_layout()
    savefig("RKDP_error.png")
    show()
end

#--------------------------------------------------------------------

# OPTIONAL 7 --------------------------------------------------------

function ex7_optional()
    u_prime(t,u,v)=v
    v_prime(t,u,v)=-u*(1+v)^3
    u0s=[0.1,0.5,0.75,0.95]
    v0=0
    a=0
    b=4*pi
    fs=[u_prime,v_prime]
    for k in 1:4
        fig,ax=subplots(ncols=3,figsize=(12,4))
        u0=u0s[k]
        uis,ti=my.RKDP(fs,[u0,v0],a,b,tol=1e-8)
        yis=uis[1]
        yis_prime=uis[2]
        fig.suptitle("\$y_0=$u0\$   \$y'_0=$v0\$",fontsize=15)
        ax[1].set_title("\$y_i\$",fontsize=15)
        ax[1].plot(ti,yis,color="indigo")
        ax[1].set_xlabel("x",fontsize=15)
        ax[1].set_ylabel("y",fontsize=15)
        ax[1].grid()
        ax[2].set_title("\$y'_i\$",fontsize=15)
        ax[2].plot(ti,yis_prime,color="indigo")
        ax[2].set_xlabel("x",fontsize=15)
        ax[2].set_ylabel("y",fontsize=15)
        ax[2].grid()
        hi=ti[3:end-1]-ti[2:end-2]
        ax[3].set_title("\$h(t)\$",fontsize=15)
        ax[3].plot(ti[3:end-1],hi,color="indigo")
        ax[3].set_xlabel("t",fontsize=15)
        ax[3].set_ylabel("h",fontsize=15)
        ax[3].grid()
        tight_layout()
        savefig("RKDP_2_$k")
        show()
        println("min(h)= ",minimum(hi))
        println("mean(h)= ",sum(hi)/length(hi))
    end
end

#--------------------------------------------------------------------

# EXECUTION ---------------------------------------------------------

function main()
    println("--- MAIN ---")
    println("Exercise 6.1.1") #LINES 8-71
    ex1()
    println("-"^30)
    println("\nExercise 6.1.2") #LINES 77-123
    ex2()
    println("-"^30)
    println("\nExercise 6.2.1") #LINES 129-174
    ex3()
    println("-"^30)
    println("\nExercise 6.2.2") #LINES 180-275
    ex4()
    println("-"^30)
    println("\nExercise 6.2.3") #LINES 281-303
    ex5()
    println("-"^30)
end

function optional()
    println("--- OPTIONAL ---")
    println("Exercise F.1") #LINES 309-396
    ex1_optional()
    println("-"^30)
    println("\nExercise F.2") #LINES 402-424
    ex2_optional()
    println("-"^30)
    println("\nExercise F.3.1") #LINES 430-479
    ex3_optional()
    println("-"^30)
    println("\nExercise F.3.2") #LINES 485-522
    ex4_optional()
    println("-"^30)
    println("\nExercise F.3.3") #LINES 528-551
    ex5_optional()
    println("-"^30)
    println("\nExercise F.4.1") #LINES 557-606
    ex6_optional()
    println("-"^30)
    println("\nExercise F.4.2") #LINES 612-649
    ex7_optional()
    println("-"^30)
end

#=
 IT IS POSSIBLE TO ANALYZE SINGLE EXERCISES BY
 REPLACING THE NEXT FEW LINE WITH THE FUNCTION
 RELATED TO THE EXERCISE YOU WANT TO ANALYZE.
=#

main() 
optional()
