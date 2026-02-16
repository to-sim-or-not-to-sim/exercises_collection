using Plots
pyplot()

#=
 THIS CODE IS NOT POLISHED AS THE OTHERS,
 BUT FOLLOWS THE SAME BASIC STRUCTURE:
 THERE ARE FUNCTIONS LINKED TO SINGULAR
 EXERCISES, WHICH ARE CALLED AT THE END OF
 THE FILE. YOU CAN MODIFY WHAT YOU SEE BY
 SIMPLY COMMENTING OR NOT COMMENTING THE 
 PARTS YOU DON'T WANT OR WANT TO SEE. I
 DECIDED TO COPY THE USEFUL FUNCTIONS 
 INSTEAD OF IMPORTING THE LIBRARIES BECAUSE
 I USED "PYPLOT" FOR OTHER CODES AND IT IS 
 IN CONFLICT WITH "PLOTS", AND I THOUGHT IT 
 WAS TIME TO TRY TO USE THE NATIVE TOOL.
=#

# USEFUL FUNCTIONS ----------------------------------------------------------------------

function diagm(diagonal::Vector{T}) where T<:Number
    n=length(diagonal)
    A=zeros(Float64,n,n)
    for i in 1:n
        A[i,i]=diagonal[i]
    end
    return A
end

function norm(x::Vector{T}) where T<:Number
    n=length(x)
    sum=0
    for i in 1:n
        sum+=x[i]^2
    end
    return sqrt(sum)
end

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

function secant(f::Function,a::Number,b::Number;k::Number=0,xtol::Float64=100*eps(),ftol::Float64=100*eps(),bracketing::Bool=true,rep_lim::Int64=100)
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
        println(x)
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
    return all_xs[end]
end

function gaussian(x::Number,x_0::Number,sigma::Number,k_0::Number)
    return (1/(2*pi*sigma^2)^(0.25))*exp(-(x-x_0)^2/(4*sigma^2))*exp((0+1im)*k_0*x)
end

#----------------------------------------------------------------------------------------

function ex1_2(x_0,sigma,k_0,a_x,b_x,a_t,b_t)
    delta_x=0.1
    x=collect(a_x:delta_x:b_x)
    delta_t=0.01
    t=collect(a_t:delta_t:b_t)

    psi_0=gaussian.(x,x_0,sigma,k_0)

    lenx=length(x)
    lent=length(t)
    psi=zeros(ComplexF64,lent,lenx)
    psi[1,:]=psi_0
    
    M=zeros(ComplexF64,lenx,lenx)
    a=(0+1im)/(24*delta_x^2)
    g=(0-5im)/(4*delta_x^2)

    norms=zeros(lent)
    norms[1]=delta_x*sum(abs.(psi[1,:]).^2)
    
    #m=1, V=0, \hbar=1
    for i in 1:lenx
        M[i,i]=g
        if i<=lenx-1
            M[i+1,i]=16*a
        end
        if i>=2
            M[i-1,i]=16*a
        end
        if i<=lenx-2
            M[i+2,i]=-a
        end
        if i>=3
            M[i-2,i]=-a
        end
    end

    for j in 2:lent
        k1=M*psi[j-1,:]
        k2=M*(psi[j-1,:].+(delta_t/2).*k1)
        k3=M*(psi[j-1,:].+(delta_t/2).*k2)
        k4=M*(psi[j-1,:].+delta_t*k3)
        psi[j,:]=psi[j-1,:].+(delta_t/6).*(k1.+2*k2.+2*k3.+k4)
        norms[j]=delta_x*sum(abs.(psi[j,:]).^2)
    end
    anim=@animate for i in 1:10:lent
        println(i)
        plot(x,abs.(psi[i,:]).^2,xlim=(a_x,b_x),ylim=(0,abs(gaussian(x_0,x_0,sigma,k_0))^2+0.1),title="t=$(t[i])",titlefontsize=15,legend=false)
    end
    gif(anim,"../ex_$(x_0).gif",fps=10)
    display(plot(t,norms,label="norm")) 

    plot(x,abs.(psi[1,:]).^2,xlim=(a_x,b_x),ylim=(0,abs(gaussian(x_0,x_0,sigma,k_0))^2+0.1),title="t=$(t[1])",titlefontsize=15,legend=false)
    savefig("frame_1_$x_0.png")

    plot(x,abs.(psi[201,:]).^2,xlim=(a_x,b_x),ylim=(0,abs(gaussian(x_0,x_0,sigma,k_0))^2+0.1),title="t=$(t[201])",titlefontsize=15,legend=false)
    savefig("frame_2_$x_0.png")

    plot(x,abs.(psi[301,:]).^2,xlim=(a_x,b_x),ylim=(0,abs(gaussian(x_0,x_0,sigma,k_0))^2+0.1),title="t=$(t[301])",titlefontsize=15,legend=false)
    savefig("frame_3_$x_0.png")

    plot(x,abs.(psi[401,:]).^2,xlim=(a_x,b_x),ylim=(0,abs(gaussian(x_0,x_0,sigma,k_0))^2+0.1),title="t=$(t[401])",titlefontsize=15,legend=false)
    savefig("frame_4_$x_0.png")

    plot(x,abs.(psi[501,:]).^2,xlim=(a_x,b_x),ylim=(0,abs(gaussian(x_0,x_0,sigma,k_0))^2+0.1),title="t=$(t[501])",titlefontsize=15,legend=false)
    savefig("frame_5_$x_0.png")

    plot(x,abs.(psi[1001,:]).^2,xlim=(a_x,b_x),ylim=(0,abs(gaussian(x_0,x_0,sigma,k_0))^2+0.1),title="t=$(t[1001])",titlefontsize=15,legend=false)
    savefig("frame_6_$x_0.png")

    plot(x,abs.(psi[end,:]).^2,xlim=(a_x,b_x),ylim=(0,abs(gaussian(x_0,x_0,sigma,k_0))^2+0.1),title="t=$(t[end])",titlefontsize=15,legend=false)
    savefig("frame_final_$x_0.png")

    plot(t,norms.-1,color=:navy,title="\$|\\psi(x,t)|^2-1\$",titlefontsize=15,legend=false)
    xlabel!("t",fontsize=15)
    savefig("probability_conservation_$x_0.png")

end

function ex3(Vmax,width)
    function V(x::Number)
        if x<-width/2 || x>width/2
            return 0
        end
        return Vmax       
    end
    x_0=-20
    k0=2
    sigma=2
    delta_x=0.2
    x=collect(-100:delta_x:100)
    delta_t=0.005
    t=collect(0:delta_t:40)

    psi_0=gaussian.(x,x_0,sigma,k0)

    lenx=length(x)
    lent=length(t)
    psi=zeros(ComplexF64,lent,lenx)
    psi[1,:]=psi_0
    
    M=zeros(ComplexF64,lenx,lenx)
    a=(0+1im)/(24*delta_x^2)
    g0=(0-5im)/(4*delta_x^2)

    #m=1, V=0, \hbar=1
    for i in 1:lenx
        M[i,i]=g0+(0-1im)*V(x[i])
        if i<=lenx-1
            M[i+1,i]=16*a
        end
        if i>=2
            M[i-1,i]=16*a
        end
        if i<=lenx-2
            M[i+2,i]=-a
        end
        if i>=3
            M[i-2,i]=-a
        end
    end

    for j in 2:lent
        k1=M*psi[j-1,:]
        k2=M*(psi[j-1,:].+(delta_t/2).*k1)
        k3=M*(psi[j-1,:].+(delta_t/2).*k2)
        k4=M*(psi[j-1,:].+delta_t*k3)
        psi[j,:]=psi[j-1,:].+(delta_t/6).*(k1.+2*k2.+2*k3.+k4)
    end
    anim=@animate for i in 1:100:lent
        println(i)
        plot(x,abs.(psi[i,:]).^2,xlim=(-100,100),ylim=(0,0.4),title="t=$(t[i])",titlefontsize=15,legend=false)
        vline!([-width/2,width/2],label=false)
    end
    gif(anim,"ex_3_$Vmax-$width.gif",fps=30)

    plot(x,abs.(psi[1,:]).^2,title="t=$(t[1])",titlefontsize=15,legend=false)
    vline!([-width/2,width/2],label=false)
    savefig("ex3_frame_1_$Vmax-$width.png")
    plot(x,abs.(psi[2001,:]).^2,title="t=$(t[2001])",titlefontsize=15,legend=false)
    vline!([-width/2,width/2],label=false)
    savefig("ex3_frame_2_$Vmax-$width.png")
    plot(x,abs.(psi[6001,:]).^2,title="t=$(t[6001])",titlefontsize=15,legend=false)
    vline!([-width/2,width/2],label=false)
    savefig("ex3_frame_3_$Vmax-$width.png")
    plot(x,abs.(psi[end,:]).^2,title="t=$(t[end])",titlefontsize=15,legend=false)
    vline!([-width/2,width/2],label=false)
    savefig("ex3_frame_4_$Vmax-$width.png")

    psi_end=psi[end,:]
    psi_30=psi[Int(div(30,delta_t)+1),:]

    P_total_30=0
    P_total_end=0
    for i in 1:lenx
        P_total_end+=abs(psi_end[i])^2
        P_total_30+=abs(psi_30[i])^2
    end
    P_total_30*=delta_x
    P_total_end*=delta_x
    println("P_total_end= ",P_total_end,"  P_total_30= ",P_total_30)

    R_end=0
    T_end=0
    R_30=0
    T_30=0
    
    x_sx=collect(-100:delta_x:-width/2)
    len_sx=length(x_sx)
    x_dx=collect(width/2:delta_x:100)
    len_dx=length(x_dx)
    for i in 1:len_sx
        R_end+=abs(psi_end[i])^2
        R_30+=abs(psi_30[i])^2
    end
    R_end*=delta_x
    R_30*=delta_x
    println("R_end= ",R_end/P_total_end,"  R_30= ",R_30/P_total_30)
    for i in 1:len_dx
        T_end+=abs(psi_end[end-i])^2
        T_30+=abs(psi_30[end-i])^2
    end
    T_end*=delta_x
    T_30*=delta_x
    println("T_end= ",T_end/P_total_end,"  T_30= ",T_30/P_total_30)
    println("R_end+T_end= ",(R_end+T_end)/P_total_end)
    println("R_30+T_30= ",(R_30+T_30)/P_total_30)
end

function ex4(V0,omega,x0,sigma,k0,t_final,frames)
    V(x)=V0/(cosh(x/omega)^2)
    
    delta_x=0.2
    x=collect(-100:delta_x:100)
    delta_t=0.005
    t=collect(0:delta_t:t_final)

    psi_0=gaussian.(x,x0,sigma,k0)

    lenx=length(x)
    lent=length(t)
    psi=zeros(ComplexF64,lent,lenx)
    psi[1,:]=psi_0
    
    M=zeros(ComplexF64,lenx,lenx)
    a=(0+1im)/(24*delta_x^2)
    g0=(0-5im)/(4*delta_x^2)

    #m=1, V=0, \hbar=1
    for i in 1:lenx
        M[i,i]=g0+(0-1im)*V(x[i])
        if i<=lenx-1
            M[i+1,i]=16*a
        end
        if i>=2
            M[i-1,i]=16*a
        end
        if i<=lenx-2
            M[i+2,i]=-a
        end
        if i>=3
            M[i-2,i]=-a
        end
    end

    for j in 2:lent
        k1=M*psi[j-1,:]
        k2=M*(psi[j-1,:].+(delta_t/2).*k1)
        k3=M*(psi[j-1,:].+(delta_t/2).*k2)
        k4=M*(psi[j-1,:].+delta_t*k3)
        psi[j,:]=psi[j-1,:].+(delta_t/6).*(k1.+2*k2.+2*k3.+k4)
    end

    anim=@animate for i in 1:100:lent
        println(i)
        plot(x,abs.(psi[i,:]).^2,xlim=(-100,100),ylim=(-0.1,0.4),title="t=$(t[i])",titlefontsize=15,legend=false)
        plot!(x,V.(x),label=false)
    end
    gif(anim,"ex_4_$V0-$omega.gif",fps=30)

    for tk in frames
        plot(x,abs.(psi[200*tk+1,:]).^2,ylim=(-0.1,maximum(abs.(psi[200*tk+1,:]).^2)+0.1),title="t=$(t[200*tk+1])",titlefontsize=15,legend=false)
        plot!(x,V.(x),label=false)
        savefig("ex4_$tk-$V0-$omega.png")
    end

    psi_end=psi[end,:]

    R_end=0
    T_end=0
    P_total_end=0
    index_dx=div(lenx,2)
    index_sx=div(lenx,2)
    while V(x[index_dx])>=0.01*V0 &&  V(x[index_sx])>=0.01*V0
        index_sx-=1
        index_dx+=1
    end
    x_sx=collect(-100:delta_x:x[index_sx])
    len_sx=length(x_sx)
    x_dx=collect(x[index_dx]:delta_x:100)
    len_dx=length(x_dx)

    for i in 1:len_sx
        P_total_end+=abs(psi_end[i])^2
    end
    for i in len_dx:lenx
        P_total_end+=abs(psi_end[i])^2
    end
    P_total_end*=delta_x
    println("P_total_end= ",P_total_end)

    for i in 1:len_sx
        R_end+=abs(psi_end[i])^2
    end
    R_end*=delta_x
    println("R_end= ",R_end/P_total_end)
    for i in 1:len_dx
        T_end+=abs(psi_end[end-i])^2
    end
    T_end*=delta_x
    println("T_end= ",T_end/P_total_end)
    println("R_end+T_end= ",(R_end+T_end)/P_total_end)
end

function ex5()
    function V(x::Number)
        if x<-1 || x>1
            return 0
        end
        return -10     
    end
    
    f_odd(E)=tan(sqrt((10+E)*2))+sqrt(-(10+E)/E)
    f_even(E)=tan(sqrt((10+E)*2))-sqrt(-E/(10+E))
    
    delta_x=0.05
    x=collect(-5:delta_x:5)
    lenx=length(x)
    H=zeros(lenx,lenx)
    b=-1/(24*delta_x^2)
    h0=5/(4*delta_x^2)
    for i in 1:lenx
        H[i,i]=h0+V(x[i])
        if i<=lenx-1
            H[i+1,i]=16*b
        end
        if i>=2
            H[i-1,i]=16*b
        end
        if i<=lenx-2
            H[i+2,i]=-b
        end
        if i>=3
            H[i-2,i]=-b
        end
    end
    diag,eigenvectors=eigenvalues(H,k=1000)
    indexes=[]
    for k in 1:size(diag)[1]
        if diag[k,k]>-10 && diag[k,k]<0
            push!(indexes,k)
        end
    end
    E_n=zeros(length(indexes))
    for i in 1:length(indexes)
        E_n[i]=diag[indexes[i],indexes[i]]
    end
    phi_n=eigenvectors[:,indexes]
    for i in 1:length(E_n)
        eps=E_n[i]
        plt=plot(x,(phi_n[:,i]),color=:navy,label="\$\\varphi_n(x)\$",grid=true,legend=true,legendfontsize=15,title="\$E_n=$eps\$",titlefontsize=15)
        hline!([0],label="",color=:slategray,linestyle=:dash)
        vline!([1,-1],label="",color=:slategray,linestyle=:dash)
        xlabel!("x",fontsize=15)
        ylabel!("\$\\varphi_n(x)\$",fontsize=15)
        savefig("eigenfunction_$i.png")
        E_exp=0
        if -8<eps<-4
            #odd solution => cot(ka)=-rho/k
            E_exp=secant(f_odd,eps-0.1,eps+0.1)
        else
            #even solution => tan(ka)=rho/k
            E_exp=secant(f_even,eps-0.1,eps+0.1)
        end
        println("Exp E= ",E_exp," Found E= ",eps)
    end

end

# IT MAKE TAKE A LONG TO EXECUTE ALL THE EXERCISES TOGETHER

ex1_2(0,1,0,-10,10,0,20)
ex1_2(20,2,2,0,100,0,20)
println("V=2, -1<x<1")
ex3(2,2) 
println("V=5, -1<x<1")
ex3(5,2)
println("V=2, -0.5<x<0.5")
ex3(2,1)
ex4(2,1,-20,2,2,40,[0,10,30])
ex4(-2,0.5,-20,4,1,80,[0,20,60])
ls=collect(1:10)
for l in ls
    V0=-2*l*(l+1)
    ex4(V0,0.5,-20,4,1,80,[0,20,60])
end
ex5()