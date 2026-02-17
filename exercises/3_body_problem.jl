using Plots

pyplot()

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
    k=1
    while ti[end]<b
        if ti[end]>=k*b/100
            println("$k %")
            k+=1
        end
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

M1=1
M2=1
M3=1

fx(x1,y1,x2,y2)=-(x1-x2)/((sqrt((x1-x2)^2+(y1-y2)^2+1e-10))^3)
fy(x1,y1,x2,y2)=-(y1-y2)/((sqrt((x1-x2)^2+(y1-y2)^2+1e-10))^3)

dot_x1(t,x1,y1,x2,y2,x3,y3,vx1,vy1,vx2,vy2,vx3,vy3)=vx1
dot_vx1(t,x1,y1,x2,y2,x3,y3,vx1,vy1,vx2,vy2,vx3,vy3)=M2*fx(x1,y1,x2,y2)+M3*fx(x1,y1,x3,y3)
dot_y1(t,x1,y1,x2,y2,x3,y3,vx1,vy1,vx2,vy2,vx3,vy3)=vy1
dot_vy1(t,x1,y1,x2,y2,x3,y3,vx1,vy1,vx2,vy2,vx3,vy3)=M2*fy(x1,y1,x2,y2)+M3*fy(x1,y1,x3,y3)

dot_x2(t,x1,y1,x2,y2,x3,y3,vx1,vy1,vx2,vy2,vx3,vy3)=vx2
dot_vx2(t,x1,y1,x2,y2,x3,y3,vx1,vy1,vx2,vy2,vx3,vy3)=M1*fx(x2,y2,x1,y1)+M3*fx(x2,y2,x3,y3)
dot_y2(t,x1,y1,x2,y2,x3,y3,vx1,vy1,vx2,vy2,vx3,vy3)=vy2
dot_vy2(t,x1,y1,x2,y2,x3,y3,vx1,vy1,vx2,vy2,vx3,vy3)=M1*fy(x2,y2,x1,y1)+M3*fy(x2,y2,x3,y3)

dot_x3(t,x1,y1,x2,y2,x3,y3,vx1,vy1,vx2,vy2,vx3,vy3)=vx3
dot_vx3(t,x1,y1,x2,y2,x3,y3,vx1,vy1,vx2,vy2,vx3,vy3)=M2*fx(x3,y3,x2,y2)+M1*fx(x3,y3,x1,y1)
dot_y3(t,x1,y1,x2,y2,x3,y3,vx1,vy1,vx2,vy2,vx3,vy3)=vy3
dot_vy3(t,x1,y1,x2,y2,x3,y3,vx1,vy1,vx2,vy2,vx3,vy3)=M2*fy(x3,y3,x2,y2)+M1*fy(x3,y3,x1,y1)

fs=[dot_x1,dot_y1,dot_x2,dot_y2,dot_x3,dot_y3,dot_vx1,dot_vy1,dot_vx2,dot_vy2,dot_vx3,dot_vy3]
#1->u0s=[0,0,10,5,10,-5,0.03,0.04,-0.05,0,0,0.05]
#2->u0s=[0,0,5*sqrt(3),5,5*sqrt(3),-5,0.03,-0.04,-0.05,0,0,0.05]
#3->u0s=[0,0,0,5,0,10,1.25,0,1,0,0.75,0]
#4->u0s=[1,0,0,5,-1,10,0.25,0,0.25,0.25,0.0,0.25]
#5->u0s=[1,0,0,5,-1,10,0.25,-0.25,0.25,0.25,-0.25,0.25]
#u0s=[1,0,0,5,-1,10,-0.25,0.25,0.5,0.5,0.25,-0.25]
u0s=[0,0,0,10,0,20,0.7,0,1,0,1.3,0]
a=0
b=1000
frames=30
println("Calculating...")
uis,tis=RKDP(fs,u0s,a,b)
#=
println("Animation...")
i=1
len=length(tis)
anim=@animate while i<len
    scatter([uis[1][i]],[uis[2][i]],color=:deepskyblue,legend=false)
    scatter!([uis[3][i]],[uis[4][i]],color=:crimson,legend=false)
    scatter!([uis[5][i]],[uis[6][i]],color=:goldenrod,legend=false)
    plot!(uis[1][1:i],uis[2][1:i],color=:deepskyblue,legend=false)
    plot!(uis[3][1:i],uis[4][1:i],color=:crimson,legend=false)
    plot!(uis[5][1:i],uis[6][1:i],color=:goldenrod,legend=false)
    p=round(100*(i-1)/(length(uis[1])),digits=3)
    println("$p %")
    i_start=i
    while tis[i]<tis[i_start]+10
        global i+=1
        if i>len
            break
        end
    end
end
println("Saving...")
gif(anim,"3bp.gif",fps=frames)
println("Done")
=#

step=div(len,100)
for i in 1:step:len
    scatter([uis[1][i]],[uis[2][i]],color=:deepskyblue,legend=false)
    scatter!([uis[3][i]],[uis[4][i]],color=:crimson,legend=false)
    scatter!([uis[5][i]],[uis[6][i]],color=:goldenrod,legend=false)
    plot!(uis[1][1:1000:i],uis[2][1:1000:i],color=:deepskyblue,legend=false)
    plot!(uis[3][1:1000:i],uis[4][1:1000:i],color=:crimson,legend=false)
    plot!(uis[5][1:1000:i],uis[6][1:1000:i],color=:goldenrod,legend=false)
    savefig("3bodyproblem_$i.png")
end