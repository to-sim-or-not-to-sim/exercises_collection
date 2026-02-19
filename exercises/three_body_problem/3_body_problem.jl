using Plots, JSON

pyplot()

# --- CAN BE OPTIMIZED A LOT ---#

function RKDP(fs::AbstractVector,u0s::Vector{T},a::Number,b::Number;tol::Float64=1e-12) where T<:Number
    step=(b-a)/5000
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
    last_us=zeros(len)
    last_us[:]=u0s[:]
    last_t=a
    last_saved_t=a
    u_4=zeros(len)
    u_5=zeros(len)
    u_step=zeros(len)
    ks=zeros(7,len) 
    for i in 1:len
        push!(uis,Float64[u0s[i]])
    end
    ti=Float64[]
    push!(ti,a)
    k=1
    while last_t<b
        if last_t>=k*b/100
            println("$k %")
            k+=1
        end
        for i in 1:7
            t_step=last_t+C[i]*h
            for l in 1:len
                if i==1
                    u_step[l]=last_us[l]
                else
                    u_step[l]=sum(A[i,1:i].*ks[1:i,l]')+last_us[l] 
                end                  
            end
            for l in 1:len
                f=fs[l]
                ks[i,l]=h*f(t_step,u_step...)
            end
        end

        for l in 1:len
            u_4[l]=sum(B_4.*ks[:,l]')+last_us[l]
            u_5[l]=sum(B_5.*ks[:,l]')+last_us[l]
        end
        errs=abs.(u_4-u_5)
        err=maximum(errs)
        if err<tol
            add=false
            if last_t>=last_saved_t+step
                add=true
            end     
            for l in 1:len
                last_us[l]=u_5[l]
                if add
                    push!(uis[l],u_5[l])
                end
            end
            if add
                push!(ti,last_t) 
                last_saved_t=last_t
            end        
            last_t+=h
        end
        
        h=0.8*h*((tol/err)^(1/5))
        h=minimum([h,h_max])
        if last_t+h==last_t
            println("Error during execution: step zero.")
            return uis,ti
        end 
        if last_t+h>b
            h=b-last_t
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
#u0s=[0,0,0,10,0,20,0.7,0,1,0,1.3,0] #<---- E>0
#u0s=[0,0,0.5,-sqrt(3)/2,-0.5,-sqrt(3)/2,1,0,-1/2,-sqrt(3)/2,-1/2,sqrt(3)/2] #<---- Lagrange
#u0s=[0,0,0.5,-sqrt(3)/2,-0.5,-sqrt(3)/2,sqrt(2),0,-1/sqrt(2),-sqrt(3/2),-1/sqrt(2),sqrt(3/2)] #<---- E=0
u0s=[0,0,1,0,-1,0,0,0,0,sqrt(5)/2,0,-sqrt(5)/2] #<---- Euler
#u0s=[0.97000436,-0.24308753,-0.97000436,0.24308753,0.0,0.0,0.4662036850,0.4323657300,0.4662036850,0.4323657300,-0.93240737,-0.86473146] #<---- 8-shape
a=0
b=50
frames=30
println("Calculating...")
uis,tis=RKDP(fs,u0s,a,b)
println("Done.\nSaving...")
data=Dict("x1"=>uis[1],"y1"=>uis[2],
          "x2"=>uis[3],"y2"=>uis[4],
          "x3"=>uis[5],"y3"=>uis[6],
          "vx1"=>uis[7],"vy1"=>uis[8],
          "vx2"=>uis[9],"vy2"=>uis[10],
          "vx3"=>uis[11],"vy3"=>uis[12],
          "t"=>tis,)
open("3_body_problem.json","w") do file
    JSON.print(file,data,2)
end
println("Done.")