using JSON,Plots
pyplot()
data=JSON.parsefile("3_body_problem.json")
x1=data["x1"]
y1=data["y1"]
x2=data["x2"]
y2=data["y2"]
x3=data["x3"]
y3=data["y3"]
vx1=data["vx1"]
vy1=data["vy1"]
vx2=data["vx2"]
vy2=data["vy2"]
vx3=data["vx3"]
vy3=data["vy3"]
tis=data["t"]

function Energy(x1,y1,x2,y2,x3,y3,vx1,vy1,vx2,vy2,vx3,vy3)
    K=(vx1^2+vy1^2+vx2^2+vy2^2+vx3^2+vy3^2)/2
    U=-1/sqrt((x1-x2)^2+(y1-y2)^2)-1/sqrt((x1-x3)^2+(y1-y3)^2)-1/sqrt((x3-x2)^2+(y3-y2)^2)
    return K+U
end

len=length(x1)
step=div(len,20)
for i in 1:step:len
    E=Energy(x1[i],y1[i],x2[i],y2[i],x3[i],y3[i],vx1[i],vy1[i],vx2[i],vy2[i],vx3[i],vy3[i])
    println("E = ",E)
end
E_start=Energy(x1[1],y1[1],x2[1],y2[1],x3[1],y3[1],vx1[1],vy1[1],vx2[1],vy2[1],vx3[1],vy3[1])
E_end=Energy(x1[end],y1[end],x2[end],y2[end],x3[end],y3[end],vx1[end],vy1[end],vx2[end],vy2[end],vx3[end],vy3[end])
println("E_i = ",E_start)
println("E_f = ",E_end)