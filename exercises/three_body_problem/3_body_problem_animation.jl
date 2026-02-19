using JSON,GLMakie

data=JSON.parsefile("3_body_problem.json")
x1=data["x1"]
y1=data["y1"]
x2=data["x2"]
y2=data["y2"]
x3=data["x3"]
y3=data["y3"]
tis=data["t"]
len=length(tis)

fig=Figure(size=(600,600))
ax=Axis(fig[1,1])

pts1=Observable(Point2f[])
pts2=Observable(Point2f[])
pts3=Observable(Point2f[])

trail1=Observable(Point2f[])
trail2=Observable(Point2f[])
trail3=Observable(Point2f[])

blue=RGBf(0,0,1)
green=RGBf(0,1,0)
red=RGBf(1,0,0)
fading_blue=Observable(RGBAf[])
fading_green=Observable(RGBAf[])
fading_red=Observable(RGBAf[])

scatter!(ax,pts1,color=blue,markersize=10)
scatter!(ax,pts2,color=green,markersize=10)
scatter!(ax,pts3,color=red,markersize=10)
lines!(ax,trail1,color=fading_blue)
lines!(ax,trail2,color=fading_green)
lines!(ax,trail3,color=fading_red)
xlims!(ax,minimum(vcat(x1,x2,x3)),maximum(vcat(x1,x2,x3)))
ylims!(ax,minimum(vcat(y1,y2,y3)),maximum(vcat(y1,y2,y3)))

record(fig,"three_body_problem.gif",1:10:len;framerate=30) do i  
    a=max(1,i-200)
    fading=@. Float32(1/(1+exp(-((a:i)-i+50)/25)))
    fading_blue[]=RGBAf.(blue.r,blue.g,blue.b,fading)
    fading_green[]=RGBAf.(green.r,green.g,green.b,fading)
    fading_red[]=RGBAf.(red.r,red.g,red.b,fading)
    ax.title[]="t=$(round(tis[i], digits=2))"
    pts1[]=[Point2f(x1[i],y1[i])]
    pts2[]=[Point2f(x2[i],y2[i])]
    pts3[]=[Point2f(x3[i],y3[i])]
    trail1[]=Point2f.(x1[a:i],y1[a:i])
    trail2[]=Point2f.(x2[a:i],y2[a:i])
    trail3[]=Point2f.(x3[a:i],y3[a:i])
    println(round(i/len*100, digits=2), " %")
end


