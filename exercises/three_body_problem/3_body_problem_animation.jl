using JSON,GLMakie

data=JSON.parsefile("./fisica computazionale/0/3bp_1_unbound/3_body_problem.json")
x1=data["x1"]
y1=data["y1"]
x2=data["x2"]
y2=data["y2"]
x3=data["x3"]
y3=data["y3"]
tis=data["t"]
len = length(tis)

fig=Figure(resolution=(600,600))
ax=Axis(fig[1,1])

pts1=Observable(Point2f[])
pts2=Observable(Point2f[])
pts3=Observable(Point2f[])

trail1=Observable(Point2f[])
trail2=Observable(Point2f[])
trail3=Observable(Point2f[])

alpha_obs=Observable(Float32[])

scatter!(ax, pts1, color=:deepskyblue, markersize=10)
scatter!(ax, pts2, color=:crimson, markersize=10)
scatter!(ax, pts3, color=:goldenrod, markersize=10)
lines!(ax, trail1, color=:deepskyblue, transparency=true)
lines!(ax, trail2, color=:crimson, transparency=true)
lines!(ax, trail3, color=:goldenrod, transparency=true)

title_obs=Observable("t=0.0")
ax.title=title_obs

record(fig,"three_body_problem.gif",1:10:len;framerate=30) do i
    println(round(i/len*100, digits=2), " %")
    title_obs[]="t=$(round(tis[i], digits=2))"
    pts1[]=[Point2f(x1[i],y1[i])]
    pts2[]=[Point2f(x2[i],y2[i])]
    pts3[]=[Point2f(x3[i],y3[i])]
    a=max(1,i-200)
    trail1[]=Point2f.(x1[a:i],y1[a:i])
    trail2[]=Point2f.(x2[a:i],y2[a:i])
    trail3[]=Point2f.(x3[a:i],y3[a:i])
end