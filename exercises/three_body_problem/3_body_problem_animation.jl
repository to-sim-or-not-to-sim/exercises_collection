using JSON,Plots
pyplot()
data=JSON.parsefile("3_body_problem.json")
x1=data["x1"]
y1=data["y1"]
x2=data["x2"]
y2=data["y2"]
x3=data["x3"]
y3=data["y3"]
tis=data["t"]

len=length(tis)
start=1
anim=@animate for i in 51:len
    t_round=round(tis[i],digits=2)
    scatter([x1[i]],[y1[i]],color=:deepskyblue,title="t=$t_round",titlefontsize=15,legend=false)
    scatter!([x2[i]],[y2[i]],color=:crimson,legend=false)
    scatter!([x3[i]],[y3[i]],color=:goldenrod,legend=false)
    plot!(x1[start:i],y1[start:i],color=:deepskyblue,legend=false)
    plot!(x2[start:i],y2[start:i],color=:crimson,legend=false)
    plot!(x3[start:i],y3[start:i],color=:goldenrod,legend=false)
    start+=1
end
gif(anim,"three_body_problem.gif",fps=30)