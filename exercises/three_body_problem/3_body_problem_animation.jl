using JSON,Plots
pyplot()
data=JSON.parsefile("./fisica computazionale/0/3bp_1_unbound/3_body_problem.json")
x1=data["x1"]
y1=data["y1"]
x2=data["x2"]
y2=data["y2"]
x3=data["x3"]
y3=data["y3"]
tis=data["t"]

len=length(tis)
anim=@animate for i in 1:10:len
    t_round=round(tis[i],digits=2)
    scatter([x1[i]],[y1[i]],color=:deepskyblue,title="t=$t_round",titlefontsize=15,legend=false)
    scatter!([x2[i]],[y2[i]],color=:crimson,legend=false)
    scatter!([x3[i]],[y3[i]],color=:goldenrod,legend=false)
    x=collect(1:i)
    alpha_val=@. 1/(1+exp(-(x-i+50)/25))
    plot!(x1[1:i],y1[1:i],alpha=alpha_val,color=:deepskyblue,legend=false)
    plot!(x2[1:i],y2[1:i],alpha=alpha_val,color=:crimson,legend=false)
    plot!(x3[1:i],y3[1:i],alpha=alpha_val,color=:goldenrod,legend=false)
end
gif(anim,"three_body_problem.gif",fps=30)