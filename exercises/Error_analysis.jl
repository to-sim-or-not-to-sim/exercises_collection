using PyPlot
using PrettyTables

#=
 DISCLAIMER:
 BEING THESE MY FIRST EXERCISES WITH JULIA I DIDN'T HAVE CONFIDENCE
 WITH SYNTAX AND FUNCTIONS ALREADY IMPLEMENTED IN JULIA, SO SOME 
 EXERCISES MAY SEEM MORE COMPLICATED THAN THEY SHOULD BE; FOR EXAMPLE
 I USED "linspace(a,b,len,type)" INSTEAD OF "collect(a:step:b)".
=#

# USEFUL FUNCTIONS --------------------------------------------------

function linspace(a::Number,b::Number,len::Int;type::Type=Float64)
        x=zeros(type,len) 
        step=(b-a)/len
        for i in 1:len
            x[i]=a+i*step
        end
        return x
    end

#--------------------------------------------------------------------

# EX 1.1 ------------------------------------------------------------

function ex1()

    function my_exp(x::Number,approx::Int)
        result=1
        for i in 1:approx
            result+=(x^i)/factorial(i)
        end
        return result
    end

    function expected_delta(x::Number,approx::Int)
        return x^(approx+1)/factorial(approx+1)
    end

    approx=[1,2,3,4]
    fig,ax=plt.subplots(ncols=2,figsize=(11,5))
    colors=["red","forestgreen","darkcyan","blue"]
    x=linspace(0,1,50)
    ax[1].set_xlabel("x",fontsize=15)
    ax[2].set_xlabel("x",fontsize=15)
    ax[1].grid()
    ax[2].grid()
    ax[1].set_title("Value of \$|e^x-g(x,N)|\$",fontsize=20)
    ax[2].set_title("Value of \$\\frac{|\\Delta_{exp}-\\Delta_{found}|}{\\Delta_{exp}}\$",fontsize=15)
    for i in 1:4
        ax[1].plot(x,abs.(exp.(x).-my_exp.(x,approx[i])),color=colors[i],label="N = $i")
        ax[1].plot(x,expected_delta.(x,approx[i]),color=colors[i],linestyle="--")
        y=@. (abs(exp(x)-my_exp(x,approx[i]))-expected_delta(x,approx[i]))/expected_delta(x,approx[i])
        ax[2].plot(x,y,color=colors[i],label="N = $i")
    end
    ax[1].legend(fontsize=15)
    ax[2].legend(fontsize=15)
    savefig("result",dpi=500)
    show()

end

#--------------------------------------------------------------------

# EX 1.2 ------------------------------------------------------------

function ex2(type::Int64)
    elements=zeros(Float64,10^4)
    if type==32
        elements=zeros(Float32,10^4)
    end        

    len_el=length(elements)

    for i in 1:len_el
        elements[i]=1/(i^2)
    end

    function sum_one2N(N::Int)
        sum=0
        for i in 1:N
            sum+=elements[i]
        end
        return sum
    end

    function sum_N2one(N::Number)
        sum=0
        for i in N:-1:1
            sum+=elements[i]
        end
        return sum
    end

    function sim(from::Number,to::Number)
        expected_value=(pi^2)/6
        len=Int(to-from)
        N=linspace(from,to,len,type=Int64)
        sum1n=zeros(len)
        sumn1=zeros(len)
        sum_1n_old=sum_one2N(N[1])
        sum_n1_old=sum_N2one(N[1])
        sum1n[1]=sum_1n_old
        sumn1[1]=sum_n1_old
        found=false
        for i in 2:len
            sum1n_new=sum_one2N(N[i])
            sumn1_new=sum_N2one(N[i])
            if found==false
                if sum_1n_old==sum1n_new
                    found=true
                    println(N[i])
                end
            end
            sum1n[i]=sum1n_new
            sumn1[i]=sumn1_new
            sum_1n_old=sum1n_new
            sum_n1_old=sumn1_new
        end
        grid()
        title("Float$type",fontsize=15)
        plot(N,abs.(expected_value.-sum1n),label="\$\\Delta_{S_{1-N}}\$")
        plot(N,abs.(expected_value.-sumn1),label="\$\\Delta_{S_{N-1}}\$")
        xlabel("N",fontsize=15)
        legend(fontsize=15)
        savefig("N$type.png")
    show() 
    end

    sim(10^3,10^4)
end

#--------------------------------------------------------------------

# EX 1.3.1 ----------------------------------------------------------

function ex3(;type::Type=Float64)
    function mean(x::Vector{T} where {T<:Number})
        sum::type=0
        len=length(x)
        for i in 1:len
            sum+=x[i]
        end
        return sum/len
    end

    function variance(x::Vector{T} where {T<:Number})
        quad_sum::type=0
        len=length(x)
        m=mean(x)
        for i in 1:len
            quad_sum+=(m-x[i])^2
        end
        return quad_sum/(len-1)
    end

    function speed_variance(x::Vector{T} where {T<:Number})
        len=length(x)
        u::type=0
        v::type=0
        for i in 1:len
            u+=x[i]^2
            v+=x[i]
        end
        return (u-(v^2/len))/(len-1)
    end

    function final_variance(x::Vector{T} where {T<:Number})
        len=length(x)
        Mk=zeros(type,len)
        Qk=zeros(type,len)
        Mk[1]=x[1]
        Qk[1]=0
        for i in 2:len
            Mk[i]=Mk[i-1]+((x[i]-Mk[i-1])/i)
            Qk[i]=Qk[i-1]+(((i-1)*(x[i]-Mk[i-1])^2)/i)
        end
        return Qk[end]/(len-1)
    end

    xs= [[ 1e3, 1+1e3, 2+1e3 ],
        [ 1e6, 1+1e6, 2+1e6 ],
        [ 1e7, 1+1e7, 2+1e7 ],
        [ 1e8, 1+1e8, 2+1e8 ],
        rand(1:9, 25)]
    len=length(xs)
    for i in 1:len
        x=xs[i]
        println(x)
        println("Mean: ", mean(x),
                " Variance: ",variance(x),
                " Variance sped up: ",speed_variance(x),
                " Last variance: ", final_variance(x))
    end
end

#--------------------------------------------------------------------

# EX 1.3.2 ----------------------------------------------------------

function ex4()
    function naive(x::Number)
        return (exp(x)-1)/x
    end

    function k_f(x)
        return abs(((exp(x)*x)/(exp(x)-1))-1)
    end

    function taylor(x::Number,N::Int64)
        result=1
        den=1
        for i in 1:N
            den*=i+1
            result+=x^i/den
        end
        return result
    end

    function f_log(x::Number)
        y=exp(x)
        return (y-1)/log(y)
    end

    x=linspace(-1,1,1000)
    grid()
    title("\$\\kappa_f(x),\\ \\ \\ \\ f(x)=\\dfrac{e^x-1}{x}\$",fontsize=15)
    plot(x,k_f.(x),color="navy")
    xlabel("x",fontsize=15)
    ylabel("\$\\kappa_f(x)\$",fontsize=15)
    savefig("k_f.png",dpi=200)
    show()

    xs=[1/10^i for i in 3:16]
    data=hcat(xs, naive.(xs), taylor.(xs,5), f_log.(xs))
    column_labels=["x","f(x)","f(x,n=5)","g(x)"]
    style = TextTableStyle(first_line_column_label = crayon"yellow bold")
    table_format = TextTableFormat(borders = text_table_borders__unicode_rounded)
    hg=TextHighlighter(
            (data, i, j) -> (j == 1),
            crayon"green bold"
        )
    hr=TextHighlighter(
            (data,i,j) -> (i==14),
            crayon"fg:white bold bg:red bold"
        )

    pretty_table(
        data;
        column_labels = column_labels,
        style         = style,
        highlighters  = [hr,hg],
        table_format  = table_format_format  = TextTableFormat(borders = text_table_borders__unicode_rounded)
    )
end

#--------------------------------------------------------------------

# OPTIONAL 1 --------------------------------------------------------

function ex1_optional()
    function x_k(k::Int64)
        xs=zeros(Float64,k+1)
        for i in 0:k
            if i==0
                x=11.0/2.0
            elseif i==1
                x=61.0/11.0
            else
                x=111-((1130-(3000/(xs[i-1])))/xs[i])
            end
            println("$i: $x")
            xs[i+1]=x
        end
        return xs
    end

    vals=x_k(34)
end

#--------------------------------------------------------------------

# OPTIONAL 2 --------------------------------------------------------

function ex2_optional()
    e_minus_one=1.7182818284590452
    for j in 1:1:16
        start_value::Float64=round(e_minus_one,digits=j)
        time=Int32[]
        money=Float64[]
        push!(time,0)
        push!(money,start_value)
        for i in 1:25
            new_val=i*money[i]-1
            if new_val<0
                break
            end
            push!(time,i)
            push!(money,new_val)
        end
        plot(time,money,color="navy",label="e rounded to $j digit(s)")
        grid()
        yscale("log")
        scatter(time,money,color="navy",s=30)	
        legend(fontsize=15)
        savefig("graph$j.png",dpi=200)
        show()	
    end
end

#--------------------------------------------------------------------




# EXECUTION ---------------------------------------------------------

function main()
    println("--- MAIN ---")
    println("Exercise 1.1")
    ex1()
    println("-"^30)
    println("Exercise 1.2")
    ex2(32)
    ex2(64)
    println("-"^30)
    println("Exercise 1.3.1")
    println("--- FLOAT32 ---")
    ex3(type=Float32)
    println("--- FLOAT64 ---")
    ex3(type=Float64)
    println("-"^30)
    println("Exercise 1.3.1")
    ex4()
    println("-"^30)
end

function optional()
    println("--- OPTIONAL ---")
    println("Exercise A.1.1")
    ex1_optional()
    println("-"^30)
    println("Exercise A.1.2")
    ex2_optional()
    println("-"^30)
end

#=
 IT IS POSSIBLE TO ANALYZE SINGLE EXERCISES BY
 REPLACING THE NEXT FEW LINE WITH THE FUNCTION
 RELATED TO EXERCISE YOU WANT TO ANALYZE.
=#

main()
optional()