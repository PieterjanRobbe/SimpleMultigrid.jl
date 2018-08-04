# test_utils.jl : utilities for testing SimpleMultigrid

# log resnorm
function log(mg::MultigridIterable)
    println(string("+","-"^30,"+"))
    println(string("|"," "^11,@sprintf("n = %4i",sqrt(size(mg.grids[1].A,1))+1)," "^11,"|"))
    println(string("+","-"^9,"+","-"^20,"+"))
    println(string("| V-cycle |  |rʰ|"," "^7,"ratio","  |"))
    println(string("+","-"^9,"+","-"^20,"+"))
    println(string(@sprintf("|   %4i  |",0),@sprintf("  %6.3e",mg.resnorm[1])," "^8," |"))
    for i in 2:length(mg.resnorm)
        println(string(@sprintf("|   %4i  |",i-1),@sprintf("  %6.3e",mg.resnorm[i]),@sprintf("  %3.2f",mg.resnorm[i]/mg.resnorm[i-1]),"   |"))
    end
    println(string("+","-"^9,"+","-"^20,"+"))
end

# plot 
@require PyPlot begin
    function mplot(x...)
        n = maximum(length.(x))+1
        figure(figsize=(5,5))
        for i in 1:length(x)
            factor = round(Int64,n/(length(x[i])+1))
            plot(factor:factor:n-1,x[i],"-o")
        end
    end

    function msurf(s...)
        n = maximum(size.(s,1))+1
        m = maximum(size.(s,2))+1
        figure(figsize=(5,5))
        alphas = linspace(0.2,1,length(s))
        for i in 1:length(s)
            nfactor = round(Int64,n/(size(s[i],1)+1))
            mfactor = round(Int64,m/(size(s[i],2)+1))
            x = nfactor:nfactor:n-1
            y = mfactor:mfactor:m-1
            xgrid = repmat(x,1,length(y))
            ygrid = repmat(y',length(x),1)
            surf(xgrid,ygrid,s[i],alpha=alphas[i])
        end
    end
end
