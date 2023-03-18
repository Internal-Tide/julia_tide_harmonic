include("t_astron.jl")

function fourpad(conin)
    conin = map(x -> string(x), conin)
    for i in eachindex(conin)
        conin[i] = rpad(conin[i], 4)
    end
    conin
end

