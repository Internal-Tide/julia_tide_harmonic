include("t_getconsts.jl")
using LinearAlgebra
using Dates
using SparseArrays

ltype = 1
ctime = DateTime(2015, 1, 1, 0, 0, 0)
lat = 5
ju = 1


function t_vuf(ltype, ctime, ju, lat=nothing)

    astro, ader = t_astron(ctime)

    consts, sat, shallow = t_getconsts(ctime)
    v = mod.(consts["doodson"]' * astro .+ consts["semi"], 1)

    if lat !== nothing
        if abs(lat) < 5
            lat = sign(lat) * 5
        end
        slat = sind(lat)

        rr = sat["amprat"]
        if isfinite(lat)
            j = findall(x -> x == 1, sat["ilatfac"])
            rr[j] .= rr[j] .* 0.36309 .* (1.0 .- 5.0 .* slat .^ 2) ./ slat

            j = findall(x -> x == 2, sat["ilatfac"])
            rr[j] .= rr[j] .* 2.59808 .* slat
        else
            rr[sat["ilatfac"].>0] .= 0
        end
        uu = mod.(sat["deldood"]' * astro[4:6, :] .+ sat["phcorr"], 1)
        nsat = maximum(size(sat["iconst"]))
        nfreq = maximum(size(consts["isat"]))

        fsum = 1 .+ sum(spdiagm(0 => vec(rr .* exp.(2im * pi * uu)))[:, vec(sat["iconst"])], dims=1)

        f = abs.(fsum)
        u = angle.(fsum) / (2 * pi)

        shallow_m1 = consts["ishallow"] .- 1
        iname_m1 = shallow["iname"] .- 1
        coefs = shallow["coef"]
        range_cache = [0:consts["nshallow"][k]-1 for k in 1:length(consts["nshallow"])]
        for k in findall(isfinite.(consts["ishallow"]))
            ik = Int(shallow_m1[k]) .+ range_cache[consts["nshallow"][k]]
            iname = iname_m1[ik.+1]
            @info iname
            coef = coefs[k]
            @info size(f)
            f[k] .= prod(f[iname+1] .^ coef)
            u[k] .= dot(u[iname+1], coef)
            v[k] .= dot(v[iname+1], coef)
        end

        f .= f[ju.+1]
        u .= u[ju.+1]
        v .= v[ju.+1]

    else

        for k in findall(isfinite.(consts["ishallow"]))
            ik = consts["ishallow"][k] - 1 + (0:consts["nshallow"][k]-1)
            v[k] .= dot(v[shallow["iname"][ik]], shallow["coef"][ik])
        end

        v .= v[ju.+1]
        f .= ones(length(v))
        u .= zeros(length(v))
    end


    return u, v, f
end


t_vuf(ltype, ctime, ju, lat)




findall(isfinite.(_const["ishallow"]))