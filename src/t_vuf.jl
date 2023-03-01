include("t_getconsts.jl")

using LinearAlgebra


function t_vuf(ltype, ctime, ju, lat=nothing)
    """T_VUF Computes nodal modulation corrections.
     [V,U,F]=T_VUF(TYPE,DATE,JU,LAT) returns the astronomical phase V, the
     nodal phase modulation U, and the nodal amplitude correction F at
     a decimal date DATE for the components specified by index JU
     at a latitude LAT.

     TYPE is either 'full' for the 18.6 year set of constitunets, or 'nodal'
     for the 1-year set with satellite modulations.

     If LAT is not specified, then the Greenwich phase V is computed with
     U=0 and F=1.

     Note that V and U are in 'cycles', not degrees or radians (i.e.,
     multiply by 360 to get degrees).

     If LAT is set to NaN, then the nodal corrections are computed for all
     satellites that do *not* have a "latitude-dependent" correction
     factor. This is for compatibility with the ways things are done in
     the xtide package. (The latitude-dependent corrections were zeroed
     out there partly because it was convenient, but this was rationalized
     by saying that since the forcing of tides can occur at latitudes
     other than where they are observed, the idea that observations have
     the equilibrium latitude-dependence is possibly bogus anyway).
     Get all the info about constituents.
     Calculate astronomical arguments at mid-point of data time series.
    """
    astro, ader = t_astron(ctime)

    if ltype == "full"
        consts = t_get18consts(ctime)
        v = mod.(consts.doodson * astro .+ consts.semi, 1)
        v = v[ju]
        u = zeros(size(v))
        f = ones(size(v))
    else
        consts, sat, shallow = t_getconsts(ctime)
        v = mod.(consts["doodson"] * astro .+ consts["semi"], 1)

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
        end
    end

    return u, v, f
end