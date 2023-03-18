using NCDatasets

include("t_astron.jl")


# set the base directory for the data
_base_dir = joinpath(dirname(@__FILE__), "..", "data")
const has_const = isfile(joinpath(_base_dir, "t_constituents_const.nc"))
const has_sat = isfile(joinpath(_base_dir, "t_constituents_sat.nc"))
const has_shallow = isfile(joinpath(_base_dir, "t_constituents_shallow.nc"))
if (has_const & has_sat & has_shallow)
    _const = Dict()
    _sat = Dict()
    _shallow = Dict()

    constituents_const = Dataset(joinpath(_base_dir, "t_constituents_const.nc"))
    for key in keys(constituents_const)
        _const[key] = constituents_const[key][:]
    end

    constituents_sat = Dataset(joinpath(_base_dir, "t_constituents_sat.nc"))
    for key in keys(constituents_sat)
        _sat[key] = constituents_sat[key][:]
    end

    constituents_shallow = Dataset(joinpath(_base_dir, "t_constituents_shallow.nc"))
    for key in keys(constituents_shallow)
        _shallow[key] = constituents_shallow[key][:]
    end

    close(constituents_const)
    close(constituents_sat)
    close(constituents_shallow)

else
    println("You do not have t_constituents_*.nc "$
    "check that package installation is correct.")
    _const = Dict()
    _sat = Dict()
    _shallow = Dict()
end



function t_getconsts(ctime)
    """
    t_getconsts - Gets constituent data structures holding
                  information for tidal analyses
    Variables are loaded from 't_constituents_*.npy'
    on init and a copy is made now.
    When ctime is specified t_getconsts recomputes the frequencies from
    the rates-of-change of astronomical parameters at the matlab TIME given.

     :Parameters:
        ctime: a datetime, the start time of the data input into t_tide.

    Note:
        Not sure if a copy has to be made here or if they can be used directly.
        For now a copy should be much fast then a load.
    """
    consts = deepcopy(_const)
    sat = deepcopy(_sat)
    shallow = deepcopy(_shallow)

    if ctime !== nothing
        # If no time, just take the "standard" frequencies, otherwise
        # compute them from derivatives of astro parameters. This is
        # probably a real overkill - the diffs are in the
        # 10th decimal place (9th sig fig).
        astro, ader = t_astron(ctime)

        ii = .!isnan.(consts["ishallow"])
        # println(size(ii))
        # println(size(consts["doodson"]))
        # println(size(ader))
        consts["freq"][.!ii] .= (consts["doodson"][:, .!ii]' * ader) / 24

        shallow_m1 = consts["ishallow"] .- 1
        iname_m1 = shallow["iname"] .- 1
        range_cache = Dict(n => collect(0:n-1) for n in 0:maximum(consts["nshallow"]))
        # for k in findall(ii)
        #     ik = Int(shallow_m1[k]) .+ range_cache[consts["nshallow"][k]]
        #     # println(typeof(shallow_m1))
        #     # println(ik)
        #     # println(size(shallow_m1))
        #     # println(size(iname_m1))
        #     # println(range_cache)
        #     consts["freq"][k] .= dot(consts["freq"][iname_m1[ik]], shallow["coef"][ik])
        # end

        for k in findall(x -> x != 0, ii)
            println(shallow_m1)
            ik = Int(shallow_m1[k]) .+ range_cache[consts["nshallow"][k]]
            println(ik)
            consts["freq"][k] = sum(consts["freq"][iname_m1[ik.+1]] .* shallow["coef"][ik.+1])
        end
    end

    return consts, sat, shallow
end

# _const["freq"][30]

# iname_m1 = _shallow["iname"] .- 1
# _const["freq"][]
# shallow_m1 = _const["ishallow"] .- 1
# iname_m1 = _shallow["iname"] .- 1
# range_cache = Dict(n => collect(0:n-1) for n in 0:maximum(_const["nshallow"]))
# ik = Int(shallow_m1[30]) .+ range_cache[_const["nshallow"][30]]
# _const["freq"][iname_m1[ik.+1]]
# _shallow["coef"][ik.+1]