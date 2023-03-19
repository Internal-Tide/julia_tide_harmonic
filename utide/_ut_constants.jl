# ut_constants["const"]["nsat"]
using MAT
_base_dir = joinpath(dirname(@__FILE__), "..", "data")
_ut_constants_fname = joinpath(_base_dir, "ut_constants.mat")
ut_constants = matread(_ut_constants_fname)
constit_names = ut_constants["const"]["name"]

# Make a dictionary for index lookups.
constit_index_dict = Dict(name => i for (i, name) in enumerate(constit_names))
_uc = ut_constants["const"]
cycles_per_hour = Dict(zip(_uc["name"], _uc["freq"]))
hours_per_cycle = Dict(zip(_uc["name"], 1 ./ _uc["freq"]))

