const xsec_info = Dict(
    :ttbar => 396.87 + 332.97, # nonallhad + allhad, keep same x-sec for all
    :single_top_s_chan => 2.0268 + 1.2676,
    :single_top_t_chan => (36.993 + 22.175)/0.252,  # scale from lepton filter to inclusive
    :single_top_tW => 37.936 + 37.906,
    :wjets => 61457 * 0.252,  # e/mu+nu final states
    :data => 1.0
)
const NJSON = joinpath(dirname(@__DIR__), "nanoaod_inputs_small.json") |> read |> JSON3.read;
const TAGS = keys(NJSON)

const BASE_PATH = Ref(joinpath(dirname(@__DIR__), "data"))
const MAX_N_FILES_PER_SAMPLE = Ref(1)
const LUMI = 3378.0 # /pb
