module LHC_AGC

using UnROOT, FHist, LorentzVectorHEP, JSON3
using LorentzVectorHEP: fromPxPyPzM
using Combinatorics: Combinations
using Distributions

include("constants.jl")
include("syst_utils.jl")
include("main_loop.jl")
include("visuals.jl")

function nevts_total(process_tag, variation=:nominal)
    NJSON[process_tag][variation][:nevts_total]
end

"""
    Convert xrd path from JSON to local path
"""
function xrd_to_local(url)
    joinpath(BASE_PATH[], last(split(url, '/')))
end

const TAG_PATH_DICT =
    Dict(
        k => Dict(
            var => xrd_to_local.(getindex.(NJSON[k][var][:files], :path)) for var in keys(NJSON[k])
        ) for k in LHC_AGC.TAGS
    )

"""
    download_data(N; process_tags=[:ttbar])

Download `N` files for each of the process tags.
"""
function download_data(N = MAX_N_FILES_PER_SAMPLE[]; process_tags = [:ttbar], variation_tags = [:nominal])
    for key in process_tags
        for var in variation_tags
            Threads.@threads for n in first(NJSON[key][var][:files], N)
                url = n[:path]
                fname = last(split(url, '/'))
                dir = BASE_PATH[]
                if !isdir(dir)
                    mkdir(dir)
                end
                out = joinpath(dir, fname)
                isfile(out) && continue
                download(url, out)
            end
        end
    end
    nothing
end

end
