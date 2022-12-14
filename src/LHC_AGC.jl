module LHC_AGC

using UnROOT, FHist, LorentzVectorHEP, JSON3
using LorentzVectorHEP: fromPxPyPzM
using Combinatorics: Combinations
using Distributions

include("constants.jl")
include("main_loop.jl")
include("syst_utils.jl")

function nevts_total(tag, variation=:nominal)
    NJSON[tag][variation][:nevts_total]
end

"""
    Convert xrd path from JSON to local path
"""
function xrd_to_local(url)
    joinpath(BASE_PATH[], last(split(url, '/')))
end

const TAG_PATH_DICT = Dict(k=>xrd_to_local.(getindex.(NJSON[k][:nominal][:files], :path))
       for k in LHC_AGC.TAGS
      )

function download_data(N = MAX_N_FILES_PER_SAMPLE[])
    for key in keys(NJSON)
        Threads.@threads for n in first(NJSON[key][:nominal][:files], N)
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
    nothing
end

end
