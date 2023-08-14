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

"""
    generate_workspace_file(all_hists::Dict, filename, real_data; rebin2=true, systematics=false)::Dict


Generates a workspace dictionary, writes it to a JSON file and returns


`all_hists` should be a dictionary of the form evt_type => hists_dict, where hists_dict is what get_histo(evt_type, ...) would return.

`real_data` should contain real data that we are going to put into `:observations` (essentially vectors of numbers). we put `real_data[1]` for `"4j1b CR"` and `real_data[2]` for `"4j2b SR"`.
"""
function generate_workspace_file(all_hists::Dict, filename, real_data; rebin2=true, systematics=false)::Dict
    all_hists = (rebin2 ? Dict(evt_type => Dict(var => all_hists[evt_type][var] |> restrict(120, Inf) |> rebin(2) for var in (systematics ? keys(all_hists[evt_type]) : [:HT_4j1b_nominal, :mbjj_4j2b_nominal])) for evt_type in keys(all_hists)) : all_hists)

    workspace = Dict(
        :channels => [
            Dict(
                :name => "4j1b CR",
                :samples => [
                    Dict(
                        :data => bincounts(all_hists[evt_type][:HT_4j1b_nominal]), 
                        :modifiers => [ # systematics basically
                            Dict(
                                :data => binerrors(all_hists[evt_type][:HT_4j1b_nominal]),
                                :name => "staterrors_4j1b-CR",
                                :type => "staterror"
                            ),
                            Dict(
                                :data => nothing,
                                :name => "ttbar_norm",
                                :type => "normfactor"
                            )
                        ], 
                        :name => String(evt_type),
                    ) for evt_type in keys(all_hists)
                ]
            ),
            Dict(
                :name => "4j2b SR",
                :samples => [
                    Dict(
                        :data => bincounts(all_hists[evt_type][:mbjj_4j2b_nominal]),
                        :modifiers => [ # systematics basically
                            Dict(
                                :data => binerrors(all_hists[evt_type][:mbjj_4j2b_nominal]),
                                :name => "staterrors_4j2b-SR",
                                :type => "staterror"
                            )
                        ], 
                        :name => String(evt_type)
                    ) for evt_type in keys(all_hists)
                ]
            )
        ],
        :measurements => [ # we leave it so for now
            Dict(
                :config => Dict(
                    :parameters => [
                        Dict(
                            :bounds => [
                                [0, 10] # allowed values
                            ],
                            :inits => [
                                1.0 # initial value
                            ],
                            :name => "ttbar_norm"
                        )
                    ],
                    :poi => "ttbar_norm"
                ),
                :name => "CMS_ttbar"
            )
        ],
        :observations => [
            Dict(
                :name => "4j1b CR",
                :data => real_data[1]
            ),
            Dict(
                :name => "4j2b SR",
                :data => real_data[2]
            )
        ],
        :version => "1.0.0"
    )

    open(filename, "w") do f
        JSON3.pretty(f, workspace)
    end

    workspace
end

end
