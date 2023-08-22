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
function generate_workspace_file(all_hists::DictT, filename, real_data; rebin2=true, systematics=false)::Dict where DictT
    all_hists::DictT = (rebin2 ? Dict(evt_type => Dict(var => all_hists[evt_type][var] |> restrict(120, Inf) |> rebin(2) for var in (systematics ? keys(all_hists[evt_type]) : [:HT_4j1b_nominal, :mbjj_4j2b_nominal])) for evt_type in keys(all_hists)) : copy(all_hists))

    region_names = Dict(
        :HT_4j1b => "4j1b CR",
        :mbjj_4j2b => "4j2b SR",
    )

    if systematics
        for evt_type in keys(all_hists)
            for x in keys(region_names)
                # symmetric reconstruction
                for var in [:_ME_var, :_PS_var, :_pt_res]
                    all_hists[evt_type][Symbol(x, var, :_sym)] = all_hists[evt_type][Symbol(x, :_nominal)]*3 - all_hists[evt_type][Symbol(x, var)]*2
                end
            end
        end
        modifiers_dict = Dict(
            #Symbol("Luminosity") => (:_lumi_up, :_lumi_down),
            Symbol("ME variation") => (:_ME_var, :_ME_var_sym),
            Symbol("PS variation") => (:_PS_var, :_PS_var_sym),
            Symbol("Jet energy resolution") => (:_pt_res, :_pt_res_sym),
            Symbol("Jet energy scale") => (:_pt_scale_up, :_pt_scale_down),
            Symbol("Scale variations") => (:_scale_var_up, :_scale_var_down),
            Symbol("b-tag NP 1") => (:_btag_var_0_up, :_btag_var_0_down),
            Symbol("b-tag NP 2") => (:_btag_var_1_up, :_btag_var_1_down),
            Symbol("b-tag NP 3") => (:_btag_var_2_up, :_btag_var_2_down),
            Symbol("b-tag NP 4") => (:_btag_var_3_up, :_btag_var_3_down),
        )
    end

    workspace = Dict(
        :channels => [
            Dict(
                :name => region_names[region],
                :samples => [
                    Dict(
                        :data => bincounts(all_hists[evt_type][Symbol(region, :_nominal)]),
                        :modifiers => cat( # systematics basically
                            Dict{Symbol, Any}[
                                Dict(
                                    :data => binerrors(all_hists[evt_type][Symbol(region, :_nominal)]),
                                    :name => "staterrors_"*replace(region_names[region], " " => "-"),
                                    :type => "staterror"
                                ),
                                Dict(
                                    :data => nothing,
                                    :name => "ttbar_norm",
                                    :type => "normfactor"
                                )
                            ],
                            reduce(append!, (systematics ? [
                                Dict{Symbol, Any}[Dict(
                                    :data => Dict(
                                        :hi => integral(all_hists[evt_type][Symbol(region, modifiers_dict[name][1])])/integral(all_hists[evt_type][Symbol(region, :_nominal)]),
                                        :lo => integral(all_hists[evt_type][Symbol(region, modifiers_dict[name][2])])/integral(all_hists[evt_type][Symbol(region, :_nominal)]),
                                    ),
                                    :name => name,
                                    :type => "normsys"
                                ),
                                Dict(
                                    :data => Dict(
                                        :hi_data => bincounts(all_hists[evt_type][Symbol(region, modifiers_dict[name][1])])*(integral(all_hists[evt_type][Symbol(region, :_nominal)])/integral(all_hists[evt_type][Symbol(region, modifiers_dict[name][1])])),
                                        :lo_data => bincounts(all_hists[evt_type][Symbol(region, modifiers_dict[name][2])])*(integral(all_hists[evt_type][Symbol(region, :_nominal)])/integral(all_hists[evt_type][Symbol(region, modifiers_dict[name][2])]))
                                    ),
                                    :name => name,
                                    :type => "histosys"
                                )] for name in keys(modifiers_dict)
                            ] : []), init=Dict{Symbol, Any}[]),
                        dims=1), 
                        :name => String(evt_type)
                    ) for evt_type in keys(all_hists)
                ]
            ) for region in keys(region_names)
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
