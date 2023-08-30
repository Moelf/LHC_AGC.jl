"""
    get_all_hists(; do_file_variations::Bool=true, wgt = 0.0, n_files_max_per_sample = MAX_N_FILES_PER_SAMPLE[], tags = LHC_AGC.TAGS, histo_getter=get_histo)

Produces the `all_hists` dictionary that may be required for building `workspace.json` and plotting.

`histo_getter` should either be `get_histo` or `get_histo_distributed`.
"""
function get_all_hists(; do_file_variations::Bool=true, wgt = 0.0, n_files_max_per_sample = MAX_N_FILES_PER_SAMPLE[], tags = LHC_AGC.TAGS, histo_getter=get_histo)
    Dict(tag => histo_getter(tag; do_file_variations, wgt, n_files_max_per_sample) for tag in tags)
end

"""
    get_histo(process_tag::Symbol; do_file_variations::Bool=true, wgt = 0.0, n_files_max_per_sample = MAX_N_FILES_PER_SAMPLE[])
"""
function get_histo(process_tag::Symbol; do_file_variations::Bool=true, wgt = 0.0, n_files_max_per_sample = MAX_N_FILES_PER_SAMPLE[])
    N = n_files_max_per_sample
    if iszero(wgt)
        wgt = LUMI * xsec_info[process_tag] / nevts_total(process_tag)
    end

    file_variation_tags = (do_file_variations ? keys(TAG_PATH_DICT[process_tag]) : [:nominal])

    all_hists = reduce(merge, [
        mapreduce(mergewith(+), @view TAG_PATH_DICT[process_tag][variation_tag][begin:N]) do path
            get_histo(LazyTree(path, "Events"), wgt, file_variation=variation_tag)
        end for variation_tag in file_variation_tags
    ])
    all_hists
end

"""
    get_histo_distributed(process_tag::Symbol; do_file_variations::Bool=true, wgt = 0.0, n_files_max_per_sample = MAX_N_FILES_PER_SAMPLE[])
"""
function get_histo_distributed(process_tag::Symbol; do_file_variations::Bool=true, wgt = 0.0, n_files_max_per_sample = MAX_N_FILES_PER_SAMPLE[])
    N = n_files_max_per_sample
    if iszero(wgt)
        wgt = LUMI * xsec_info[process_tag] / nevts_total(process_tag)
    end

    file_variation_tags = (do_file_variations ? keys(TAG_PATH_DICT[process_tag]) : [:nominal])

    files = Tuple{String, Symbol}[]
    for variation_tag in file_variation_tags
        append!(files, Tuple{String, Symbol}[(path, variation_tag) for path in @view TAG_PATH_DICT[process_tag][variation_tag][begin:N]])
    end

    mainloop = function (tuple)
        path, variation_tag = tuple
        get_histo(LazyTree(path, "Events"), wgt,file_variation=variation_tag)
    end
    dicts = pmap(mainloop, files)

    return reduce(mergewith(+), dicts)
end

function generate_hists(file_variation::Symbol)
    nedges=26 #25 bins
    start=50
    stop=550

    all_keys = (file_variation,)
    if file_variation == :nominal
        scale_base = keys(SCALE_VARS)
        scale_names = map(splat(Symbol), Iterators.product(scale_base, (:_up, :_down)))
        all_keys = (keys(SHAPE_VARS)..., scale_names...)
    end

    hists = merge(
        Dict(
            Symbol(:mbjj_4j2b_, k) => Hist1D(Float64; bins = range(; start, stop, length=nedges)) for k in all_keys
        ),
        Dict(
            Symbol(:HT_4j1b_, k) => Hist1D(Float64; bins = range(; start, stop, length=nedges)) for k in all_keys
        )
    )
    return hists
end

"""
    get_histo(tree, wgt; file_variation::Symbol=:nominal, evts=nothing)

    `evts` is used to track the events processed for each histogram type and should be a dictionary of the format histogram_type => Vector{Int}. the dictionary gets mutated and is not returned.
"""
function get_histo(tree, wgt; file_variation::Symbol=:nominal, evts=nothing)
    is_nominal_file = (:nominal == file_variation)
    hists = generate_hists(file_variation)
    Threads.@threads for evt in tree
        # single lepton requirement
        (; Electron_pt, Muon_pt) = evt
        (count(>(25), Electron_pt) + count(>(25), Muon_pt) != 1) && continue

        # get pt
        (; Jet_pt) = evt
        is_nominal_file && (Jet_pt_nominal = Jet_pt)

        for hist_type in (is_nominal_file ? keys(SHAPE_VARS) : (:nominal,))
            # modify pt
            is_nominal_file && (Jet_pt = SHAPE_VARS[hist_type].(Jet_pt_nominal))

            scale_info = (; Jet_pt, wgt)

            # at least 4 jets
            jet_pt_mask = Jet_pt .> 25
            if count(jet_pt_mask) >= 4
                jet_btag = @view evt.Jet_btagCSVV2[jet_pt_mask]

                btag_count = count(>(0.5), jet_btag)
                # MASS HISTOGRAM
                if btag_count >= 2 # at least 2 btag
                    if evts !== nothing
                        if hist_type in keys(evts)
                            push!(evts[hist_type], evt.event)
                        end
                    end

                    (; Jet_eta, Jet_phi, Jet_mass) = evt

                    # construct jet lorentz vector
                    jet_p4 = @views LorentzVectorCyl.(Jet_pt[jet_pt_mask], Jet_eta[jet_pt_mask], Jet_phi[jet_pt_mask], Jet_mass[jet_pt_mask])
                    #jet_p4 = @views fromPxPyPzM.(Jet_px[jet_pt_mask], Jet_py[jet_pt_mask], Jet_pz[jet_pt_mask], Jet_mass[jet_pt_mask])

                    Njets = length(jet_btag) 
                    # Njets == length(jet_p4) || error("impossible reached")

                    # tri jet combinatorics
                    max_pt = -Inf
                    local best_mass
                    
                    # 1. all tri-jet combinations
                    for comb in Combinations(Njets, 3)
                        p4s = @view jet_p4[comb]
                        btags = @view jet_btag[comb]
                        # 2. keep those maximum(btags1,2,3) > 0.5
                        maximum(btags) <= 0.5 && continue
                        tri = sum(p4s)
                        _pt = pt(tri)
                        # 3. pick the tri-p4 with highest tri-pt
                        if _pt > max_pt
                            max_pt = _pt
                            best_mass = mass(tri)
                        end
                    end

                    # tri-p4 with highest tri-pt first
                    push!(hists[Symbol(:mbjj_4j2b_, (is_nominal_file ? hist_type : file_variation))], best_mass, wgt)
                    if is_nominal_file && (hist_type == :nominal)
                        @scale_var_loop :mbjj_4j2b best_mass
                    end
                # HT HISTOGRAM
                elseif btag_count == 1 # no more than 1 btag
                    HT = @views sum(Jet_pt[jet_pt_mask])
                    push!(hists[Symbol(:HT_4j1b_, (is_nominal_file ? hist_type : file_variation))], HT, wgt)
                    if is_nominal_file && (hist_type == :nominal)
                        @scale_var_loop :HT_4j1b HT
                    end
                end
            end
        end
    end

    hists
end
