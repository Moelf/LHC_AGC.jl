"""
    get_histo(process_tag::Symbol; wgt = 0.0, n_files_max_per_sample = MAX_N_FILES_PER_SAMPLE[])
"""
function get_histo(process_tag::Symbol; variation_tag::Symbol=:nominal, wgt = 0.0, n_files_max_per_sample = MAX_N_FILES_PER_SAMPLE[])
    N = n_files_max_per_sample
    if iszero(wgt)
        wgt = LUMI * xsec_info[process_tag] / nevts_total(process_tag)
    end
    fs = @view TAG_PATH_DICT[process_tag][variation_tag][begin:N]
    println(fs)
    hists = mapreduce(mergewith(mergewith(+)), fs) do path #
        println(path)
        tree = LazyTree(path, "Events")
        println(length(tree))
        get_histo(tree, wgt)
    end
    hists
end

function generate_hists()
    nbins=26
    start=50
    stop=550

    scale_base = keys(SCALE_VARS)
    scale_names = map(splat(Symbol), Iterators.product(scale_base, (:_up, :_down)))
    hists = Dict(k => Dict(
            "4j2b" => Hist1D(Float64; bins = range(; start, stop, length=nbins)),
            "4j1b" => Hist1D(Float64; bins = range(; start, stop, length=nbins))
            )
            for k in (keys(SHAPE_VARS)..., scale_names...)
        )
    return hists
end


"""
    get_histo(tree, wgt; evts=nothing)

    `evts` is used to track the events processed for each histogram type and should be a dictionary of the format histogram_type => Vector{Int}. the dictionary gets mutated and is not returned.
"""
function get_histo(tree, wgt; evts=nothing)
    hists = generate_hists()
    for evt in tree
        # single lepton requirement
        (; Electron_pt, Muon_pt) = evt
        (count(>(25), Electron_pt) + count(>(25), Muon_pt) != 1) && continue

        # get pt
        (; Jet_pt) = evt
        Jet_pt_nominal = Jet_pt

        for hist_type in keys(SHAPE_VARS)
            # modify pt
            Jet_pt = SHAPE_VARS[hist_type].(Jet_pt_nominal)

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
                    push!(hists[hist_type]["4j2b"], best_mass, wgt)
                    if hist_type == :nominal
                        @scale_var_loop "4j2b" best_mass
                    end
                # HT HISTOGRAM
                elseif btag_count == 1 # no more than 1 btag
                    HT = @views sum(Jet_pt[jet_pt_mask])
                    push!(hists[hist_type]["4j1b"], HT, wgt)
                    if hist_type == :nominal
                        @scale_var_loop "4j1b" HT
                    end
                end
            end
        end
    end

    hists
end
