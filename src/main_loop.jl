function get_histo(process_tag::Symbol; wgt = 0.0, n_files_max_per_sample = MAX_N_FILES_PER_SAMPLE[], nbins=26, start=0, stop=375)
    N = n_files_max_per_sample
    if iszero(wgt)
        wgt = LUMI * xsec_info[process_tag] / nevts_total(process_tag)
    end
    fs = @view TAG_PATH_DICT[process_tag][begin:N]
    println(fs)
    hists = mapreduce(mergewith(mergewith(+)), fs) do path #
        println(path)
        tree = LazyTree(path, "Events")
        get_histo(tree, wgt; start=start, stop=stop, nbins=nbins)
    end
    h_nom = hists["nominal"]
    hists["luminocity_up"] = Dict(k => h_nom[k]*1.03 for k in keys(h_nom))
    hists["luminocity_down"] = Dict(k => h_nom[k]*0.97 for k in keys(h_nom))
    hists
end

"""
    get_histo(tree, wgt; nbins=26, start=0, stop=375)
"""
function get_histo(tree, wgt; nbins=26, start=0, stop=375)
    hists = Dict(
        "nominal" => Dict(
            "4j2b" => Hist1D(Float64; bins = range(; start=start, stop=stop, length=nbins)),
            "4j1b" => Hist1D(Float64; bins = range(; start=start, stop=stop, length=nbins))
        ),

        "pt_scale_up" => Dict(
            "4j2b" => Hist1D(Float64; bins = range(; start=start, stop=stop, length=nbins)),
            "4j1b" => Hist1D(Float64; bins = range(; start=start, stop=stop, length=nbins))
        ),

        "pt_scale_down" => Dict(
            "4j2b" => Hist1D(Float64; bins = range(; start=start, stop=stop, length=nbins)),
            "4j1b" => Hist1D(Float64; bins = range(; start=start, stop=stop, length=nbins))
        ),

        "pt_res" => Dict(
            "4j2b" => Hist1D(Float64; bins = range(; start=start, stop=stop, length=nbins)),
            "4j1b" => Hist1D(Float64; bins = range(; start=start, stop=stop, length=nbins))
        ),
    )
    pt_var = Dict(
        "nominal" => identity,
        "pt_scale_up" => (pt)->1.03pt,
        "pt_scale_down" => (pt)->0.97pt,
        "pt_res" => jet_pt_resolution
    )
    
    for evt in tree
        # single lepton requirement
        (; Electron_pt, Muon_pt) = evt
        (count(>(25), Electron_pt) + count(>(25), Muon_pt) != 1) && continue

        # get pt
        (; Jet_pt) = evt
        Jet_pt_nominal = Jet_pt

        for hist_type in keys(hists)
            # modify pt
            Jet_pt = pt_var[hist_type].(Jet_pt_nominal)

            # at least 4 jets
            jet_pt_mask = Jet_pt .> 25
            if count(jet_pt_mask) >= 4
                jet_btag = @view evt.Jet_btagCSVV2[jet_pt_mask]

                # MASS HISTOGRAM
                if count(>=(0.5), jet_btag) >= 2 # at least 2 btag
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
                        maximum(btags) < 0.5 && continue
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
                end

                # HT HISTOGRAM
                if count(>=(0.5), jet_btag) == 1 # no more than 1 btag
                    HT = @views sum(Jet_pt[jet_pt_mask])
                    push!(hists[hist_type]["4j1b"], HT, wgt)
                end
            end
        end
    end

    hists
end