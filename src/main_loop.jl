function get_histo(process_tag::Symbol; wgt = 0.0, n_files_max_per_sample = MAX_N_FILES_PER_SAMPLE[])
    N = n_files_max_per_sample
    if iszero(wgt)
        wgt = LUMI * xsec_info[process_tag] / nevts_total(process_tag)
    end
    fs = @view TAG_PATH_DICT[process_tag][begin:N]
    mapreduce(mergewith(+), fs) do path #
        tree = LazyTree(path, "events")
        get_histo(tree, wgt)
    end
end

"""
    get_histo(tree, wgt; nbins=26, start=0, stop=375)
"""
function get_histo(tree, wgt; nbins=26, start=0, stop=375)
    hist_mass = Hist1D(Float64; bins = range(; start=start, stop=stop, length=nbins))
    hist_HT = Hist1D(Float64; bins = range(; start=start, stop=stop, length=nbins))
    
    for evt in tree
        # single lepton requirement
        (; Electron_pt, Muon_pt) = evt
        if count(>(25), Electron_pt) + count(>(25), Muon_pt) != 1
            continue
        end

        # at least 4 jets
        (; Jet_pt) = evt
        jet_pt_mask = Jet_pt .> 25
        count(jet_pt_mask) < 4 && continue

        jet_btag = @view evt.Jet_btagCSVV2[jet_pt_mask]

        # MASS HISTOGRAM
        if count(>=(0.5), jet_btag) >= 2 # at least 2 btag
            (; Jet_px, Jet_py, Jet_pz, Jet_mass) = evt

            # construct jet lorentz vector
            jet_p4 = @views fromPxPyPzM.(Jet_px[jet_pt_mask], Jet_py[jet_pt_mask], Jet_pz[jet_pt_mask], Jet_mass[jet_pt_mask])

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
            push!(hist_mass, best_mass, wgt)
        end

        # HT HISTOGRAM
        if count(>=(0.5), jet_btag) == 1 # no more than 1 btag
    
            HT = @views sum(Jet_pt[jet_pt_mask])
    
            # HT
            push!(hist_HT, HT, wgt)
        end
    end
    return Dict("4j2b" => hist_mass, "4j1b" => hist_HT) # the labels might change
end