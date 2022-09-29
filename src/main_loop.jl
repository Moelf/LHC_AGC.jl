function get_histo(tag::Symbol; wgt = 0.0, n_files_max_per_sample = MAX_N_FILES_PER_SAMPLE[])
    lumi = 3378 # /pb
    N = n_files_max_per_sample
    if iszero(wgt)
        wgt = lumi * xsec_info[tag] / nevts_total(tag)
    end
    fs = @view TAG_PATH_DICT[tag][begin:N]
    mapreduce(+, fs) do path
        tree = LazyTree(path, "events")
        get_histo(tree, wgt)
    end
end

function get_histo(tree, wgt; nbins=26, bin_low=50, bin_high=550)
    hist_mass = Hist1D(Float64; bins = range(; start=0, stop=375, length=nbins))
    
    for evt in tree

        # single lepton requirement
        (; electron_pt, muon_pt) = evt
        if count(>(25), electron_pt) + count(>(25), muon_pt) != 1
            continue
        end

        # at least 4 jets
        (; jet_pt) = evt
        jet_pt_mask = jet_pt .> 25
        count(jet_pt_mask) < 4 && continue

        # at least 2 btag
        jet_btag = @view evt.jet_btag[jet_pt_mask]
        count(>=(0.5), jet_btag) < 2 && continue
        (; jet_px, jet_py, jet_pz, jet_mass) = evt

        # construct jet lorentz vector
        jet_p4 = @views fromPxPyPzM.(jet_px[jet_pt_mask], jet_py[jet_pt_mask], jet_pz[jet_pt_mask], jet_mass[jet_pt_mask])

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
    return hist_mass
end
