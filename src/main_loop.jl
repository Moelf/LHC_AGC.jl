abstract type AbstractRegion end
struct Region4j1b <: AbstractRegion end
struct Region4j2b <: AbstractRegion end
# these constants should be passed to the function as region labels for now
# this way we utilise the multiple dispatch by having labels of different types
const reg4j1b = Region4j1b() 
const reg4j2b = Region4j2b()

function get_histo(tag::Symbol, region::AbstractRegion; wgt = 0.0, n_files_max_per_sample = MAX_N_FILES_PER_SAMPLE[])
    lumi = 3378 # /pb
    N = n_files_max_per_sample
    if iszero(wgt)
        wgt = lumi * xsec_info[tag] / nevts_total(tag)
    end
    fs = @view TAG_PATH_DICT[tag][begin:N]
    mapreduce(+, fs) do path
        tree = LazyTree(path, "events")
        get_histo(tree, wgt, region)
    end
end

"""
    get_histo(tree, wgt, region::Region4j2b; nbins=26, start=0, stop=375)
"""
function get_histo(tree, wgt, region::Region4j2b; nbins=26, start=0, stop=375)
    hist_mass = Hist1D(Float64; bins = range(; start=start, stop=stop, length=nbins))
    
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

"""
    get_histo(tree, wgt, region::Region4j1b; nbins=26, start=0, stop=375)
"""
function get_histo(tree, wgt, region::Region4j1b; nbins=26, start=0, stop=375)
    hist = Hist1D(Float64; bins = range(; start=start, stop=stop, length=nbins))
    
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

        # no more than 1 btag
        jet_btag = @view evt.jet_btag[jet_pt_mask]
        count(>=(0.5), jet_btag) >= 2 && continue

        HT = @views sum(jet_pt[jet_pt_mask])

        push!(hist, HT, wgt)
    end
    return hist
end