module LHC_AGC

using UnROOT, FHist, LorentzVectorHEP, JSON3
using LorentzVectorHEP: fromPxPyPzM
using Combinatorics: combinations

const xsec_info = Dict(
    :ttbar => 396.87 + 332.97, # nonallhad + allhad, keep same x-sec for all
    :single_top_s_chan => 2.0268 + 1.2676,
    :single_top_t_chan => (36.993 + 22.175)/0.252,  # scale from lepton filter to inclusive
    :single_top_tW => 37.936 + 37.906,
    :wjets => 61457 * 0.252,  # e/mu+nu final states
    :data => 1.0
)
const NJSON = JSON3.read(read("./ntuples.json"));
const TAGS = keys(NJSON)

const BASE_PATH = Ref("/data/jiling/Analysis_Grand_Challenge/")

function xrd_to_local(url)
    joinpath(BASE_PATH[], last(split(url, '/')))
end


const TAG_PATH_DICT = Dict(k=>LHC_AGC.xrd_to_local.(getindex.(LHC_AGC.NJSON[k][:nominal][:files], :path))
       for k in LHC_AGC.TAGS
      )

function get_histo(tag::Symbol, N)
    x_sec = xsec_info[tag]
    fs = @view TAG_PATH_DICT[tag][begin:N]
    mapreduce(+, fs) do path
        tree = LazyTree(path, "events")
        get_histo(tree, x_sec)
    end
end

function get_histo(tree, x_sec; nbins=25, bin_low=50, bin_high=550)
    nevts_total = length(tree)
    lumi = 3378 # /pb
    wgt = x_sec * lumi / nevts_total
    
    hist_mass = Hist1D(Float64; bins = 50:20:550)
    
    for evt in tree

        # single lepton req.
        (; electron_pt, muon_pt) = evt
        electron_pt_mask, muon_pt_mask = electron_pt .> 25, muon_pt .> 25
        if count(electron_pt_mask) + count(muon_pt_mask) != 1
            continue
        end

        # 4 jets req.
        (; jet_pt) = evt
        jet_pt_mask = jet_pt .> 25
        count(jet_pt_mask) < 4 && continue

        # btag req.
        jet_btag = @view evt.jet_btag[jet_pt_mask]
        count(>=(0.5), jet_btag) < 2 && continue
        (; jet_px, jet_py, jet_pz, jet_mass) = evt

        # construct jet lorentz vector
        jet_p4 = @views fromPxPyPzM.(jet_px[jet_pt_mask], jet_py[jet_pt_mask], jet_pz[jet_pt_mask], jet_mass[jet_pt_mask])

        length(jet_btag) == length(jet_p4) || error("impossible reached")

        # tri jet combinatorics
        max_pt = -Inf
        local best_mass
        
        # 1. all tri-jet combinations
        for (p4s, btags) in zip(combinations(jet_p4, 3), combinations(jet_btag, 3))
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

end
