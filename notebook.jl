### A Pluto.jl notebook ###
# v0.19.12

using Markdown
using InteractiveUtils

# ╔═╡ 524225b9-d165-47ab-9c49-dffbd93393aa
# hideall
begin
	using Pkg
	Pkg.activate(@__DIR__)
    Pkg.instantiate()
	using CairoMakie, Revise, LHC_AGC, JSON3, FHist, UnROOT
end

# ╔═╡ 97e9c1f4-d444-4db6-9037-28a68f2deccd
using BenchmarkTools, PlutoUI # hide

# ╔═╡ 5f403020-c4ba-4814-926d-389e2347e9d2
# by default download 1 file per sample
LHC_AGC.download_data()

# ╔═╡ c82d8346-fa92-4159-b4ac-85ecfe035297
md"""
## Introduction: exploring $t\bar{t}$
"""

# ╔═╡ c0175af3-2903-4def-8599-da934fc28803
const tt_tree = LazyTree(joinpath(LHC_AGC.BASE_PATH[], "00DF0A73-17C2-E511-B086-E41D2D08DE30.root"), "events"); # hide

# ╔═╡ 87017bd9-8622-43c5-aa44-f24a6c1b39b1
njet_hist = Hist1D((length(x) for x in tt_tree.jet_pt), 0:16; overflow=true)

# ╔═╡ 72be7e4c-d85a-4795-ba36-4e83771ba160
# chop off the overflow bin
integral(restrict(njet_hist, 0, 15))

# ╔═╡ 92cd3a46-10bf-4b52-ad1a-ae07459830c8
with_terminal() do
	# this is a lazy generator, evident since it has 0 allocation
	global _c = @time (count(>(25), x) for x in tt_tree.jet_pt)
end

# ╔═╡ 33b8830a-8365-4a37-ad12-180113101700
with_terminal() do
	# and the loop over this generator is fully inferred, which means it's fast
	@code_warntype Hist1D(_c, 2:10)
end

# ╔═╡ dbfd0d3c-2b83-48ba-b2df-ba93a4a7d34d
Hist1D(_c, 2:10)

# ╔═╡ f80e56a0-2ffb-44c8-a11d-6922c96d3188
md"""
## $t\bar{t}$ kinematic reconstruction of the top mass.

You can find the full function in the `main_looper.jl` file, and look below this cell to look at the result of calling the function. Here we present a snippt view of the main three parts of the main looper:

### Basic cuts
```julia
	# single lepton requirement
	# Avoid allocation:
	# - bad example: `length(electron_pt .> 25)`
	if count(>(25), electron_pt) + count(>(25), muon_pt) != 1
		continue
	end

	# at least 4 jets
	(; jet_pt) = evt
	# can't avoid allocation since we need this mask later
	jet_pt_mask = jet_pt .> 25
	count(jet_pt_mask) < 4 && continue

	# at least 2 btag
	# use view to avoid allocation
	jet_btag = @view evt.jet_btag[jet_pt_mask]
	count(>=(0.5), jet_btag) < 2 && continue
	(; jet_px, jet_py, jet_pz, jet_mass) = evt
```

### Tri-jet combinatorics
```julia
	# construct jet lorentz vector
	jet_p4 = @views fromPxPyPzM.(jet_px[jet_pt_mask], jet_py[jet_pt_mask], jet_pz[jet_pt_mask], jet_mass[jet_pt_mask])

	# tri jet combinatorics
	max_pt = -Inf
	local best_mass
	
	# 1. all tri-jet combinations
	# again it's lazy, avoid allocating a bunch of 3-tuples
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

```
"""

# ╔═╡ e484f0c5-73fc-40a3-b875-ee8a9e939e17
# let's not weight this one
res = @time LHC_AGC.get_histo(:ttbar; wgt=1.0, n_files_max_per_sample=1);

# ╔═╡ 516e76f1-cc21-49b5-8a02-be47318e78ad
begin
	stairs(res, 
		axis=(xlabel="mass (GeV)", ylabel="# events", 
		limits=(0,375,0,1000),
		yticks=0:200:1000,
		)
	)
	errorbars!(res)
	vlines!(175, color=:grey, linestyle=:dash)
    text!(180, 20, text = L"m_t = 175 \mathrm{GeV}")
	current_figure()
end

# ╔═╡ Cell order:
# ╠═524225b9-d165-47ab-9c49-dffbd93393aa
# ╠═97e9c1f4-d444-4db6-9037-28a68f2deccd
# ╠═5f403020-c4ba-4814-926d-389e2347e9d2
# ╟─c82d8346-fa92-4159-b4ac-85ecfe035297
# ╠═c0175af3-2903-4def-8599-da934fc28803
# ╠═87017bd9-8622-43c5-aa44-f24a6c1b39b1
# ╠═72be7e4c-d85a-4795-ba36-4e83771ba160
# ╠═92cd3a46-10bf-4b52-ad1a-ae07459830c8
# ╠═33b8830a-8365-4a37-ad12-180113101700
# ╠═dbfd0d3c-2b83-48ba-b2df-ba93a4a7d34d
# ╟─f80e56a0-2ffb-44c8-a11d-6922c96d3188
# ╠═e484f0c5-73fc-40a3-b875-ee8a9e939e17
# ╠═516e76f1-cc21-49b5-8a02-be47318e78ad
