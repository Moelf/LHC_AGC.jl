### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ 524225b9-d165-47ab-9c49-dffbd93393aa
begin
	using Pkg
	Pkg.activate(@__DIR__)
	using CairoMakie, Revise, LHC_AGC, JSON3, FHist, ProfileCanvas
end

# ╔═╡ 97e9c1f4-d444-4db6-9037-28a68f2deccd
using BenchmarkTools

# ╔═╡ e484f0c5-73fc-40a3-b875-ee8a9e939e17
res = @time LHC_AGC.get_histo(:ttbar, 100);

# ╔═╡ 516e76f1-cc21-49b5-8a02-be47318e78ad
stairs(res, axis=(xlabel="mass (GeV)", ylabel="# events"))

# ╔═╡ Cell order:
# ╠═524225b9-d165-47ab-9c49-dffbd93393aa
# ╠═97e9c1f4-d444-4db6-9037-28a68f2deccd
# ╠═e484f0c5-73fc-40a3-b875-ee8a9e939e17
# ╠═516e76f1-cc21-49b5-8a02-be47318e78ad
