# LHC_AGC.jl

Analysis Grand Challenge @ LHC in pure Julia, this repo is set-up to render the Pluto notebook automatically, so feel free to contribute!

note: for RDataFrame implementation, see: https://github.com/andriiknu/RDF

# Instruction for running locally
1. download Julia, for example, through [juliaup](https://github.com/JuliaLang/juliaup), preferably same version as:
https://github.com/Moelf/LHC_AGC.jl/blob/master/Manifest.toml#L1-L5

but any version within the same LTS should be fine.

2. in your shell, run:
```bash
mkdir -p ~/.julia/environments/
git clone https://github.com/Moelf/LHC_AGC.jl/ ~/.julia/environments/LHC_AGC
```
3. `julia --project=@LHC_AGC`
4. run:
```julia
using Pkg
Pkg.instantiate() #only needed first time
```
