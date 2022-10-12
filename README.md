# LHC_AGC.jl [(Link to rendered notebook)](https://moelf.github.io/LHC_AGC.jl/dev/notebook/)

Analysis Grand Challenge @ LHC in pure Julia, this repo is set-up to render the Pluto notebook automatically, so feel free to contribute!

note: for RDataFrame implementation, see: https://github.com/andriiknu/RDF

# Instruction for running locally
1. download Julia somehow, preferably same version as:
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

using Pluto
Pluto.run()
```
5. Open up the `notebook.jl` and edit files in `/src` in your favourite editor, `Revise.jl` will auto hot reload the code for you.

The entry point to looper is this function:
https://github.com/Moelf/LHC_AGC.jl/blob/2106790a8e417ade474d37993e1d059d29528af0/src/LHC_AGC.jl#L38
