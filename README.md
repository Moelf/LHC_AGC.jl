# LHC_AGC.jl
Analysis Grand Challenge @ LHC in pure Julia

For ROOT RDataFrame implementation, see: https://github.com/andriiknu/RDF

# Running instruction
**Note: currently only works if you're on af.uchicago.edu because that's where files are.**

1. download Julia somehow
2. in your shell, run:
```bash
mkdir -p ~/.julia/environments/
git clone https://github.com/Moelf/LHC_AGC.jl/ ~/.julia/environments/LHC_AGC
```
3. `julia --project=@LHC_AGC`
4. run:
```julia
using Pkg
Pkg.instantiate()

using Pluto
Pluto.run()
```
5. Open up the `notebook.jl` and edit files in `/src` in your favourite editor, `Revise.jl` will auto hot reload the code for you.

The entry point to looper is this function:
https://github.com/Moelf/LHC_AGC.jl/blob/2106790a8e417ade474d37993e1d059d29528af0/src/LHC_AGC.jl#L38
