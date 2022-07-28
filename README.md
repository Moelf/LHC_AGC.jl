# LHC_AGC.jl
Analysis Grand Challenge @ LHC in pure Julia

For ROOT RDataFrame implementation, see: https://github.com/andriiknu/RDF

# Running instruction
**Note: currently only works if you're on af.uchicago.edu because that's where files are.**

1. download Julia somehow
2. Open Julia REPL, hit `]` once to enter pkg mode, run
```julia
(v1.7) pkg> dev https://github.com/Moelf/LHC_AGC.jl/
```
3. `cd ~/.julia/dev/LHC_AGC`
4. `julia --project=.`
5. run:
```julia
using Pluto
Pluto.run()
```
6. Open up the `notebook.jl` and edit fiels in `/src` in your favourite editor.

The entry point to looper is this function:
https://github.com/Moelf/LHC_AGC.jl/blob/2106790a8e417ade474d37993e1d059d29528af0/src/LHC_AGC.jl#L38
