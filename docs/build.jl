tutorials_dir = dirname(@__DIR__)

"""Run all Pluto notebooks (".jl" files) in `tutorials_dir` and write outputs to HTML files."""
function build()
    println("Building notebooks")
    # Evaluate notebooks in the same process to avoid having to recompile from scratch each time.
    # This is similar to how Documenter and Franklin evaluate code.
    # Note that things like method overrides and other global changes may leak between notebooks!
    use_distributed = false
    output_format = documenter_output
    bopts = BuildOptions(tutorials_dir; use_distributed, output_format)
    build_notebooks(bopts)
    return nothing
end

"Return Markdown file links which can be passed to Documenter.jl."
function markdown_files()
    md_files = map(notebooks) do notebook
        file = lowercase(replace(notebook, " " => '_'))
        fmd = "$file.md"
        mv(joinpath(tutorials_dir, fmd), joinpath(@__DIR__, "src", fmd); force=true)
        return fmd
    end
    return md_files
end
