using CairoMakie

"""
    plot_mbjj_4j2b_stack(all_hists; evt_types::Union{Nothing, Vector{Symbol}}=nothing, syst_variation::Symbol=:mbjj_4j2b_nominal, color=Makie.wong_colors())

creates a plot for the stack of histograms in the 4j2b region for the mbjj variable.
if `evt_types` is `nothing` `keys(all_hists)` are used instead.
`color` should be a vector of colors of a meaningful size.
"""
function plot_mbjj_4j2b_stack(all_hists; evt_types::Union{Nothing, Vector{Symbol}}=nothing, syst_variation::Symbol=:mbjj_4j2b_nominal, color=Makie.wong_colors())
    with_theme(ATLASTHEME) do
        if evt_types === nothing
            evt_types = keys(all_hists)
        end
        f, ax, p = stackedhist([all_hists[evt_type][syst_variation] for evt_type in evt_types]; errors=true, color) # errors defaults to `true`
        ax.title = "4j2b"
        ax.xlabel = "mbjj"
        #a = Axis(f[1,1], title="4j2b", xlabel="mbjj")

        # make legend
        labels = String.(evt_types)
        elements = [PolyElement(polycolor = p.attributes.color[][i]) for i in 1:length(labels)]
        title = "legend"
        Legend(f[1,2], elements, labels, title)
        f
    end
end