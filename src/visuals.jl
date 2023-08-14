using CairoMakie

"""
    plot_mbjj_4j2b_stack(all_hists; evt_types::Union{Nothing, Vector{Symbol}}=nothing, syst_variation::Symbol=:mbjj_4j2b_nominal, color=Makie.wong_colors())


Creates a plot for the stack of histograms in the 4j2b region for the mbjj variable.


`all_hists` should be a dictionary of the form evt_type => hists_dict, where hists_dict is what get_histo(evt_type, ...) would return.
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

        # make legend
        labels = String.(evt_types)
        elements = [PolyElement(polycolor = color[i]) for i in 1:length(labels)]
        Legend(f[1,2], elements, labels)
        f
    end
end

"""
    plot_variation(hists, variable_region::Symbol, variations; color=Makie.wong_colors())

Usage example:
    ```
    plot_variation(all_hists[:ttbar], :mbjj_4j2b, [:nominal, :pt_scale_down, :pt_scale_up, :pt_res])
    ```
"""
function plot_variation(hists, variable_region::Symbol, variations; color=Makie.wong_colors())
	f = Figure()
	ax = Axis(f[1,1], title=String(variable_region)*" variations")
	for i=1:length(variations)
		v = variations[i]
		full_name = Symbol(variable_region, :_, v)
		stairs!(ax, hists[full_name], label=String(v), color=color[i])
	end

    # make legend
	labels = String.(variations)
    elements = [PolyElement(polycolor = color[i]) for i in 1:length(labels)]
    Legend(f[1,2], elements, labels)
	f
end