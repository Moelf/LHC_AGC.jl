function btag_weight_variation(jet_pts)
    pt_var = 0.075 * jet_pts / 50
    return 1 + pt_var, 1 - pt_var
end

function jet_pt_resolution(pt::T; distr=Normal(1, 0.05)) where T
    T(rand(distr))*pt
end

const SCALE_VARS = (
    scale_var = evt -> (1.025f0*evt.wgt, 0.975f0*evt.wgt),
    # the following code repetition is artificial
    # they essentially act as different sys, AGC just picked same base name
    btag_var_0 = evt -> btag_weight_variation(evt.Jet_pt[1]),
    btag_var_1 = evt -> btag_weight_variation(evt.Jet_pt[2]),
    btag_var_2 = evt -> btag_weight_variation(evt.Jet_pt[3]),
    btag_var_3 = evt -> btag_weight_variation(evt.Jet_pt[4])
)

const SHAPE_VARS = (
    nominal = identity,
    pt_scale_up = (pt)->1.03f0 * pt,
    pt_scale_down = (pt)->0.97f0 * pt,
    pt_res = jet_pt_resolution
)

macro scale_var_loop(region, phys_var)
    exs = Expr[]
    for scale_name in keys(SCALE_VARS)
        sym = QuoteNode(scale_name)
        ex = quote
            up, down = SCALE_VARS[$sym](scale_info)
            push!(hists[Symbol($sym, :_up)][$region], $phys_var, up*wgt)
            push!(hists[Symbol($sym, :_down)][$region], $phys_var, down*wgt)
        end

        push!(exs, ex)
    end

    esc(Expr(:block, exs...))
end

