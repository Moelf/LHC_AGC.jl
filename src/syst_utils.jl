#=
function flat_variation(ones)
    # 2.5% weight variations
    return 1 + 0.025, 1 - 0.025
end
=#

function btag_weight_variation(jet_pts)
    pt_var = 0.075 * jet_pts / 50
    (1 .+ pt_var), (1 .- pt_var)
end

function jet_pt_resolution(pt::T; distr=Normal(1, 0.05)) where T
    T(rand(distr))*pt
end
