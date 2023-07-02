#=
function flat_variation(ones)
    # 2.5% weight variations
    return 1 + 0.025, 1 - 0.025
end
=#

function btag_weight_variation(jet_pt)
    pt_var = 0.025 * jet_pt / 50 #2.5% per 50 GeV
    1 + pt_var, 1 - pt_var
end

function jet_pt_resolution(pt::T; distr=Normal(1, 0.05)) where T
    T(rand(distr))*pt
end
