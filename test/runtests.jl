using LHC_AGC
using Test, UnROOT, FHist

"""
A function for test comparison
"""
function AGC_quicktest(filepath, _bincounts; eps=0.001)
    tt_tree = LazyTree(filepath, "Events")
    res = LHC_AGC.get_histo(tt_tree, 1.0, nbins=26, start=50, stop=550)

    for hist_type in keys(_bincounts)
        for region in keys(_bincounts[hist_type])
            if sum(abs.(bincounts(res[hist_type][region]) - _bincounts[hist_type][region])) > eps
                return false
            end
        end
    end

    true
end

@testset "LHC_AGC.jl" begin
    #=@test AGC_quicktest(
        "data/cmsopendata2015_ttbar_19981_PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext4-v1_80000_0007.root",

        Dict(
            "nominal" => Dict( 
                "4j1b" => [0.0, 0.0, 0.0, 10.0, 71.0, 256.0, 568.0, 977.0, 1253.0, 1488.0, 1700.0, 1684.0, 1717.0, 1626.0, 1451.0, 1340.0, 1161.0, 1115.0, 977.0, 808.0, 671.0, 620.0, 473.0, 429.0, 388.0],
                "4j2b" => [217.0, 801.0, 1766.0, 2884.0, 4175.0, 5671.0, 5927.0, 3930.0, 3088.0, 2400.0, 1897.0, 1562.0, 1355.0, 1098.0, 949.0, 826.0, 686.0, 624.0, 535.0, 468.0, 463.0, 392.0, 358.0, 315.0, 285.0]
            )
        )
    )=#
end
