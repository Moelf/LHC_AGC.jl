using LHC_AGC
using Test, UnROOT, FHist

"""
A function for test comparison
"""
function AGC_quicktest(filepath, _bincounts, wgt; eps=0.0001)
    tt_tree = LazyTree(filepath, "Events")
    res = LHC_AGC.get_histo(tt_tree, wgt, nbins=26, start=50, stop=550)

    for hist_type in keys(_bincounts)
        for region in keys(_bincounts[hist_type])
            @test (maximum(abs.(bincounts(res[hist_type][region]) - _bincounts[hist_type][region])) <= eps)
        end
    end

    true
end

@testset "LHC_AGC.jl" begin
    #=AGC_quicktest(
        "data/cmsopendata2015_ttbar_19981_PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext4-v1_80000_0007.root",

        Dict(
            "nominal" => Dict( 
                "4j1b" => [ 0.,    0.,    0.,   10.,   71.,  256.,  568.,  977., 1253.,
                            1488., 1700., 1684., 1717., 1626., 1451., 1340., 1161., 1115.,
                            977.,  808.,  671.,  620.,  473.,  429.,  388. ],
                "4j2b" => [ 217.,  801., 1766., 2884., 4174., 5671., 5927., 3931., 3086.,
                            2402., 1896., 1563., 1355., 1098.,  949.,  826.,  686.,  624.,
                            535.,  468.,  463.,  392.,  358.,  316.,  285. ]
            ),
            "pt_scale_up" => Dict( 
                "4j1b" => [ 0.,    0.,    0.,   10.,   71.,  256.,  568.,  977., 1253.,
                            1488., 1700., 1684., 1717., 1626., 1451., 1340., 1161., 1115.,
                            977.,  808.,  671.,  620.,  473.,  429.,  388. ],
                "4j2b" => [ 0.,    0.,    0.,    7.,   60.,  185.,  461.,  830., 1102.,
                            1379., 1594., 1608., 1701., 1661., 1480., 1375., 1220., 1131.,
                            1028.,  871.,  736.,  660.,  551.,  477.,  397.]
            )
        )

        1.0
        #LHC_AGC.LUMI * LHC_AGC.xsec_info[:ttbar] / 225000
    ),=#
end
