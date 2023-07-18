using LHC_AGC
using Test, UnROOT, FHist

function load_coffea_events(hist_type)
    coffea_evts = Int[]

    open("processed_coffea_"*hist_type, "r") do f
        for l in readlines(f)
            push!(coffea_evts, parse(Int, l))
        end
    end

    coffea_evts
end

"""
A function for test comparison
"""
function AGC_quicktest(filepath, _bincounts, wgt; eps=0.01)
    tt_tree = LazyTree(filepath, "Events")
    evts = Dict(k => Int[] for k in keys(_bincounts))
    res = LHC_AGC.get_histo(tt_tree, wgt; nbins=26, start=50, stop=550, evts=evts)

    for hist_type in keys(_bincounts)
        # test bincounts
        for region in keys(_bincounts[hist_type])
            @test (maximum(abs.(bincounts(res[hist_type][region]) - _bincounts[hist_type][region])) <= eps)
        end
        # test processed events
        @test isempty(setdiff(evts[hist_type], load_coffea_events(hist_type)))
    end

    nothing
end

@testset "LHC_AGC.jl" begin
    AGC_quicktest(
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
                "4j1b" => [ 0.,    0.,    0.,    7.,   60.,  185.,  461.,  830., 1102.,
                            1379., 1594., 1608., 1701., 1661., 1480., 1375., 1220., 1131.,
                            1028.,  871.,  736.,  660.,  551.,  477.,  397.],
                "4j2b" => [ 197.,  773., 1679., 2730., 3895., 5290., 6104., 4306., 3190.,
                            2560., 2015., 1627., 1421., 1176.,  971.,  860.,  739.,  630.,
                            586.,  517.,  462.,  435.,  371.,  342.,  301.]
            )
        ),

        1.0,
        #LHC_AGC.LUMI * LHC_AGC.xsec_info[:ttbar] / 225000

        eps=2
    )
end
