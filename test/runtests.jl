using LHC_AGC
using Test, UnROOT, FHist

function load_coffea_events(hist_type)
    coffea_evts = Int[]

    open("processed_coffea_"*String(hist_type), "r") do f
        for l in readlines(f)
            push!(coffea_evts, parse(Int, l))
        end
    end

    coffea_evts
end

"""
A function for test comparison
"""
function AGC_quicktest(filepath, _bincounts, wgt, hist_types_to_match; eps=0.01) # this one runs the get_histo and compares the processed events by id as well as the bincounts
    tt_tree = LazyTree(filepath, "Events")
    evts = Dict(k => Int[] for k in hist_types_to_match)
    res = LHC_AGC.get_histo(tt_tree, wgt; evts=evts)

    for k in keys(_bincounts)
        # test bincounts
        @test (maximum(abs.(bincounts(res[k]) - _bincounts[k])) <= eps)
    end

    for hist_type in hist_types_to_match
        # test processed events
        @test isempty(setdiff(evts[hist_type], load_coffea_events(hist_type)))
    end

    res
end
function AGC_quicktest(res, _bincounts; eps=0.01) # this one doesn't run the get_histo and expects the result to be already given to it and only compares the bincounts
    for k in keys(_bincounts)
        # test bincounts
        @test (maximum(abs.(bincounts(res[k]) - _bincounts[k])) <= eps)
    end
    res
end

@testset "LHC_AGC.jl" begin
    # basic tests
    res = AGC_quicktest(
        "data/cmsopendata2015_ttbar_19981_PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext4-v1_80000_0007.root",

        Dict(
            :HT_4j1b_nominal => [ 0.,    0.,    0.,   10.,   71.,  256.,  568.,  978., 1253.,
                            1489., 1700., 1684., 1718., 1627., 1452., 1341., 1161., 1115.,
                            977.,  808.,  671.,  620.,  473.,  428.,  388.],
            :mbjj_4j2b_nominal => [ 217.,  801., 1766., 2884., 4174., 5670., 5925., 3930., 3086.,
                            2401., 1895., 1562., 1354., 1098.,  949.,  825.,  686.,  624.,
                            535.,  468.,  463.,  392.,  358.,  315.,  285.],
            :HT_4j1b_pt_scale_up => [ 0.,    0.,    0.,    7.,   60.,  185.,  461.,  831., 1102.,
                            1379., 1595., 1608., 1702., 1661., 1482., 1375., 1221., 1131.,
                            1028.,  871.,  736.,  660.,  551.,  477.,  396.],
            :mbjj_4j2b_pt_scale_up => [ 197.,  773., 1679., 2730., 3895., 5290., 6101., 4305., 3189.,
                            2560., 2015., 1625., 1420., 1176.,  971.,  859.,  739.,  630.,
                            586.,  517.,  462.,  435.,  371.,  342.,  300.]
        ),

        1.0, # using wgt=1.0 for clarity, should be (LHC_AGC.LUMI * LHC_AGC.xsec_info[:ttbar] / 225000) realistically

        (:nominal, :pt_scale_up),

        eps=2
    )

    # btag_var tests
    AGC_quicktest(
        res,

        Dict(
            :HT_4j1b_btag_var_0_up => [0.        ,    0.        ,    0.        ,   10.59873437,
                                    76.028     ,  277.56775003,  622.17878134, 1082.2217344 ,
                                    1398.25474995, 1675.03745329, 1928.10035937, 1922.50023447,
                                    1975.27395312, 1880.8956252 , 1691.31810932, 1572.80507809,
                                    1371.94771919, 1327.58604702, 1171.27085932,  978.92481253,
                                    817.11846835,  760.50706245,  585.77740631,  533.144     ,
                                    487.09825007],
            :mbjj_4j2b_btag_var_0_up => [243.69878125,  909.68610937, 2023.68795313, 3314.6561875 ,
                                    4845.58923438, 6684.84834375, 7084.89196875, 4681.40925   ,
                                    3681.8915    , 2873.84603125, 2283.42003125, 1891.38645313,
                                    1653.494875  , 1335.7666875 , 1158.01885937, 1014.36698438,
                                    845.81323437,  766.87603125,  660.23575   ,  583.7158125 ,
                                    576.48690625,  495.98346875,  444.51517188,  391.0558125 ,
                                    362.30765625],
            :HT_4j1b_btag_var_1_up => [0.        ,    0.        ,    0.        ,   10.46457812,
                                    74.84644533,  271.89008596,  607.97017187, 1053.09464075,
                                    1357.95214067, 1623.34060946, 1864.53940638, 1856.67709394,
                                    1906.95574994, 1812.89167182, 1627.94301566, 1511.32832811,
                                    1315.22090626, 1268.93871893, 1117.22909359,  932.51687491,
                                    777.49596882,  720.43999993,  553.24404685,  502.61187514,
                                    460.00675001],
            :mbjj_4j2b_btag_var_1_up => [237.11509375,  882.45391406, 1960.73227344, 3210.09683594,
                                    4680.17078125, 6430.480875  , 6788.63749219, 4487.83682812,
                                    3520.92209375, 2745.87507813, 2174.20491406, 1798.95577344,
                                    1567.99946875, 1265.7365625 , 1097.27790625,  957.22725   ,
                                    799.5483125 ,  724.09279687,  621.98246875,  547.52294531,
                                    540.5395    ,  461.30384375,  417.65064062,  366.9429375 ,
                                    337.039875  ],
            :HT_4j1b_btag_var_0_down => [0.        ,    0.        ,    0.        ,    9.40126563,
                                      65.972     ,  234.43224997,  513.82121866,  873.7782656 ,
                                      1107.74525005, 1302.96254671, 1471.89964063, 1445.49976553,
                                      1460.72604688, 1373.1043748 , 1212.68189068, 1109.19492191,
                                      950.05228081,  902.41395298,  782.72914068,  637.07518747,
                                      524.88153165,  479.49293755,  360.22259369,  322.856     ,
                                      288.90174993],
            :mbjj_4j2b_btag_var_0_down => [190.30121875,  692.31389063, 1508.31204688, 2453.3438125 ,
                                      3502.41076562, 4655.15165625, 4765.10803125, 3178.59075   ,
                                      2490.1085    , 1928.15396875, 1506.57996875, 1232.61354688,
                                      1054.505125  ,  860.2333125 ,  739.98114063,  635.63301562,
                                      526.18676562,  481.12396875,  409.76425   ,  352.2841875 ,
                                      349.51309375,  288.01653125,  271.48482812,  238.9441875 ,
                                      207.69234375],
            :HT_4j1b_btag_var_1_down => [  0.        ,    0.        ,    0.        ,    9.53542188,
                                        67.15355467,  240.10991404,  528.02982813,  902.90535925,
                                        1148.04785933, 1354.65939054, 1535.46059362, 1511.32290606,
                                        1529.04425006, 1441.10832818, 1276.05698434, 1170.67167189,
                                        1006.77909374,  961.06128107,  836.77090641,  683.48312509,
                                        564.50403118,  519.56000007,  392.75595315,  353.38812486,
                                        315.99324999],
            :mbjj_4j2b_btag_var_1_down => [  196.88490625,  719.54608594, 1571.26772656, 2557.90316406,
                                        3667.82921875, 4909.519125  , 5061.36250781, 3372.16317188,
                                        2651.07790625, 2056.12492188, 1615.79508594, 1325.04422656,
                                        1140.00053125,  930.2634375 ,  800.72209375,  692.77275   ,
                                        572.4516875 ,  523.90720313,  448.01753125,  388.47705469,
                                        385.4605    ,  322.69615625,  298.34935938,  263.0570625 ,
                                        232.960125  ]
        ),

        eps=0.0001
    )
end
