import ROOT
from inputFiles_cfi import get_list_of_files

# make sure ROOT.TFile.Open(fileURL) does not seg fault when $ is in sys.argv (e.g. $ passed in as argument)
ROOT.PyConfig.IgnoreCommandLineOptions = True
# make plots faster without displaying them
ROOT.gROOT.SetBatch(ROOT.kTRUE)
# Tell ROOT not to be in charge of memory, fix issue of histograms being deleted when ROOT file is closed:
ROOT.TH1.AddDirectory(False)

def plotChain(chain, name, var, cuts=""):
    plot_dir  = "plots/"
    if cuts:
        plot_name = "{0}hist_{1}_{2}_{3}.png".format(plot_dir, name, var, cuts)
    else:
        plot_name = "{0}hist_{1}_{2}.png".format(plot_dir, name, var)
    
    c = ROOT.TCanvas("c", "c", 800, 800)

    chain.Draw(var, cuts)

    c.Update()
    c.SaveAs(plot_name)

def run(input_directory, num_files, name, isData):
    input_files = get_list_of_files(input_directory)
    input_files = input_files[0:num_files]
    
    chain = ROOT.TChain('pixelTree')
    for f in input_files:
        chain.Add(f)
    
    
    plotChain(chain, name, "PvN", "ClN>=2000")
    plotChain(chain, name, "PvN", "ClN<2000")
    plotChain(chain, name, "ClN", "ClN>=2000")
    plotChain(chain, name, "ClN", "ClN<2000")
    
    #plotChain(chain, name, "ClN")
    #plotChain(chain, name, "ClN", "ClTkN<1")
    #plotChain(chain, name, "ClN", "ClTkN>=1")
    #plotChain(chain, name, "ClN", "ClTkN>=0")
    
    #plotChain(chain, name, "ClSize")
    #plotChain(chain, name, "ClSizeX")
    #plotChain(chain, name, "ClSizeY")
    #plotChain(chain, name, "ClCharge")
    #plotChain(chain, name, "ClChargeCorr")
    
    #plotChain(chain, name, "PvN")
    #plotChain(chain, name, "PvX")
    #plotChain(chain, name, "PvY")
    #plotChain(chain, name, "PvZ")
    #plotChain(chain, name, "PvX", "ClTkN<1")
    #plotChain(chain, name, "PvX", "ClTkN>=1")
    #plotChain(chain, name, "PvY", "ClTkN<1")
    #plotChain(chain, name, "PvY", "ClTkN>=1")
    #plotChain(chain, name, "PvZ", "ClTkN<1")
    #plotChain(chain, name, "PvZ", "ClTkN>=1")
    
    if not isData:
        plotChain(chain, name, "BsX")
        plotChain(chain, name, "BsY")
        plotChain(chain, name, "BsZ")

def runData2017B():
    # # Min Bias 2017B
    # input_directory = '/mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/user/caleb/PixelTrees/MinimumBias/crab_PixelTree_MinBias_2017B_RAW_v3/210715_152921/0000'
    # name            = "MinBias2017B"
    # isData          = True
    # num_files       = 30
    # run(input_directory, num_files, name, isData)
    
    # Single Muon 2017B
    input_directory = '/mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/user/caleb/PixelTrees/SingleMuon/crab_PixelTree_SingleMuon_2017B_RAW_v1/210810_160009/0000'
    name            = "SingleMuon_2017B_MoreEvents"
    isData          = True
    num_files       = 20
    run(input_directory, num_files, name, isData)

    # Zero Bias 2017B
    input_directory = '/mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/user/caleb/PixelTrees/ZeroBias/crab_PixelTree_ZeroBias_2017B_RAW_v2/210824_165140/0000'
    name            = "ZeroBias_2017B_MoreEvents"
    isData          = True
    num_files       = 15
    run(input_directory, num_files, name, isData)

def runData2018C():
    # # Single Muon
    # input_directory = '/mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/user/caleb/PixelTrees/SingleMuon/crab_PixelTree_SingleMuon_2018C_RAW_Run319337_v1/210403_235502/0000'
    # name            = "SingleMuon"
    # isData          = True
    # num_files       = 12
    # run(input_directory, num_files, name, isData)
    
    # # Zero Bias
    # input_directory = '/mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/user/caleb/PixelTrees/ZeroBias/crab_PixelTree_ZeroBias_2018C_RAW_AllRuns_v1/210405_171418/0000'
    # name            = "ZeroBias"
    # isData          = True
    # num_files       = 2
    # run(input_directory, num_files, name, isData)
    
    # # Min Bias 2018C
    # input_directory = '/mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/user/caleb/PixelTrees/MinimumBias2/crab_PixelTree_MinBias_2018C_RAW_v2/210707_165008/0000'
    # name            = "MinBias2018C"
    # isData          = True
    # num_files       = 2
    # run(input_directory, num_files, name, isData)
    
    # Single Muon 2018C
    input_directory = '/mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/user/caleb/PixelTrees/SingleMuon/crab_PixelTree_SingleMuon_2018C_RAW_Run319337_v1/210403_235502/0000'
    name            = "SingleMuon_2018C_MoreEvents"
    isData          = True
    num_files       = 20
    run(input_directory, num_files, name, isData)
    
    # Zero Bias 2018C
    input_directory = '/mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/user/caleb/PixelTrees/ZeroBias/crab_PixelTree_ZeroBias_2018C_RAW_AllRuns_v1/210405_171418/0000'
    name            = "ZeroBias_2018C_MoreEvents"
    isData          = True
    num_files       = 3
    run(input_directory, num_files, name, isData)

def runMC():
    # TTbar (Z smeared)
    input_directory = '/mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/user/eschmitz/PixelTrees/RelValTTbar_13TeV/crab_RelValTTbar_13TeVdesign_0p1_neg0p08_GEN/190819_222045/0000'
    name            = "TTbar_Zsmeared"
    isData          = False
    num_files       = 100
    run(input_directory, num_files, name, isData)
    
    # TTbar (not Z smeared)
    input_directory = '/mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/user/eschmitz/PixelTrees/nosmear_newZ/RelValTTbar_13TeV/crab_RelValTTbar_13TeV0p1_neg0p08_GEN_SIM/200519_042750/0000/'
    name            = "TTbar_NotZsmeared"
    isData          = False
    num_files       = 200
    run(input_directory, num_files, name, isData)
    
    # TTbar with pileup (x = 0.2)
    input_directory = '/mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/user/eschmitz/PixelTrees/RelValTTbar_13TeV/crab_RelValTTbar_13TeVRAW_eschmitzcrab_design_0p2/190930_182712/0000'
    name            = "TTbar_Pileup_0p2"
    isData          = False
    num_files       = 100
    run(input_directory, num_files, name, isData)

    # TTbar without pileup (x = 0.2)
    input_directory = '/mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/user/eschmitz/PixelTrees/RelValTTbar_13TeV/crab_RelValTTbar_13TeVdesign_0p2_GEN_SIM/190819_222136/0000'
    name            = "TTbar_NoPileup_0p2"
    isData          = False
    num_files       = 100
    run(input_directory, num_files, name, isData)

def main():
    runData2017B()
    #runData2018C()
    #runMC()

if __name__ == "__main__":
    main()
    
