import ROOT
from tools import makeDir, getChain, get_file_list, get_eos_file_list

# Make sure ROOT.TFile.Open(fileURL) does not seg fault when $ is in sys.argv (e.g. $ passed in as argument)
ROOT.PyConfig.IgnoreCommandLineOptions = True
# Make plots faster without displaying them
ROOT.gROOT.SetBatch(ROOT.kTRUE)
# Tell ROOT not to be in charge of memory, fix issue of histograms being deleted when ROOT file is closed:
ROOT.TH1.AddDirectory(False)

def getCleanName(name):
    result  = name
    # map: keys replaced by values in name
    clean_map = {
        "=="    : "eq",
        "!="    : "neq",
        ">="    : "gte",
        "<="    : "lte",
        ">"     : "gt",
        "<"     : "lt",
        "&&"    : "and",
        "||"    : "or",
        " "     : "_"
    }
    # Order matters! Define order here:
    keys = ["==", "!=", ">=", "<=", ">", "<", "&&", "||", " "]
    for key in keys:
        value = clean_map[key]
        result = result.replace(key, value)
    return result

def plotChain(chain, name, var, cuts=""):
    plot_dir  = "plots/"
    makeDir(plot_dir)
    if cuts:
        cuts_name = getCleanName(cuts)
        plot_name = "{0}hist_{1}_{2}_{3}.pdf".format(plot_dir, name, var, cuts_name)
    else:
        plot_name = "{0}hist_{1}_{2}.pdf".format(plot_dir, name, var)
    
    c = ROOT.TCanvas("c", "c", 800, 800)

    chain.Draw(var, cuts)

    c.Update()
    c.SaveAs(plot_name)

# run over local files
def run(input_directory, num_files, name, isData):
    input_files = get_file_list(input_directory)
    chain       = getChain(input_files, num_files)
    makePlots(chain, name, isData)

# run over eos files
def runEOS(input_directory, num_files, name, isData):
    input_files = get_eos_file_list(input_directory)
    chain       = getChain(input_files, num_files)
    makePlots(chain, name, isData)

def makePlots(chain, name, isData):
    #plotChain(chain, name, "PvN", "ClN>=2000")
    #plotChain(chain, name, "PvN", "ClN<2000")
    #plotChain(chain, name, "ClN", "ClN>=2000")
    #plotChain(chain, name, "ClN", "ClN<2000")
    #plotChain(chain, name, "ClSize", "ClN>=2000")
    #plotChain(chain, name, "ClSize", "ClN>=2000 && ClSize>=2")
    
    plotChain(chain, name, "ClN")
    plotChain(chain, name, "ClN", "ClN>=2000")
    plotChain(chain, name, "ClN", "ClN>=2000 && ClSize>=2")
    #plotChain(chain, name, "ClN", "ClTkN<1")
    #plotChain(chain, name, "ClN", "ClTkN>=1")
    #plotChain(chain, name, "ClN", "ClTkN>=0")
    
    plotChain(chain, name, "ClSize")
    plotChain(chain, name, "ClSize", "ClSize>=2")
    plotChain(chain, name, "ClSize", "ClN>=2000 && ClSize>=2")
    
    #plotChain(chain, name, "ClSizeX")
    #plotChain(chain, name, "ClSizeY")
    #plotChain(chain, name, "ClCharge")
    #plotChain(chain, name, "ClChargeCorr")
    
    plotChain(chain, name, "PvN")
    plotChain(chain, name, "PvX")
    plotChain(chain, name, "PvY")
    plotChain(chain, name, "PvZ")
    plotChain(chain, name, "PvN", "ClN>=2000")
    plotChain(chain, name, "PvX", "ClN>=2000")
    plotChain(chain, name, "PvY", "ClN>=2000")
    plotChain(chain, name, "PvZ", "ClN>=2000")
    plotChain(chain, name, "PvN", "ClN>=2000 && ClSize>=2")
    plotChain(chain, name, "PvX", "ClN>=2000 && ClSize>=2")
    plotChain(chain, name, "PvY", "ClN>=2000 && ClSize>=2")
    plotChain(chain, name, "PvZ", "ClN>=2000 && ClSize>=2")
    #plotChain(chain, name, "PvX", "ClTkN<1")
    #plotChain(chain, name, "PvX", "ClTkN>=1")
    #plotChain(chain, name, "PvY", "ClTkN<1")
    #plotChain(chain, name, "PvY", "ClTkN>=1")
    #plotChain(chain, name, "PvZ", "ClTkN<1")
    #plotChain(chain, name, "PvZ", "ClTkN>=1")
    
    #plotChain(chain, name, "ClN", "run==297050")
    #plotChain(chain, name, "PvN", "run==297050")
    #plotChain(chain, name, "PvX", "run==297050")
    #plotChain(chain, name, "PvY", "run==297050")
    #plotChain(chain, name, "PvZ", "run==297050")
    
    #plotChain(chain, name, "ClN", "run==297050 && ClN>=2000")
    #plotChain(chain, name, "PvN", "run==297050 && ClN>=2000")
    #plotChain(chain, name, "PvX", "run==297050 && ClN>=2000")
    #plotChain(chain, name, "PvY", "run==297050 && ClN>=2000")
    #plotChain(chain, name, "PvZ", "run==297050 && ClN>=2000")
    
    if not isData:
        plotChain(chain, name, "BsX")
        plotChain(chain, name, "BsY")
        plotChain(chain, name, "BsZ")

def runData2017B():
    # # Min Bias 2017B
    # input_directory = '/home/t3-ku/calebsmith/store/PixelTrees/MinimumBias/crab_PixelTree_MinBias_2017B_RAW_v3/210715_152921/0000'
    # name            = "MinBias2017B"
    # isData          = True
    # num_files       = 30
    # run(input_directory, num_files, name, isData)
    
    # Zero Bias 2017B
    # input_directory = '/home/t3-ku/calebsmith/store/PixelTrees/ZeroBias/crab_PixelTree_ZeroBias_2017B_RAW_v2/210824_165140/0000'
    # name            = "ZeroBias_2017B_MoreEvents"
    # isData          = True
    # num_files       = 15
    # run(input_directory, num_files, name, isData)
    
    # Single Muon 2017B
    # input_directory = '/home/t3-ku/calebsmith/store/PixelTrees/SingleMuon/crab_PixelTree_SingleMuon_2017B_RAW_v1/210810_160009/0000'
    # name            = "SingleMuon_2017B_MoreEvents"
    # isData          = True
    # num_files       = 20
    # run(input_directory, num_files, name, isData)
    
    # Zero Bias 2017B v2: using golden json
    input_directory = '/home/t3-ku/calebsmith/store/PixelTrees/ZeroBias/crab_PixelTree_ZeroBias_2017B_RAW_v4/210909_195432/0000'
    name            = "ZeroBias_2017B_v2_MoreEvents"
    isData          = True
    #num_files       = 20
    num_files       = 2
    run(input_directory, num_files, name, isData)
    
    # Single Muon 2017B v2: using golden json
    input_directory = '/home/t3-ku/calebsmith/store/PixelTrees/SingleMuon/crab_PixelTree_SingleMuon_2017B_RAW_v2/210909_000356/0000'
    name            = "SingleMuon_2017B_v2_MoreEvents"
    isData          = True
    #num_files       = 20
    num_files       = 2
    run(input_directory, num_files, name, isData)

def runData2018C():
    # # Single Muon
    # input_directory = '/home/t3-ku/calebsmith/store/PixelTrees/SingleMuon/crab_PixelTree_SingleMuon_2018C_RAW_Run319337_v1/210403_235502/0000'
    # name            = "SingleMuon"
    # isData          = True
    # num_files       = 12
    # run(input_directory, num_files, name, isData)
    
    # # Zero Bias
    # input_directory = '/home/t3-ku/calebsmith/store/PixelTrees/ZeroBias/crab_PixelTree_ZeroBias_2018C_RAW_AllRuns_v1/210405_171418/0000'
    # name            = "ZeroBias"
    # isData          = True
    # num_files       = 2
    # run(input_directory, num_files, name, isData)
    
    # # Min Bias 2018C
    # input_directory = '/home/t3-ku/calebsmith/store/PixelTrees/MinimumBias2/crab_PixelTree_MinBias_2018C_RAW_v2/210707_165008/0000'
    # name            = "MinBias2018C"
    # isData          = True
    # num_files       = 2
    # run(input_directory, num_files, name, isData)
    
    # Single Muon 2018C
    input_directory = '/home/t3-ku/calebsmith/store/PixelTrees/SingleMuon/crab_PixelTree_SingleMuon_2018C_RAW_Run319337_v1/210403_235502/0000'
    name            = "SingleMuon_2018C_MoreEvents"
    isData          = True
    num_files       = 20
    run(input_directory, num_files, name, isData)
    
    # Zero Bias 2018C
    input_directory = '/home/t3-ku/calebsmith/store/PixelTrees/ZeroBias/crab_PixelTree_ZeroBias_2018C_RAW_AllRuns_v1/210405_171418/0000'
    name            = "ZeroBias_2018C_MoreEvents"
    isData          = True
    num_files       = 3
    run(input_directory, num_files, name, isData)

def runData2021():
    # directory on root://cmseos.fnal.gov
    # use eos functions for eos files
    input_directory = "/store/user/lpcsusylep/PixelTrees/ExpressPhysics/crab_PixelTree_Express_2021_Run346512_v1/211110_190014/0000"
    name            = "ExpressPhysics_2021_Run346512"
    isData          = True
    num_files       = -1
    runEOS(input_directory, num_files, name, isData)

def runData2022():
    # directory on root://cmseos.fnal.gov
    # use eos functions for eos files
    input_directory = "/store/user/lpcsusylep/PixelTrees/ZeroBias/crab_PixelTree_ZeroBias_2022F_RAW_v1_Run362154_v2/221206_141708/0000"
    name            = "ZeroBias_2022F_NumFiles1"
    isData          = True
    num_files       = 1
    runEOS(input_directory, num_files, name, isData)

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
    #runData2017B()
    #runData2018C()
    #runData2021()
    runData2022()
    #runMC()

if __name__ == "__main__":
    main()
 
