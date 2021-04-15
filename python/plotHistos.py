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
    plotChain(chain, name, "ClN")
    plotChain(chain, name, "ClN", "ClTkN<1")
    plotChain(chain, name, "ClN", "ClTkN>=1")
    plotChain(chain, name, "PvN")
    plotChain(chain, name, "PvX")
    plotChain(chain, name, "PvY")
    plotChain(chain, name, "PvZ")
    if not isData:
        plotChain(chain, name, "BsX")
        plotChain(chain, name, "BsY")
        plotChain(chain, name, "BsZ")

def runData():
    input_directory = '/mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/user/caleb/PixelTrees/SingleMuon/crab_PixelTree_SingleMuon_2018C_RAW_Run319337_v1/210403_235502/0000'
    name            = "SingleMuon"
    isData          = True
    num_files       = 12
    #num_files       = 1
    run(input_directory, num_files, name, isData)
    
    input_directory = '/mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/user/caleb/PixelTrees/ZeroBias/crab_PixelTree_ZeroBias_2018C_RAW_AllRuns_v1/210405_171418/0000'
    name            = "ZeroBias"
    isData          = True
    num_files       = 2
    #num_files       = 1
    run(input_directory, num_files, name, isData)

def runMC():
    input_directory = '/mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/user/eschmitz/PixelTrees/RelValTTbar_13TeV/crab_RelValTTbar_13TeVdesign_0p1_neg0p08_GEN/190819_222045/0000'
    name            = "TTbar"
    isData          = False
    num_files       = 100
    run(input_directory, num_files, name, isData)

def main():
    #runData()
    runMC()
    

if __name__ == "__main__":
    main()
    
