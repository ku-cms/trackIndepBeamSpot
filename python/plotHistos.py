import ROOT
from inputFiles_cfi import get_list_of_files

# make sure ROOT.TFile.Open(fileURL) does not seg fault when $ is in sys.argv (e.g. $ passed in as argument)
ROOT.PyConfig.IgnoreCommandLineOptions = True
# make plots faster without displaying them
ROOT.gROOT.SetBatch(ROOT.kTRUE)
# Tell ROOT not to be in charge of memory, fix issue of histograms being deleted when ROOT file is closed:
ROOT.TH1.AddDirectory(False)

def plot(chain):
    c = ROOT.TCanvas("c", "c", 800, 800)
    chain.Draw("ClN")
    c.Update()
    c.SaveAs("plots/hist_ClN.png")

def run(input_files):
    chain = ROOT.TChain('pixelTree')
    for f in input_files:
        chain.Add(f)
    plot(chain)

def main():
    directory   = '/mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/user/caleb/PixelTrees/SingleMuon/crab_PixelTree_SingleMuon_2018C_RAW_Run319337_v1/210403_235502/0000'
    print "directory: {0}".format(directory)
    num_files   = 1
    input_files = get_list_of_files(directory)
    input_files = input_files[0:num_files]
    run(input_files)

if __name__ == "__main__":
    main()
    
