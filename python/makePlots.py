# makePlots.py

import ROOT

# Make sure ROOT.TFile.Open(fileURL) does not seg fault when $ is in sys.argv (e.g. $ passed in as argument)
ROOT.PyConfig.IgnoreCommandLineOptions = True
# Make plots faster without displaying them
ROOT.gROOT.SetBatch(ROOT.kTRUE)
# Tell ROOT not to be in charge of memory, fix issue of histograms being deleted when ROOT file is closed:
ROOT.TH1.AddDirectory(False)

def plot(input_file, plot_name):
    f = ROOT.TFile(input_file)
    h2d_occupancy = f.Get("h2d_occupancy")
    c1 = ROOT.TCanvas("c1", "c1", 800, 800)
    pad = c1.cd(1)
    h2d_occupancy.SetStats(0)
    h2d_occupancy.Draw("colz")
    c1.Update()
    c1.SaveAs(plot_name + ".pdf")
    c1.SaveAs(plot_name + ".png")

def makePlots():
    input_file = "output/TTBar_AllClusters_zsmear.root"
    plot_name  = "plots/TTBar_AllClusters_zsmear_occupancy"
    plot(input_file, plot_name)

def main():
    makePlots()

if __name__ == "__main__":
    main()
 
