# makePlots.py

import ROOT

# Make sure ROOT.TFile.Open(fileURL) does not seg fault when $ is in sys.argv (e.g. $ passed in as argument)
ROOT.PyConfig.IgnoreCommandLineOptions = True
# Make plots faster without displaying them
ROOT.gROOT.SetBatch(ROOT.kTRUE)
# Tell ROOT not to be in charge of memory, fix issue of histograms being deleted when ROOT file is closed:
ROOT.TH1.AddDirectory(False)

def plot(input_file, plot_name, title, x_title, y_title):
    f = ROOT.TFile(input_file)
    h2d_occupancy = f.Get("h2d_occupancy")
    
    c1  = ROOT.TCanvas("c1", "c1", 800, 800)
    pad = c1.cd(1)
    pad.SetLeftMargin(0.1)
    pad.SetRightMargin(0.2)
    pad.SetTopMargin(0.1)
    pad.SetBottomMargin(0.1)
    
    x_axis = h2d_occupancy.GetXaxis()
    y_axis = h2d_occupancy.GetYaxis()
    
    x_axis.SetTitle(x_title)
    y_axis.SetTitle(y_title)
    h2d_occupancy.SetTitle(title)
    h2d_occupancy.SetStats(0)
    h2d_occupancy.Draw("colz")
    
    c1.Update()
    c1.SaveAs(plot_name + ".pdf")
    c1.SaveAs(plot_name + ".png")

def makePlots():
    input_file = "output/TTBar_AllClusters_zsmear.root"
    plot_name  = "plots/TTBar_AllClusters_zsmear_occupancy"
    title      = "TTBar Occupancy"
    x_title    = "Ring"
    y_title    = "Ladder"
    plot(input_file, plot_name, title, x_title, y_title)
    
    input_file = "output/SingleMuon_2017B_Legacy_MoreEvents_ClusterSize2_NumberClusters2000_AllClusters.root"
    plot_name  = "plots/SingleMuon_2017B_Legacy_MoreEvents_ClusterSize2_NumberClusters2000_AllClusters_occupancy"
    title      = "Single Muon 2017B Occupancy"
    x_title    = "Ring"
    y_title    = "Ladder"
    plot(input_file, plot_name, title, x_title, y_title)
    
    input_file = "output/ZeroBias_2017B_Legacy_MoreEvents_ClusterSize2_NumberClusters2000_AllClusters.root"
    plot_name  = "plots/ZeroBias_2017B_Legacy_MoreEvents_ClusterSize2_NumberClusters2000_AllClusters_occupancy"
    title      = "Zero Bias 2017B Occupancy"
    x_title    = "Ring"
    y_title    = "Ladder"
    plot(input_file, plot_name, title, x_title, y_title)

def main():
    makePlots()

if __name__ == "__main__":
    main()
 
