# makePlots.py

import ROOT

# Make sure ROOT.TFile.Open(fileURL) does not seg fault when $ is in sys.argv (e.g. $ passed in as argument)
ROOT.PyConfig.IgnoreCommandLineOptions = True
# Make plots faster without displaying them
ROOT.gROOT.SetBatch(ROOT.kTRUE)
# Tell ROOT not to be in charge of memory, fix issue of histograms being deleted when ROOT file is closed:
ROOT.TH1.AddDirectory(False)

def plot2D(input_file, plot_name, hist_name, title, x_title, y_title):
    f = ROOT.TFile(input_file)
    h2d = f.Get(hist_name)
    
    c1  = ROOT.TCanvas("c1", "c1", 800, 800)
    pad = c1.cd(1)
    pad.SetLeftMargin(0.1)
    pad.SetRightMargin(0.2)
    pad.SetTopMargin(0.1)
    pad.SetBottomMargin(0.1)
    
    x_axis = h2d.GetXaxis()
    y_axis = h2d.GetYaxis()
    
    x_axis.SetTitle(x_title)
    y_axis.SetTitle(y_title)
    h2d.SetTitle(title)
    h2d.SetStats(0)
    h2d.Draw("colz")
    
    c1.Update()
    c1.SaveAs(plot_name + ".pdf")

def makePlots():
    input_file = "output/TTBar_AllClusters_zsmear.root"
    plot_name  = "plots/TTBar_AllClusters_zsmear_occupancy"
    hist_name  = "h2d_occupancy"
    title      = "TTBar Occ."
    x_title    = "Ring"
    y_title    = "Ladder"
    plot2D(input_file, plot_name, hist_name, title, x_title, y_title)
    
    input_file = "output/TTBar_AllClusters_zsmear.root"
    plot_name  = "plots/TTBar_AllClusters_zsmear_occupancy_cut"
    hist_name  = "h2d_occupancy_cut"
    title      = "TTBar Occ. (cut)"
    x_title    = "Ring"
    y_title    = "Ladder"
    plot2D(input_file, plot_name, hist_name, title, x_title, y_title)
    
    input_file = "output/TTBar_AllClusters_zsmear.root"
    plot_name  = "plots/TTBar_AllClusters_zsmear_occupancy_mask"
    hist_name  = "h2d_occupancy_mask"
    title      = "TTBar Occ. (mask)"
    x_title    = "Ring"
    y_title    = "Ladder"
    plot2D(input_file, plot_name, hist_name, title, x_title, y_title)
    
    input_file = "output/SingleMuon_2017B_Legacy_MoreEvents_ClusterSize2_NumberClusters2000_AllClusters.root"
    plot_name  = "plots/SingleMuon_2017B_Legacy_MoreEvents_ClusterSize2_NumberClusters2000_AllClusters_occupancy"
    hist_name  = "h2d_occupancy"
    title      = "Single Muon 2017B Occ."
    x_title    = "Ring"
    y_title    = "Ladder"
    plot2D(input_file, plot_name, hist_name, title, x_title, y_title)
    
    input_file = "output/SingleMuon_2017B_Legacy_MoreEvents_ClusterSize2_NumberClusters2000_AllClusters.root"
    plot_name  = "plots/SingleMuon_2017B_Legacy_MoreEvents_ClusterSize2_NumberClusters2000_AllClusters_occupancy_cut"
    hist_name  = "h2d_occupancy_cut"
    title      = "Single Muon 2017B Occ. (cut)"
    x_title    = "Ring"
    y_title    = "Ladder"
    plot2D(input_file, plot_name, hist_name, title, x_title, y_title)
    
    input_file = "output/SingleMuon_2017B_Legacy_MoreEvents_ClusterSize2_NumberClusters2000_AllClusters.root"
    plot_name  = "plots/SingleMuon_2017B_Legacy_MoreEvents_ClusterSize2_NumberClusters2000_AllClusters_occupancy_mask"
    hist_name  = "h2d_occupancy_mask"
    title      = "Single Muon 2017B Occ. (mask)"
    x_title    = "Ring"
    y_title    = "Ladder"
    plot2D(input_file, plot_name, hist_name, title, x_title, y_title)
    
    input_file = "output/ZeroBias_2017B_Legacy_MoreEvents_ClusterSize2_NumberClusters2000_AllClusters.root"
    plot_name  = "plots/ZeroBias_2017B_Legacy_MoreEvents_ClusterSize2_NumberClusters2000_AllClusters_occupancy"
    hist_name  = "h2d_occupancy"
    title      = "Zero Bias 2017B Occ."
    x_title    = "Ring"
    y_title    = "Ladder"
    plot2D(input_file, plot_name, hist_name, title, x_title, y_title)
    
    input_file = "output/ZeroBias_2017B_Legacy_MoreEvents_ClusterSize2_NumberClusters2000_AllClusters.root"
    plot_name  = "plots/ZeroBias_2017B_Legacy_MoreEvents_ClusterSize2_NumberClusters2000_AllClusters_occupancy_cut"
    hist_name  = "h2d_occupancy_cut"
    title      = "Zero Bias 2017B Occ. (cut)"
    x_title    = "Ring"
    y_title    = "Ladder"
    plot2D(input_file, plot_name, hist_name, title, x_title, y_title)

    input_file = "output/ZeroBias_2017B_Legacy_MoreEvents_ClusterSize2_NumberClusters2000_AllClusters.root"
    plot_name  = "plots/ZeroBias_2017B_Legacy_MoreEvents_ClusterSize2_NumberClusters2000_AllClusters_occupancy_mask"
    hist_name  = "h2d_occupancy_mask"
    title      = "Zero Bias 2017B Occ. (mask)"
    x_title    = "Ring"
    y_title    = "Ladder"
    plot2D(input_file, plot_name, hist_name, title, x_title, y_title)

    # weighted
    input_file = "output/TTBar_AllClusters_zsmear_weighted_avg.root"
    plot_name  = "plots/TTBar_AllClusters_zsmear_weighted_avg_occupancy"
    hist_name  = "h2d_occupancy"
    title      = "TTBar Occ. (Weighted Avg.)"
    x_title    = "Ring"
    y_title    = "Ladder"
    plot2D(input_file, plot_name, hist_name, title, x_title, y_title)
    
    input_file = "output/SingleMuon_2017B_Legacy_MoreEvents_ClusterSize2_NumberClusters2000_AllClusters_weighted_avg.root"
    plot_name  = "plots/SingleMuon_2017B_Legacy_MoreEvents_ClusterSize2_NumberClusters2000_AllClusters_weighted_avg_occupancy"
    hist_name  = "h2d_occupancy"
    title      = "Single Muon 2017B Occ. (Weighted Avg.)"
    x_title    = "Ring"
    y_title    = "Ladder"
    plot2D(input_file, plot_name, hist_name, title, x_title, y_title)
    
    input_file = "output/ZeroBias_2017B_Legacy_MoreEvents_ClusterSize2_NumberClusters2000_AllClusters_weighted_avg.root"
    plot_name  = "plots/ZeroBias_2017B_Legacy_MoreEvents_ClusterSize2_NumberClusters2000_AllClusters_weighted_avg_occupancy"
    hist_name  = "h2d_occupancy"
    title      = "Zero Bias 2017B Occ. (Weighted Avg.)"
    x_title    = "Ring"
    y_title    = "Ladder"
    plot2D(input_file, plot_name, hist_name, title, x_title, y_title)

def makePlots2():
    input_file = "output/design_0_ge_2pix.root"
    plot_name  = "plots/design_0_ge_2pix_occupancy"
    title      = "design_0_ge_2pix occupancy"
    x_title    = "Ring"
    y_title    = "Ladder"
    plot(input_file, plot_name, title, x_title, y_title)

def makePlots3():
    # using weighted average
    input_file = "output/ZeroBias_2022F_nFiles1_NoCuts_weighted_avg.root"
    plot_name  = "plots2D/ZeroBias_2022F_nFiles1_NoCuts_weighted_avg"
    hist_name  = "h2d_occupancy"
    title      = "Zero Bias 2022F Occ."
    x_title    = "Ring"
    y_title    = "Ladder"
    plot2D(input_file, plot_name, hist_name, title, x_title, y_title)
    
    input_file = "output/ZeroBias_2022F_nFiles1_NoCuts_weighted_avg.root"
    plot_name  = "plots2D/ZeroBias_2022F_nFiles1_NoCuts_weighted_avg_cut"
    hist_name  = "h2d_occupancy_cut"
    title      = "Zero Bias 2022F Occ. (after cut)"
    x_title    = "Ring"
    y_title    = "Ladder"
    plot2D(input_file, plot_name, hist_name, title, x_title, y_title)
    
    input_file = "output/ZeroBias_2022F_nFiles1_ClustSize2_nClust2000_weighted_avg.root"
    plot_name  = "plots2D/ZeroBias_2022F_nFiles1_ClustSize2_nClust2000_weighted_avg"
    hist_name  = "h2d_occupancy"
    title      = "Zero Bias 2022F Occ."
    x_title    = "Ring"
    y_title    = "Ladder"
    plot2D(input_file, plot_name, hist_name, title, x_title, y_title)
    
    input_file = "output/ZeroBias_2022F_nFiles1_ClustSize2_nClust2000_weighted_avg.root"
    plot_name  = "plots2D/ZeroBias_2022F_nFiles1_ClustSize2_nClust2000_weighted_avg_cut"
    hist_name  = "h2d_occupancy_cut"
    title      = "Zero Bias 2022F Occ. (after cut)"
    x_title    = "Ring"
    y_title    = "Ladder"
    plot2D(input_file, plot_name, hist_name, title, x_title, y_title)

def main():
    #makePlots()
    #makePlots2()
    makePlots3()

if __name__ == "__main__":
    main()
 
