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
    c1.Close()

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

    # map for plot info
    plot_info = {}

    plot_info["plot_1"] = {}
    plot_info["plot_1"]["input_file"] = "output/ZeroBias_2022F_nFiles1_NoCuts_weighted_avg.root"
    plot_info["plot_1"]["plot_name"]  = "plots2D/ZeroBias_2022F_nFiles1_NoCuts_weighted_avg"
    plot_info["plot_1"]["hist_name"]  = "h2d_occupancy"
    plot_info["plot_1"]["title"]      = "Zero Bias 2022F Occ."
    plot_info["plot_1"]["x_title"]    = "Ring"
    plot_info["plot_1"]["y_title"]    = "Ladder"
    
    plot_info["plot_2"] = {}
    plot_info["plot_2"]["input_file"] = "output/ZeroBias_2022F_nFiles1_NoCuts_weighted_avg.root"
    plot_info["plot_2"]["plot_name"]  = "plots2D/ZeroBias_2022F_nFiles1_NoCuts_weighted_avg_cut"
    plot_info["plot_2"]["hist_name"]  = "h2d_occupancy_cut"
    plot_info["plot_2"]["title"]      = "Zero Bias 2022F Occ. (after cut)"
    plot_info["plot_2"]["x_title"]    = "Ring"
    plot_info["plot_2"]["y_title"]    = "Ladder"
    
    plot_info["plot_3"] = {}
    plot_info["plot_3"]["input_file"] = "output/ZeroBias_2022F_nFiles1_ClustSize2_nClust2000_weighted_avg.root"
    plot_info["plot_3"]["plot_name"]  = "plots2D/ZeroBias_2022F_nFiles1_ClustSize2_nClust2000_weighted_avg"
    plot_info["plot_3"]["hist_name"]  = "h2d_occupancy"
    plot_info["plot_3"]["title"]      = "Zero Bias 2022F Occ."
    plot_info["plot_3"]["x_title"]    = "Ring"
    plot_info["plot_3"]["y_title"]    = "Ladder"
    
    plot_info["plot_4"] = {}
    plot_info["plot_4"]["input_file"] = "output/ZeroBias_2022F_nFiles1_ClustSize2_nClust2000_weighted_avg.root"
    plot_info["plot_4"]["plot_name"]  = "plots2D/ZeroBias_2022F_nFiles1_ClustSize2_nClust2000_weighted_avg_cut"
    plot_info["plot_4"]["hist_name"]  = "h2d_occupancy_cut"
    plot_info["plot_4"]["title"]      = "Zero Bias 2022F Occ. (after cut)"
    plot_info["plot_4"]["x_title"]    = "Ring"
    plot_info["plot_4"]["y_title"]    = "Ladder"
    
    plot_info["plot_5"] = {}
    plot_info["plot_5"]["input_file"] = "output/ZeroBias_2022F_nFiles10_NoCuts_weighted_avg.root"
    plot_info["plot_5"]["plot_name"]  = "plots2D/ZeroBias_2022F_nFiles10_NoCuts_weighted_avg"
    plot_info["plot_5"]["hist_name"]  = "h2d_occupancy"
    plot_info["plot_5"]["title"]      = "Zero Bias 2022F Occ."
    plot_info["plot_5"]["x_title"]    = "Ring"
    plot_info["plot_5"]["y_title"]    = "Ladder"
    
    plot_info["plot_6"] = {}
    plot_info["plot_6"]["input_file"] = "output/ZeroBias_2022F_nFiles10_NoCuts_weighted_avg.root"
    plot_info["plot_6"]["plot_name"]  = "plots2D/ZeroBias_2022F_nFiles10_NoCuts_weighted_avg_cut"
    plot_info["plot_6"]["hist_name"]  = "h2d_occupancy_cut"
    plot_info["plot_6"]["title"]      = "Zero Bias 2022F Occ. (after cut)"
    plot_info["plot_6"]["x_title"]    = "Ring"
    plot_info["plot_6"]["y_title"]    = "Ladder"
    
    plot_info["plot_7"] = {}
    plot_info["plot_7"]["input_file"] = "output/ZeroBias_2022F_nFiles10_ClustSize2_nClust2000_weighted_avg.root"
    plot_info["plot_7"]["plot_name"]  = "plots2D/ZeroBias_2022F_nFiles10_ClustSize2_nClust2000_weighted_avg"
    plot_info["plot_7"]["hist_name"]  = "h2d_occupancy"
    plot_info["plot_7"]["title"]      = "Zero Bias 2022F Occ."
    plot_info["plot_7"]["x_title"]    = "Ring"
    plot_info["plot_7"]["y_title"]    = "Ladder"
    
    plot_info["plot_8"] = {}
    plot_info["plot_8"]["input_file"] = "output/ZeroBias_2022F_nFiles10_ClustSize2_nClust2000_weighted_avg.root"
    plot_info["plot_8"]["plot_name"]  = "plots2D/ZeroBias_2022F_nFiles10_ClustSize2_nClust2000_weighted_avg_cut"
    plot_info["plot_8"]["hist_name"]  = "h2d_occupancy_cut"
    plot_info["plot_8"]["title"]      = "Zero Bias 2022F Occ. (after cut)"
    plot_info["plot_8"]["x_title"]    = "Ring"
    plot_info["plot_8"]["y_title"]    = "Ladder"

    for key in plot_info:
        print("Creating '{0}'.".format(key))
        input_file  = plot_info[key]["input_file"]
        plot_name   = plot_info[key]["plot_name"]
        hist_name   = plot_info[key]["hist_name"]
        title       = plot_info[key]["title"]
        x_title     = plot_info[key]["x_title"]
        y_title     = plot_info[key]["y_title"]
        plot2D(input_file, plot_name, hist_name, title, x_title, y_title)


def main():
    #makePlots()
    #makePlots2()
    makePlots3()

if __name__ == "__main__":
    main()
 
