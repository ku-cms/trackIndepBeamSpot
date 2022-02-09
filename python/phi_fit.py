#phi_fit.py

import ROOT
import os

def makeDir(dir_name):
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)
# Make sure ROOT.TFile.Open(fileURL) does not seg fault when $ is in sys.argv (e.g. $ passed in as argument)
ROOT.PyConfig.IgnoreCommandLineOptions = True
# Make plots faster without displaying them
ROOT.gROOT.SetBatch(ROOT.kTRUE)
# Tell ROOT not to be in charge of memory, fix issue of histograms being deleted when ROOT file is closed:
ROOT.TH1.AddDirectory(False)

def fit(input_file, input_dir, plot_dir, y_min, y_max):

    print("Fitting for input file {0}".format(input_file))

    makeDir(plot_dir)

    ROOT.gStyle.SetOptFit(111)
    ROOT.gStyle.SetStatW(0.1)                
    ROOT.gStyle.SetStatH(0.1)
    ROOT.gStyle.SetStatY(0.99)
    ROOT.gStyle.SetTitleFontSize(0.04)

    pi = 3.14159265

    chisq = []
    amp = []
    shift = []
    offset = []
    amp_err = []
    shift_err = []
    offset_err = []

    g = []
    c = []
    
    h_num_phi     = ROOT.TH1F ("h_num_phi",     "Number of phi points", 64, 0, 64)
    h_chisq       = ROOT.TH1F ("h_chisq",       "Fit chi squares", 64, 0, 64)
    h_amp         = ROOT.TH1F ("h_amp",         "Fit amplitudes", 64, 0, 64)
    h_shift       = ROOT.TH1F ("h_shift",       "Fit phi shifts", 64, 0, 64)
    h_offset      = ROOT.TH1F ("h_offset",      "Fit offsets", 64, 0, 64)
    h_num_phi_cut = ROOT.TH1F ("h_num_phi_cut", "Number of phi points after cut", 64, 0, 64)
    h_chisq_cut   = ROOT.TH1F ("h_chisq_cut",   "Fit chi squares after cut", 64, 0, 64)
    h_amp_cut     = ROOT.TH1F ("h_amp_cut",     "Fit amplitudes after cut", 64, 0, 64)
    h_shift_cut   = ROOT.TH1F ("h_shift_cut",   "Fit phi shifts after cut", 64, 0, 64)
    h_offset_cut  = ROOT.TH1F ("h_offset_cut",  "Fit offsets after cut", 64, 0, 64)

    a = ROOT.TFile("{0}/{1}.root".format(input_dir, input_file), "READ")
    f = ROOT.TF1("f", "[0]*sin(x - [1]) + [2]", -pi, pi)

    f.SetParNames("amp", "shift", "offset")
    f.SetParLimits(0, 0, 1e10)
    f.SetParLimits(1, -pi, pi)

    tag = "postcut"

    base_name = plot_dir + "/" + input_file + "_" + tag

    num_phi =  a.Get("gr_num_phi_per_ring")
    if (num_phi): 
        print("PASS: loaded num_phi")
    else:          
        print("FAIL: did not load num_phi")

    for i in range(64):

        f.SetParameters(4386, 1.484, 1.132e5)

        ring_name = "gr_phi_occ_ring_{0}_{1}".format(tag, i)
        file_ring_name = "{0} gr_phi_occ_ring_{1}_{2}".format(input_file, tag, i)
        
        g.append (a.Get(ring_name))
        c.append (ROOT.TCanvas(file_ring_name, file_ring_name))
        g[-1].Fit(f, "R")
        g[-1].SetTitle(ring_name)
        
        g[-1].SetMarkerStyle(20)
        g[-1].GetYaxis().SetRangeUser(y_min, y_max)
   
        c[-1].cd()
        g[-1].Draw("AP")

    c[0].Print("{0}.pdf(".format(base_name))
    for i in range(1,63):
        c[i].Print("{0}.pdf".format(base_name))
    c[63].Print("{0}.pdf)".format(base_name))

def loop():

    input_dir = "output"
    plot_dir = "phi_fit_plots"
    input_files = []
    y_min_vals = []
    y_max_vals = []

    input_files.append("SingleMuon_2017B_Legacy_MoreEvents_ClusterSize2_NumberClusters2000_AllClusters")

    y_min_vals.append(0.0)  
    y_max_vals.append(300000.0)

    for i in range(len(input_files)):
        fit(input_files[i], input_dir, plot_dir, y_min_vals[i], y_max_vals[i])

def main():

    print("Run phi fits")
    loop()

main()



