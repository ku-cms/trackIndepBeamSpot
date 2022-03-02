#phi_fit2.py

import ROOT
import os
import iminuit
import numpy as np
import matplotlib.pyplot as plt
from iminuit import Minuit
from iminuit.cost import LeastSquares
import string
from matplotlib.backends.backend_pdf import PdfPages

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

    pi = 3.14159265

    g = []
    c = []    
    
    a = ROOT.TFile("{0}/{1}.root".format(input_dir, input_file), "READ")
    #f = ROOT.TF1("f", "[0]*sin(x - [1]) + [2]", -pi, pi)

    #f.SetParNames("amp", "shift", "offset")
    #f.SetParLimits(0, 0, 1e10)
    #f.SetParLimits(1, -pi, pi)

    tag = "postcut"

    base_name = plot_dir + "/" + input_file + "_" + tag

    num_phi =  a.Get("gr_num_phi_per_ring")
    if (num_phi): 
        print("PASS: loaded num_phi")
    else:          
        print("FAIL: did not load num_phi")

    pp = PdfPages("{0}2.pdf".format(input_file))

    for i in range(64):

        #f.SetParameters(4386, 1.484, 1.132e5)

        ring_name = "gr_phi_occ_ring_{0}_{1}".format(tag, i)
        file_ring_name = "{0} gr_phi_occ_ring_{1}_{2}".format(input_file, tag, i)
        
        g.append (a.Get(ring_name))
        c.append (ROOT.TCanvas(file_ring_name, file_ring_name))

        n = g[-1].GetN()
        data_x = []
        data_y = []
        data_yerr = []
        
        for j in range(n):
            data_x.append(g[-1].GetPointX(j))
            data_y.append(g[-1].GetPointY(j))
            data_yerr.append(np.sqrt(g[-1].GetPointY(j)))

        min_y = np.min(data_y)
        max_y = np.max(data_y)
        delta_y = max_y - min_y
        #h_y = 0.8*delta_y + min_y
        h_y = 250000
        c_start = np.mean(data_y)
        least_squares = LeastSquares(data_x, data_y, data_yerr, sin)
        m = Minuit(least_squares, a=5000, b=0, c=c_start, limit_a=[0,20000], limit_b=[-pi,pi])
        #m = Minuit(least_squares, a=0, b=0)

        m.migrad() 
        m.hesse()
        print(m.values)
        #print(m.values["error_b"])
        print(m.get_param_states())
        print(m.get_fmin())
        
        #print(m.values["a"])
        a_fit = m.values["a"]
        b_fit = m.values["b"]
        c_fit = m.values["c"]
 
        a_fit_list = np.array(n*[a_fit])
        b_fit_list = np.array(n*[b_fit])
        c_fit_list = np.array(n*[c_fit])

        xf = np.arange(-pi, pi, 0.1)
        nf = len(xf)
        an = np.array(nf*[a_fit])
        bn = np.array(nf*[b_fit])
        cn = np.array(nf*[c_fit])

        yf = sin(xf, an, bn, cn)

        plt.errorbar(data_x, data_y, data_yerr, fmt="o", label="data")
        plt.plot(xf, yf, label="fit")
        #plt.plot(xf, sin(xf, a_fit_list, b_fit_list, c_fit_list), label="fit")
        plt.text(-2, h_y, "a={0:.2f}\nb={1:.2f}\nc={2:.2f}".format(a_fit, b_fit, c_fit))
        plt.legend()
        plt.title("Ring {0}".format(i))
        plt.xlabel("phi")
        plt.ylabel("occupancy")
        plt.xlim([-pi, pi])
        plt.ylim([0, 300000])
        #plt.savefig("plots/ring{0}".format(i))
        pp.savefig()
        plt.clf()
        #print("n={0}".format(n))
        #print("x={0}".format(data_x))
        #print("y={0}".format(data_y))
        #g[-1].Fit(f, "R")
        #g[-1].SetTitle(ring_name)
        
        #g[-1].SetMarkerStyle(20)
        #g[-1].GetYaxis().SetRangeUser(y_min, y_max)
   
        #c[-1].cd()
        #g[-1].Draw("AP")
        
        
    pp.close()
'''        
    c[0].Print("{0}.pdf(".format(base_name))
    for i in range(1,63):
        c[i].Print("{0}.pdf".format(base_name))
    c[63].Print("{0}.pdf)".format(base_name))
'''

def sin(x, a, b, c):
    return a*np.sin(x-b)+c

def line(x, a, b):
    return a + x * b

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
