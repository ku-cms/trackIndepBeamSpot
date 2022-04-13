#phi_fit2.py

import ROOT
import iminuit
import csv
import os
import numpy as np
import matplotlib.pyplot as plt
from iminuit import Minuit
from iminuit.cost import LeastSquares
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

def fit(input_file, input_dir, plot_dir):

    print("Fitting for input file {0}".format(input_file))

    makeDir(plot_dir)

    g = []
    c = []    
    
    a = ROOT.TFile("{0}/{1}.root".format(input_dir, input_file), "READ")

    tag = "postcut"
    base_name = plot_dir + "/" + input_file + "_" + tag

    num_phi =  a.Get("gr_num_phi_per_ring")
    if (num_phi): 
        print("PASS: loaded num_phi")
    else:          
        print("FAIL: did not load num_phi")

    #pp = PdfPages("{0}.pdf".format(input_file))

    csv_output_name = "{0}_parameters.csv".format(input_file)
    output_column_titles = ["ring", "nop", "amp", "amp error", "shift", "shift error", "offset", "offset error", "chi^2", "dof", "ssv index" , "ssv"]
    with open(csv_output_name, 'w', newline='') as output_csv:
        output_writer = csv.writer(output_csv)
        output_writer.writerow(output_column_titles)

        for i in range(64):

            ring_name = "gr_phi_occ_ring_{0}_{1}".format(tag, i)
            file_ring_name = "{0} gr_phi_occ_ring_{1}_{2}".format(input_file, tag, i)

            g.append (a.Get(ring_name))
            c.append (ROOT.TCanvas(file_ring_name, file_ring_name))

            n = g[-1].GetN()

            # cut on number of points
            #if n < 9:
                #continue

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
            h_y = 250000
            c_start = np.mean(data_y)
            
            least_squares = LeastSquares(data_x, data_y, data_yerr, sin)

            ssv = [-3, -2, -1, 0, 1, 2, 3]
            for j,v in enumerate(ssv):
                m = Minuit(least_squares, a=5000, b=v, c=c_start, limit_a=[0,20000], limit_b=[-np.pi,np.pi])
            
                m.migrad() 
                m.hesse()
                #print(m.values)
                #print(m.errors)
                #print(m.get_param_states())
                #print(m.get_fmin())

                a_fit = m.values["a"]
                b_fit = m.values["b"]
                c_fit = m.values["c"]
                a_fit_errors = m.errors["a"]
                b_fit_errors = m.errors["b"]
                c_fit_errors = m.errors["c"]

                ra_fit = float('%.2g' % m.values["a"])
                rb_fit = float('%.2g' % m.values["b"])
                rc_fit = float('%.2g' % m.values["c"])
                ra_fit_errors = float('%.2g' % m.errors["a"])
                rb_fit_errors = float('%.2g' % m.errors["b"])
                rc_fit_errors = float('%.2g' % m.errors["c"])
                
                fval = round(m.fval,1)
                dof = len(data_x)-m.nfit

                row = [i, n, ra_fit, ra_fit_errors, rb_fit, rb_fit_errors, rc_fit, rc_fit_errors, fval, dof, j, v]
                output_writer.writerow(row)
                
                """
                a_fit_list = np.array(n*[a_fit])
                b_fit_list = np.array(n*[b_fit])
                c_fit_list = np.array(n*[c_fit])

                xf = np.arange(-np.pi, np.pi, 0.1)
                nf = len(xf)
                an = np.array(nf*[a_fit])
                bn = np.array(nf*[b_fit])
                cn = np.array(nf*[c_fit])
                yf = sin(xf, an, bn, cn)
            
                plt.errorbar(data_x, data_y, data_yerr, fmt="o", label="data")
                plt.plot(xf, yf, label="fit")
                #plt.plot(xf, sin(xf, a_fit_list, b_fit_list, c_fit_list), label="fit")
                #plt.text(-2, h_y, "a={0:.2f}\nb={1:.2f}\nc={2:.2f}".format(a_fit, b_fit, c_fit))

                fit_info = [f"$\\chi^2$ / $n_\\mathrm{{dof}}$ = {m.fval:.1f} / {len(data_x) - m.nfit}"]
                fit_info.append("${0} = {1} \\pm {2}$".format("amp", ra_fit, ra_fit_errors))
                fit_info.append("${0} = {1} \\pm {2}$".format("shift", rb_fit, rb_fit_errors))
                fit_info.append("${0} = {1} \\pm {2}$".format("offset", rc_fit, rc_fit_errors))
                plt.legend(title="\n".join(fit_info))

                plt.title("Ring {0}".format(i))
                plt.xlabel("phi")
                plt.ylabel("occupancy")
                plt.xlim([-np.pi, np.pi])
                plt.ylim([0, 20000])
                #plt.savefig("plots/ring{0}".format(i))
                pp.savefig()
                plt.clf()
                #print("n={0}".format(n))
                #print("x={0}".format(data_x))
                #print("y={0}".format(data_y))
                """
        
    #pp.close()

def sin(x, a, b, c):
    return a*np.sin(x-b)+c

def line(x, a, b):
    return a + x * b

def loop():

    input_dir = "output"
    plot_dir = "phi_fit_plots"
    input_files = ["SingleMuon_2017B_Legacy_MoreEvents_ClusterSize2_NumberClusters2000_AllClusters", "TTBar_AllClusters_zsmear", "ZeroBias_2017B_Legacy_MoreEvents_ClusterSize2_NumberClusters2000_AllClusters"]

    for i in range(len(input_files)):
        fit(input_files[i], input_dir, plot_dir)

def main():

    print("Run phi fits")
    loop()

main()
