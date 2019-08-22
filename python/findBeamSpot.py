#!/usr/bin/env python

import ROOT as rt
import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict
from array import array


rt.gROOT.SetBatch()
rt.gStyle.SetOptStat(111)
rt.gStyle.SetOptFit(1111)
rt.TH1.AddDirectory(rt.kFALSE)


def preprocess_hist(hist):
    """
    Function to remove outlier bins from histogram
    :param hist: TH1 to be processed
    :return hist_clean: processed Histogram
    """
    ##  for real data
    # good_bins = [[76, 108],
    #              [122, 155],
    #              [255, 288],
    #              [301, 332],
    #              [416, 446],
    #              [460, 491],
    #              [591, 624],
    #              [638, 671],
    #              [751, 778],
    #              [791, 826],
    #              [917, 950],
    #              [967, 999]]
    ## for realistic MC
    good_bins = [[81, 118],
                 [128, 165],
                 [249, 285],
                 [295, 331],
                 [415, 451],
                 [461, 498],
                 [581, 618],
                 [628, 665],
                 [748, 785],
                 [795, 831],
                 [915, 951],
                 [961, 998]]
    # bad_bins_mods = [[]]
    nbins = hist.GetNbinsX()
    name = hist.GetName()
    hist_clean = rt.TH1F(name + '_cleaned', 'Occ vs phi', 1000, -rt.TMath.Pi(), rt.TMath.Pi())
    for ibin in xrange(nbins):
        binc = hist.GetBinContent(ibin)
        for bin_low, bin_high in good_bins:
            if bin_low < ibin < bin_high:
                hist_clean.SetBinContent(ibin, binc)
            else:
                continue
    return hist_clean


def make_xy_graph(hist=None, graph=None):
    """
    Function to convert TH1 in phi into TGraph in X and Y. Writes out graph to a TFile
    :param hist: TH1 histogram
    :return gr_new: TGraph
    """

    if hist is not None:
        x = array('d')
        y = array('d')
        gr_array = []
        nbins = hist.GetNbinsX()
        hist_name = hist.GetName()
        for ibin in xrange(nbins):
            phi = hist.GetBinCenter(ibin)
            r = hist.GetBinContent(ibin)
            if r == 0: continue
            if 'r_phi' in hist_name:
                x.append(r * np.cos(phi))
                y.append(r * np.sin(phi))
            else:
                # x.append(r * np.cos(phi))
                # y.append(r * np.sin(phi))
                x.append((1. / np.sqrt(r)) * np.cos(phi))
                y.append((1. / np.sqrt(r)) * np.sin(phi))
            gr_array.append((x, y))
            n = len(gr_array)
            gr_new = rt.TGraph(n, x, y)
            gr_new.SetName(hist.GetName()+"_xy_gr")
    elif graph is not None:
        phi = graph.GetX()
        r = graph.GetY()
        x = array('d')
        y = array('d')
        n = len(r)
        for i in xrange(n):
            x.append(r[i] * np.cos(phi[i]))
            y.append(r[i] * np.sin(phi[i]))
        gr_new = rt.TGraph(n, x, y)
        gr_new.SetName(graph.GetName()+"_xy_gr")
    # out_f = rt.TFile(hist.GetName()+'_gr', 'w')
    # out_f.cd()
    # gr_new.Write()
    # out_f.Close()

    return gr_new


def fit_ellipse(gr, data_type, mod):


    mod_start = {
        '-4' : 4.8,
        '-3' : 4.8,
        '-2' : 3.4,
        '-1' : 1.8,
        '1' : 1.8,
        '2' : 3.4,
        '3' : 4.8,
        '4' : 4.8,
    }
    # mod_start = {
    #    '-4' : 0.05,
    #    '-3' : 0.04,
    #    '-2' : 0.1,
    #    '-1' : 0.3,
    #    '1' : 0.3,
    #    '2' : 0.1,
    #    '3' : 0.04,
    #    '4' : 0.05,
    #}
    def func_ellipse(npar, gin, f, par, iflag):
        np = gr.GetN()
        chi = 0.
        x = gr.GetX()
        y = gr.GetY()

        for i in xrange(np):
            u = (x[i] - par[0])
            v = (y[i] - par[1])
            rot_sin = rt.TMath.Sin(par[4])
            rot_cos = rt.TMath.Cos(par[4])
            r = rt.TMath.Sqrt((u*rot_cos-v*rot_sin)**2/par[2]**2 + (u*rot_sin+v*rot_cos)**2/par[3]**2)

            chi += abs(r - 1.)
        f[0] = chi
    # func = rt.TF1('func', func_ellipse)
    # wf = rt.Math.WrappedTF1(func)

    fitter = rt.TVirtualFitter.Fitter(0, 5)

    fitter.SetFCN(func_ellipse)
    fitter.SetParameter(0, 'x0', 0.001, 0.001, -0.2, 0.2)
    fitter.SetParameter(1, 'y0', 0.001, 0.001, -0.2, 0.2)
    fitter.SetParameter(2, 'r0_x', mod_start[str(mod)], 0.01, 0.0, 5.)
    fitter.SetParameter(3, 'r0_y', mod_start[str(mod)], 0.01, 0.0, 5.)
    fitter.SetParameter(4, 'a', 0.761, 0.01, 0.0, 3.0)
    # fitter.SetParameter(0, 'x0', 0.0973, 0.001, -1.0, 1.0)
    # fitter.SetParameter(1, 'y0', -0.1015, 0.001, -1.0, 1.0)
    # fitter.SetParameter(2, 'r0_x', 2.733, 0.01, 0.0, 3.5)
    # fitter.SetParameter(3, 'r0_y', 2.788, 0.01, 0.0, 3.5)
    # fitter.SetParameter(4, 'a', 0.434, 0.01, -1.0, 3.0)
    #fitter.SetParameter(0, 'x0', 0.0168, 0.001, -0.1, 0.1)
    #fitter.SetParameter(1, 'y0', -0.0137, 0.001, -0.1, 0.1)
    #fitter.SetParameter(2, 'r0_x', 0.188, 0.01, 0.0, 1.0)
    #fitter.SetParameter(3, 'r0_y', 0.153, 0.01, 0.0, 1.0)
    #fitter.SetParameter(4, 'a', 0.7, 0.001, 0.0, 2.0)

    arg_list_proxy = array('d', [0, 0, 0, 0, 0])
    fitter.ExecuteCommand("SIMPLEX", arg_list_proxy, 0)
    fitter.ExecuteCommand("MIGRAD", arg_list_proxy, 0)

    fitter.PrintResults(1, 0)

    c0 = rt.TCanvas('c0', '', 800, 600)
    gr.Draw('A*')
    gr.SetMarkerStyle(2)
    c0.Update()

    x0 = fitter.GetParameter(0)
    y0 = fitter.GetParameter(1)
    xr = fitter.GetParameter(2)
    yr = fitter.GetParameter(3)
    a = fitter.GetParameter(5)

    el = rt.TEllipse(x0, y0, xr, yr, 0., 360., a*(180./rt.TMath.Pi()))

    el.SetLineColor(rt.kRed)
    el.SetLineWidth(4)
    el.SetFillStyle(0)
    el.Draw("SAME")

    c0.Update()
    gr.SaveAs("plots_4219/"+data_type+"/graph_"+gr.GetName()+"_"+data_type+"_phi_ellipse_21419.root")
    c0.SaveAs("plots_4219/"+data_type+"/fit_"+gr.GetName()+"_"+data_type+"_phi_ellipse_21419.pdf")
    c0.SaveAs("plots_4219/"+data_type+"/fit_"+gr.GetName()+"_"+data_type+"_phi_ellipse_21419.png")
    c0.SaveAs("plots_4219/"+data_type+"/fit_"+gr.GetName()+"_"+data_type+"_phi_ellipse_21419.root")


def quad_func(phi, par):
    quad_list = [[-2.7, -2.1, 1.42, 7.55, 12.79, 0.287],
                 [-1.6, -1.01, 1.465, -1.387, 3.182, -1.798],
                 [-0.6, 0.0, 1.452, 0.6599, 2.931, -0.053],
                 [0.5, 1.15, 1.4, 7.073, 11.65, 3.29],
                 [1.53, 2.10, 1.34, -5.985, 9.299, -0.4137],
                 [2.56, 3.2, 1.34, 1.633, 3.098, 3.468]]
    r = 2.77
    for quad in quad_list:
        phi_min = quad[0]
        phi_max = quad[1]
        if phi_min < phi < phi_max:
            r = quad[2]*(phi[0] - quad[5])**2 + quad[3] * (phi[0] - quad[5]) + quad[4]
        else: continue
    return (r**2/2.77**2) * (par[0] * np.sin(par[1] * phi[0] - par[2])) + par[3]


def fit_func_root_2d(inFile, histName):
    f = rt.TFile(inFile, 'r')

    h = f.Get(histName)

    func = rt.TF1('func', quad_func, -3.0, 3.0, 4)

    func.SetParameters(-0.4, -1.3, 33., 3.)

    func.SetParName(0, 'a')
    func.SetParName(1, 'b')
    func.SetParName(2, 'c')
    func.SetParName(3, 'd')
    #
    # func.SetParLimits(0, 0., 0.9)
    # func.SetParLimits(3, 3.5, 6.)
    #
    fitResult = h.Fit('func', 'SRM')

    fitResult.Print('V')

    c0 = rt.TCanvas('c0', '', 800, 600)

    h.Draw()

    func.Draw('same')

    c0.SaveAs('d_' + histName + '_102318.png')
    c0.SaveAs('d_' + histName + '_102318.root')


def fit_func_root_quad_mods(inFile, histName_template, quadList, mods):
    # f = rt.TFile(inFile, 'r')

   # h = f.Get(histName)
   # h_all = f.Get(histName)


    fit_results_mods = OrderedDict()

    c1 = rt.TCanvas('c0', '', 800, 600)
    for mod in mods:
        f = rt.TFile(inFile, 'r')
        if isinstance(mod, int) and mod < 0:
            histName = histName_template.format('neg_'+str(mod))
        else:
            histName = histName_template.format(str(abs(mod)))
        print histName
        fit_results_mods[mod] = {}
        fit_results_mods[mod][histName] = f.Get(histName)
        fit_results_mods[mod][histName+'_all'] = f.Get(histName)
        fit_results_mods[mod]['results'] = []
        fit_results_mods[mod]['funcs'] = []
        print "##############################"
        print "Fitting for mod: ", mod
        print "##############################"

        for i, quad in enumerate(quadList):
            phi_min = quad[0]
            phi_max = quad[1]
            a0 = quad[2]
            b0 = quad[3]
            c0 = quad[4]
            d0 = quad[5]

            fit_results_mods[mod]['funcs'].append(rt.TF1('func_'+str(mod)+'_'+str(i), '[0]*(x-[3])*(x-[3])+[1]*(x-[3])+[2]', phi_min, phi_max))
            fit_results_mods[mod]['funcs'][i].SetParameter(0, rt.Double(a0))
            fit_results_mods[mod]['funcs'][i].SetParameter(1, rt.Double(b0))
            fit_results_mods[mod]['funcs'][i].SetParameter(2, rt.Double(c0))
            fit_results_mods[mod]['funcs'][i].SetParameter(3, rt.Double(d0))

            fit_results_mods[mod]['funcs'][i].SetParName(0, 'a')
            fit_results_mods[mod]['funcs'][i].SetParName(1, 'b')
            fit_results_mods[mod]['funcs'][i].SetParName(2, 'c')
            fit_results_mods[mod]['funcs'][i].SetParName(3, 'd')
            #
            # func.SetParLimits(0, 0., 0.9)
            # func.SetParLimits(3, 3.5, 6.)
            #
            fit_results_mods[mod]['results'].append(fit_results_mods[mod][histName].Fit('func_'+str(mod)+'_'+str(i), 'SRM'))

            if str(mod) == '1' and i == 0: continue
            if str(mod) == '-1' and i == 3: continue
            fit_results_mods[mod]['results'][i].Print('V')

            fit_results_mods[mod][histName].Draw('A*')

            fit_results_mods[mod]['funcs'][i].Draw('same')

            # c1.SaveAs('d_' + histName + '_' + str(i) + '.png')
            # c1.SaveAs('d_' + histName + '_' + str(i) + '.root')
        rt.gStyle.SetOptFit(0)
        fit_results_mods[mod][histName + '_all'].Draw('A*')
        for i, f in enumerate(fit_results_mods[mod]['funcs']):
            f.Draw('same')
        c1.SaveAs('d_' + histName + '_quad_all' + '.png')
        c1.SaveAs('d_' + histName + '_quad_all' + '.root')
    return fit_results_mods


def fit_func_root_sine_mods(inFile, histName_template, quadList, quad_fit, mods, params):

    sine_fit_results = OrderedDict()

    for mod in mods:
        f = rt.TFile(inFile, 'r')
        if isinstance(mod, int) and mod < 0:
            hist_name = histName_template.format('neg_'+str(abs(mod)))
        else:
            hist_name = histName_template.format(str(abs(mod)))
        sine_fit_results[mod] = {}
        tmp_hist = f.Get(hist_name)
        sine_fit_results[mod]['hist'] = preprocess_hist(tmp_hist)

        new_quad_list = []
        for i, fit in enumerate(quad_fit[mod]['funcs']):
            if str(mod) == '1' and i == 0: continue
            if str(mod) == '-1' and i == 3: continue
            phi_min = quadList[i][0]
            phi_max = quadList[i][1]

            a_tmp = fit.GetParameter(0)
            b_tmp = fit.GetParameter(1)
            c_tmp = fit.GetParameter(2)
            d_tmp = fit.GetParameter(3)

            new_quad_list.append((phi_min, phi_max, a_tmp, b_tmp, c_tmp, d_tmp))

        def new_quad_func(phi, par):
            r = 2.77
            for quad in new_quad_list:
                phi_min = quad[0]
                phi_max = quad[1]
                if phi_min < phi[0] < phi_max:
                    r = quad[2]*(phi[0] - quad[5])**2 + quad[3] * (phi[0] - quad[5]) + quad[4]
                else: continue
            return (r**2/2.77**2) * (par[0] * np.sin(par[1] * phi[0] - par[2])) + par[3]

        sine_fit_results[mod]['func'] = rt.TF1('func_sine_'+str(mod), new_quad_func, -3.0, 3.0)

        sine_fit_results[mod]['func'].SetParameter(0, rt.Double(params[0]))
        sine_fit_results[mod]['func'].SetParameter(1, rt.Double(params[1]))
        sine_fit_results[mod]['func'].SetParameter(2, rt.Double(params[2]))
        sine_fit_results[mod]['func'].SetParameter(3, rt.Double(params[3]))

        sine_fit_results[mod]['func'].SetParName(0, 'a')
        sine_fit_results[mod]['func'].SetParName(1, 'b')
        sine_fit_results[mod]['func'].SetParName(2, 'c')
        sine_fit_results[mod]['func'].SetParName(3, 'd')
        #
        # func.SetParLimits(0, 0., 0.9)
        # func.SetParLimits(3, 3.5, 6.)
        #
        sine_fit_results[mod]['result'] = sine_fit_results[mod]['hist'].Fit('func_sine_'+str(mod), 'S')

        sine_fit_results[mod]['result'].Print('V')

        c0 = rt.TCanvas('c0', '', 800, 600)

        sine_fit_results[mod]['hist'].Draw()

        sine_fit_results[mod]['func'].Draw('same')

        c0.SaveAs('d_' + hist_name + '_sine_171218.png')
        c0.SaveAs('d_' + hist_name + '_sine_171218.root')


def fit_func_root_quad(inFile, histName, quadList):
    f = rt.TFile(inFile, 'r')

    h = f.Get(histName)
    h_all = f.Get(histName)

    funcs = []
    fit_results = []
    for i, quad in enumerate(quadList):

        phi_min = quad[0]
        phi_max = quad[1]
        a0 = quad[2]
        b0 = quad[3]
        c0 = quad[4]
        d0 = quad[5]
        c1 = rt.TCanvas('c0', '', 800, 600)

        funcs.append(rt.TF1('func_'+str(i), '[0]*(x-[3])*(x-[3])+[1]*(x-[3])+[2]', phi_min, phi_max))
        funcs[i].SetParameter(0, rt.Double(a0))
        funcs[i].SetParameter(1, rt.Double(b0))
        funcs[i].SetParameter(2, rt.Double(c0))
        funcs[i].SetParameter(3, rt.Double(d0))

        funcs[i].SetParName(0, 'a')
        funcs[i].SetParName(1, 'b')
        funcs[i].SetParName(2, 'c')
        funcs[i].SetParName(3, 'd')
        #
        # func.SetParLimits(0, 0., 0.9)
        # func.SetParLimits(3, 3.5, 6.)
        #
        fit_results.append(h.Fit('func_'+str(i), 'SRM'))

        fit_results[i].Print('V')


        h.Draw('A*')

        funcs[i].Draw('same')

        c1.SaveAs('d_' + histName + '_' + str(i) + '.png')
        c1.SaveAs('d_' + histName + '_' + str(i) + '.root')
    rt.gStyle.SetOptFit(0)
    c1 = rt.TCanvas('c0', '', 800, 600)
    h_all.Draw('A*')
    for i, f in enumerate(funcs):
        f.Draw('same')
    c1.SaveAs('d_' + histName + '_quad_all' + '.png')
    c1.SaveAs('d_' + histName + '_quad_all' + '.root')


def fit_func_root_sine(inFile, histName, a0, b0, c0, d0):
    f = rt.TFile(inFile, 'r')

    h = f.Get(histName)
    func = rt.TF1('func', '[0]*TMath::Sin([1]*x-[2])+[3]', -3.0, 3.0)

    func.SetParameter(0, rt.Double(a0))
    func.SetParameter(1, rt.Double(b0))
    func.SetParameter(2, rt.Double(c0))
    func.SetParameter(3, rt.Double(d0))

    func.SetParName(0, 'a')
    func.SetParName(1, 'b')
    func.SetParName(2, 'c')
    func.SetParName(3, 'd')
    #
    # func.SetParLimits(0, 0., 0.9)
    # func.SetParLimits(3, 3.5, 6.)
    #
    fitResult = h.Fit('func', 'SRM')

    fitResult.Print('V')

    c0 = rt.TCanvas('c0', '', 800, 600)

    h.Draw()

    func.Draw('same')

    c0.SaveAs('d_' + histName + '_sine.png')
    c0.SaveAs('d_' + histName + '_sine.root')


def sineFunc(x, a, b, c, d):
    return a * np.sin(x*b - c) + d


def plot_circle(a, b, c, d, name):
    x = np.arange(0, 2 * np.pi + np.pi / 64, np.pi / 32)
    y = sineFunc(x, a, b, c, d)
    plt.polar(x, y)
    plt.polar([np.sqrt(0.096178 ** 2 + 0.063320 ** 2)], [np.arctan(0.096178 / -0.063320)])
    plt.savefig(name + '.png')
    plt.close()


def plot_circle_func(a, b, c, d, name):
    x = np.arange(0, 2 * np.pi + np.pi / 64, np.pi / 32)
    y = quad_func([x], [a, b, c, d])
    plt.polar(x, y)
    # plt.polar([np.sqrt(0.096178 ** 2 + 0.063320 ** 2)], [np.arctan(0.096178 / -0.063320)])
    plt.savefig(name + '.png')
    plt.close()


def make_occ_table(inner_hist_, outer_hist_):

    phi_quadrants = [
        (0., np.pi/2., '+x+y'),
        (np.pi/2., np.pi, '-x+y'),
        (np.pi, (3*np.pi)/2., '-x-y'),
        ((3*np.pi)/2., 2*np.pi, '+x-y'),
    ]

    inner_phi_modules = {
        ('-5', -np.pi, -1.86),
        ('-3', -1.86, -0.8),
        ('-1', -0.8, 0.2),
        ('2', 0.2, 1.3),
        ('4', 1.3, 2.3),
        ('6', 2.3, np.pi),
    }

    outer_phi_modules = {
        ('-6', -np.pi, -2.37),
        ('-4', -2.37, -1.33),
        ('-2', -1.33, -0.33),
        ('1', -0.33, 0.80),
        ('3', 0.80, 1.85),
        ('5', 1.85, np.pi),
    }
    out_table = OrderedDict([
        ('-5' , OrderedDict([('+x+y' , 0.),
                            ('+x-y' , 0.),
                            ('-x+y' , 0.),
                            ('-x-y' , 0.),])),
        ('-3' , OrderedDict([('+x+y' , 0.),
                             ('+x-y' , 0.),
                             ('-x+y' , 0.),
                             ('-x-y' , 0.),])),
        ('-1' , OrderedDict([('+x+y' , 0.),
                             ('+x-y' , 0.),
                             ('-x+y' , 0.),
                             ('-x-y' , 0.),])),
        ('2' , OrderedDict([('+x+y' , 0.),
                            ('+x-y' , 0.),
                            ('-x+y' , 0.),
                            ('-x-y' , 0.),])),
        ('4' , OrderedDict([('+x+y' , 0.),
                            ('+x-y' , 0.),
                            ('-x+y' , 0.),
                            ('-x-y' , 0.),])),
        ('6' , OrderedDict([('+x+y' , 0.),
                            ('+x-y' , 0.),
                            ('-x+y' , 0.),
                            ('-x-y' , 0.),])),
        ('-6' , OrderedDict([('+x+y' , 0.),
                             ('+x-y' , 0.),
                             ('-x+y' , 0.),
                             ('-x-y' , 0.),])),
        ('-4' , OrderedDict([('+x+y' , 0.),
                             ('+x-y' , 0.),
                             ('-x+y' , 0.),
                             ('-x-y' , 0.),])),
        ('-2' , OrderedDict([('+x+y' , 0.),
                             ('+x-y' , 0.),
                             ('-x+y' , 0.),
                             ('-x-y' , 0.),])),
        ('1' , OrderedDict([('+x+y' , 0.),
                            ('+x-y' , 0.),
                            ('-x+y' , 0.),
                            ('-x-y' , 0.),])),
        ('3' , OrderedDict([('+x+y' , 0.),
                            ('+x-y' , 0.),
                            ('-x+y' , 0.),
                            ('-x-y' , 0.),])),
        ('5' , OrderedDict([('+x+y' , 0.),
                            ('+x-y' , 0.),
                            ('-x+y' , 0.),
                            ('-x-y' , 0.),])),
    ])
    for ibin in xrange(inner_hist_.GetNbinsX()):
        width = inner_hist_.GetBinWidth(ibin)
        low_edge = inner_hist_.GetBinLowEdge(ibin)
        high_edge = low_edge + width
        error = inner_hist_.GetBinError(ibin)
        content = inner_hist_.GetBinContent(ibin)
        if low_edge < 0.:
            pos_low_edge = low_edge + 2*np.pi
        else:
            pos_low_edge = low_edge
        if high_edge < 0.:
            pos_high_edge = high_edge + 2*np.pi
        else:
            pos_high_edge = high_edge

        mod = '0'
        for entry in inner_phi_modules:
            if low_edge >= entry[1] and high_edge < entry[2]:
                mod = entry[0]
            else:
                continue

        if mod == '0':
            continue

        quad = '0'
        low_quad = '0'
        high_quad = '0'
        quad_split = False
        low_div = 1.
        high_div = 1.
        for entry in phi_quadrants:
            if pos_low_edge >= entry[0] and pos_high_edge < entry[1]:
                quad = entry[2]
            elif pos_low_edge < entry[0] <= pos_high_edge:
                high_quad = entry[2]
                quad_split = True
                high_diff = pos_high_edge - entry[0]
                low_diff = entry[0] - pos_low_edge
                print mod, '1st'
                print pos_low_edge, entry[0], pos_high_edge
                print low_diff, high_diff
                print width
                low_div = abs(low_diff / width)
                high_div = abs(high_diff / width)

                for entry_2 in phi_quadrants:
                    if entry_2[0] <= pos_low_edge < entry_2[1]:
                        low_quad = entry_2[2]
                    else:
                        continue
            elif pos_low_edge < entry[1] <= pos_high_edge:
                low_quad = entry[2]
                quad_split = True
                high_diff = pos_high_edge - entry[1]
                low_diff = entry[1] - pos_low_edge
                print mod, '2nd'
                print pos_low_edge, entry[1], pos_high_edge
                print low_diff, high_diff
                print width

                low_div = abs(low_diff / width)
                high_div = abs(high_diff / width)

                for entry_2 in phi_quadrants:
                    if entry_2[0] <= pos_high_edge < entry_2[1]:
                        high_quad = entry_2[2]
                    else:
                        continue
            else:
                continue
        if quad_split is False:
            out_table[mod][quad] += content
        elif quad_split is True:
            print low_div, high_div
            out_table[mod][low_quad] += content * low_div
            out_table[mod][high_quad] += content * high_div

    for ibin in xrange(outer_hist_.GetNbinsX()):
        width = outer_hist_.GetBinWidth(ibin)
        low_edge = outer_hist_.GetBinLowEdge(ibin)
        high_edge = low_edge + width
        error = outer_hist_.GetBinError(ibin)
        content = outer_hist_.GetBinContent(ibin)
        if low_edge < 0.:
            pos_low_edge = low_edge + 2*np.pi
        else:
            pos_low_edge = low_edge
        if high_edge < 0.:
            pos_high_edge = high_edge + 2*np.pi
        else:
            pos_high_edge = high_edge

        mod = '0'
        for entry in outer_phi_modules:
            if low_edge >= entry[1] and high_edge < entry[2]:
                mod = entry[0]
            else:
                continue

        if mod == '0':
            continue
        if mod == '-6':
            continue
        if mod == '1':
            continue
        quad = '0'
        low_quad = '0'
        high_quad = '0'
        quad_split = False
        low_div = 1.
        high_div = 1.
        for entry in phi_quadrants:
            if pos_low_edge >= entry[0] and pos_high_edge < entry[1]:
                quad = entry[2]
            elif pos_low_edge < entry[0] <= pos_high_edge:
                high_quad = entry[2]
                quad_split = True
                high_diff = pos_high_edge - entry[0]
                low_diff = entry[0] - pos_low_edge
                print mod, '1st'
                print pos_low_edge, entry[0], pos_high_edge
                print low_diff, high_diff
                print width
                low_div = abs(low_diff / width)
                high_div = abs(high_diff / width)

                for entry_2 in phi_quadrants:
                    if entry_2[0] <= pos_low_edge < entry_2[1]:
                        low_quad = entry_2[2]
                    else:
                        continue
            elif pos_low_edge < entry[1] <= pos_high_edge:
                low_quad = entry[2]
                quad_split = True
                high_diff = pos_high_edge - entry[1]
                low_diff = entry[1] - pos_low_edge
                print mod, '2nd'
                print pos_low_edge, entry[1], pos_high_edge
                print low_diff, high_diff
                print width

                low_div = abs(low_diff / width)
                high_div = abs(high_diff / width)

                for entry_2 in phi_quadrants:
                    if entry_2[0] <= pos_high_edge < entry_2[1]:
                        high_quad = entry_2[2]
                    else:
                        continue
            else:
                continue
        if quad_split is False:
            out_table[mod][quad] += content
        elif quad_split is True:
            print low_div, high_div
            out_table[mod][low_quad] += content * low_div
            out_table[mod][high_quad] += content * high_div

    return out_table


def make_missing_table(inner_hist_, outer_hist_):
    out_table = OrderedDict([
        ('inner' , OrderedDict([('+x+y' , 0.),
                             ('+x-y' , 0.),
                             ('-x+y' , 0.),
                             ('-x-y' , 0.),])),
        ('outer' , OrderedDict([('+x+y' , 0.),
                             ('+x-y' , 0.),
                             ('-x+y' , 0.),
                             ('-x-y' , 0.),])),
    ])

    quad_map = { '1' : '+x+y',
                 '2' : '-x+y',
                 '3' : '-x-y',
                 '4' : '+x-y'}

    bins = [1, 2, 3, 4]
    total_pix_inner_map = {
        '1' : 36864 + (0.07*36864) + (0.019*36864),
        '2' : 73728 - (0.07*36864) - (0.019*36864),
        '3' : 36864 + (0.07*36864) + (0.019*36864),
        '4' : 73728 - (0.07*36864) - (0.019*36864),
    }
    total_pix_outer_map = {
        '2' : 36864,
        '1' : 36864,
        '4' : 36864,
        '3' : 36864
    }

    for bin in bins:
        in_tmp = inner_hist_.GetBinContent(bin)
        out_tmp = outer_hist_.GetBinContent(bin)
        pix_total_inner = total_pix_inner_map[str(bin)]
        pix_total_outer = total_pix_outer_map[str(bin)]
        out_table['inner'][quad_map[str(bin)]] = (pix_total_inner - in_tmp) / pix_total_inner
        out_table['outer'][quad_map[str(bin)]] = (pix_total_outer - out_tmp) / pix_total_outer

    return out_table


if __name__ == "__main__":
    quad_list = [[-2.7, -2.1, 1.42, 7.55, 12.79, 0.287],
                 [-1.6, -1.01, 1.465, -1.387, 3.182, -1.798],
                 [-0.6, 0.0, 1.452, 0.6599, 2.931, -0.053],
                 [0.5, 1.15, 1.4, 7.073, 11.65, 3.29],
                 [1.53, 2.10, 1.34, -5.985, 9.299, -0.4137],
                 [2.56, 3.2, 1.34, 1.633, 3.098, 3.468]]
    mods = [-1, 1, -2, 2, -3, 3, -4, 4]
    mods_pairs = ['-4-3', '-3-2', '-2-1', '-11', '12', '23', '34']
    # in_file_name = './PixelTree_out_rereco_size4_29119.root'
    # in_file_name = './pixelTree_out_design_12419.root'
    # in_file_name = './pixelTree_out_design_allevts_13519.root'
    # in_file_name = './pixelTree_out_realistic_3000evts_1pix_with_outer_23519.root'
    # in_file_name = './pixelTree_out_design_10000evts_1pix_0p0_15719.root'
    in_file_name = './pixelTree_out_realistic_5000evts_1pix_0p01_15719.root'
    # in_file_name = './pixelTree_out_design_3sig_5000_12419.root'
    # in_file_name = './pixelTree_out_realistic_12419.root'
    # in_file_name = './pixelTree_out_realistic_sig3_5000_12419.root'
    gr_template = 'gr_r_phi_{}'
    hist_template = 'h_clust_phi_{}'
    hist_outer_template = 'h_clust_outer_phi_{}'
    params = [-0.4, -1.3, 1.2, 3.]

    # gr_file_name = './out_hists_xy_171218.root'
    # gr_hist_name = 'hclustOccIn_scaled_cleaned_gr'
    # gr_file_name = './out_hists_xy_noscale_171218.root'
    # gr_hist_name = 'hClustOccIn_cleaned_gr'
    # in_file = rt.TFile(gr_file_name, 'r')
    # in_gr = in_file.Get(gr_hist_name)
    # it_ellipse(in_gr)

    # hist_name_1 = 'hClustOccIn_1'
    # # gr_name_1 = 'gr_r_phi_4'
    # # gr_name_1 = 'gr_r_phi_neg_-1'
    # # # hist_name_1 = 'h_clust_phi_-4-3'
    # in_file = rt.TFile(in_file_name, 'r')
    # # gr = in_file.Get(gr_name_1)
    # hist1 = in_file.Get(hist_name_1)
    # hist1_clean = preprocess_hist(hist1)
    # # xy_gr = make_xy_graph(None, gr)
    # xy_gr = make_xy_graph(hist1_clean, None)
    # # out_file = rt.TFile('./out_hists_xy_occ_4219.root', 'recreate')
    # # out_file.cd()
    # # xy_gr.Write()
    # # gr.Write()
    # # hist1.Write()
    # # hist1_clean.Write()
    # # out_file.Close()
    # fit_ellipse(xy_gr)

    in_file = rt.TFile(in_file_name, 'r')
    mods_pairs = ['-4-3', '-3-2', '-2-1', '-11', '12', '23', '34']
    # mods_pairs = ['-4-3', '-3-2', '-2-1', '-11', '12', '23', '34']
    # mods = [1, -2]
    # mods = [-4,-3,-2,-1]
    total_table = OrderedDict()
    # total_table_err = OrderedDict()
    for mod in mods:
        if isinstance(mod, int) and mod < 0:
            hist_name = hist_template.format('neg_'+str(abs(mod)))
            hist_outer_name = hist_outer_template.format('neg_'+str(abs(mod)))
            # hist_name_2 = hist_template.format(str(abs(mod)))
            # hist_name = gr_template.format('neg_'+str(mod))
        else:
            hist_name = hist_template.format(str(mod))
            hist_outer_name = hist_outer_template.format(str(mod))
            # hist_name = gr_template.format(str(mod))
        hist_tmp = in_file.Get(hist_name)
        hist_outer_tmp = in_file.Get(hist_outer_name)
        # hist_tmp_2 = in_file.Get(hist_name_2)
        # hist_tmp.Add(hist_tmp_2)
        hist_tmp_clean = preprocess_hist(hist_tmp)
        # gr_tmp = make_xy_graph(hist_tmp_clean, None)
        # gr_tmp = make_xy_graph(None, hist_tmp)
        print '################'
        print 'mod', mod
        print '################'
        tmp_table = make_occ_table(hist_tmp, hist_outer_tmp)
        total_table[str(mod)] = tmp_table
        # fit_ellipse(gr_tmp, 'design_single_mods', mod)
    print '#####################'
    out_lines = []
    for factor, mod in enumerate(total_table):
        out_lines.append('\nmodule section: {:<35}\n'.format(mod))
        for iline, lad, in enumerate(total_table[mod]):
            if iline == 0:
                out_lines.append('{:35}'.format(' '))
            out_lines.append('{:^35}'.format(lad))
            for quad, value in total_table[mod][lad].items():
                if iline == 0:
                    out_lines[0+1+factor*14] += '{:^35}'.format(quad)
                out_lines[(iline+2)+factor*14] += '{:>15.2f} +/- {:<15.2f}'.format(value, np.sqrt(value))
            out_lines[(iline+2)+factor*14] += '\n'
        out_lines[0+1+factor*14] += '\n'


    # with open('quadrant_occupancies_realistic_3000evts_5pix_with_outer.txt', 'w') as f:
    #     f.writelines(out_lines)
    #     f.close()

    #hist_template = 'h_clust_missing_{}'
    #hist_outer_template = 'h_clust_outer_missing_{}'
    #total_table = OrderedDict()
    #for mod in mods:
    #    if isinstance(mod, int) and mod < 0:
    #        hist_name = hist_template.format('neg_'+str(abs(mod)))
    #        hist_outer_name = hist_outer_template.format('neg_'+str(abs(mod)))
    #        # hist_name_2 = hist_template.format(str(abs(mod)))
    #        # hist_name = gr_template.format('neg_'+str(mod))
    #    else:
    #        hist_name = hist_template.format(str(mod))
    #        hist_outer_name = hist_outer_template.format(str(mod))
    #        # hist_name = gr_template.format(str(mod))
    #    hist_tmp = in_file.Get(hist_name)
    #    hist_outer_tmp = in_file.Get(hist_outer_name)
    #    # hist_tmp_2 = in_file.Get(hist_name_2)
    #    # hist_tmp.Add(hist_tmp_2)
    #    # hist_tmp_clean = preprocess_hist(hist_tmp)
    #    # gr_tmp = make_xy_graph(hist_tmp_clean, None)
    #    # gr_tmp = make_xy_graph(None, hist_tmp)
    #    print '################'
    #    print 'mod', mod
    #    print '################'
    #    tmp_table = make_missing_table(hist_tmp, hist_outer_tmp)
    #    total_table[str(mod)] = tmp_table
    #    # fit_ellipse(gr_tmp, 'realistic_single_mods', mod)
    #print '#####################'
    #out_lines = []
    #for factor, mod in enumerate(total_table):
    #    out_lines.append('\nmodule section: {:<25}\n'.format(mod))
    #    for iline, lad, in enumerate(total_table[mod]):
    #        if iline == 0:
    #            out_lines.append('{:25}'.format(' '))
    #        out_lines.append('{:^25}'.format(lad))
    #        for quad, value in total_table[mod][lad].items():
    #            if iline == 0:
    #                out_lines[0+1+factor*4] += '{:^25}'.format(quad)
    #            out_lines[(iline+2)+factor*4] += '{:>25.8f}'.format(value)
    #        out_lines[(iline+2)+factor*4] += '\n'
    #    out_lines[0+1+factor*4] += '\n'

    # with open('quadrant_design0p0_10000evts_1pix_with_outer.txt', 'w') as f:
    with open('quadrant_realistic0p01_5000evts_1pix_with_outer.txt', 'w') as f:
        f.writelines(out_lines)
        f.close()


    #gr_name_1 = 'gr_r_phi_neg_-1'
    #gr = in_file.Get(gr_name_1)
    #xy_gr = make_xy_graph(None, gr)
    #c0 = rt.TCanvas('c0', '', 800, 600)
    #xy_gr.Draw('A*')
    #xy_gr.SetMarkerStyle(2)
    #c0.Update()

    #e2 = rt.TEllipse(9.7414e-2, -1.01673e-1, 2.73255, , 0., 360., 0.434*(180./rt.TMath.Pi()))
    #e2.SetLineColor(rt.kRed)
    #e2.SetLineWidth(4)
    #e2.SetFillStyle(0)
    # el = rt.TEllipse(1.05043e-01, -9.38342e-2, 2.75181, 2.8033, 0., 360., 0.812466*(180./rt.TMath.Pi()))
    # el.SetLineColor(rt.kBlack)
    # el.SetLineWidth(4)
    # el.SetFillStyle(0)
    # el.Draw("SAME")
    # e2.Draw("SAME")

    #c0.Update()
    #c0.SaveAs("plots_4219/fit_compare_"+gr.GetName()+"_design_r_phi_ellipse_26419.pdf")
    # c0.SaveAs("plots_4219/fit_compare_"+gr.GetName()+"_rereco_r_phi_ellipse_4219.png")
    # c0.SaveAs("plots_4219/fit_compare_"+gr.GetName()+"_rereco_r_phi_ellipse_4219.root")
    # out_file = rt.TFile('./out_hists_xy_-4-3_171218.root', 'recreate')
    # out_file.cd()
    # hist1_clean.Write()
    # hist1_gr.Write()
    # out_file.Close()

    # quad_fits = fit_func_root_quad_mods(in_file_name, gr_template, quad_list, mods)
    # fit_func_root_sine_mods(in_file_name, hist_template, quad_list, quad_fits, mods, params)
