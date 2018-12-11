#!/usr/bin/env python

import ROOT as rt
import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict

rt.gROOT.SetBatch()
rt.gStyle.SetOptStat(111)
rt.gStyle.SetOptFit(1111)


def quad_func(phi, par):
    quad_list = [[-2.7, -2.1, 1.42, 7.55, 12.79, 0.287],
                 [-1.6, -1.01, 1.465, -1.387, 3.182, -1.798],
                 [-0.6, 0.0, 1.452, 0.6599, 2.931, -0.053],
                 [0.5, 1.15, 1.4, 7.073, 11.65, 3.29],
                 [1.53, 2.10, 1.34, -5.985, 9.299, -0.4137],
                 [2.56, 3.2, 1.34, 1.633, 3.098, 3.468]]
    r = 2.7
    for quad in quad_list:
        phi_min = quad[0]
        phi_max = quad[1]
        if phi_min < phi < phi_max:
            r = quad[2]*(phi[0] - quad[5])**2 + quad[3] * (phi[0] - quad[5]) + quad[4]
        else: continue
    return r * (par[0] * np.sin(par[1] * phi[0] - par[2])) + par[3]


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


def fit_func_root_quad_mods(inFile, histName_template, quadList):
    # f = rt.TFile(inFile, 'r')

   # h = f.Get(histName)
   # h_all = f.Get(histName)

    mods = [-4, -3, -2, -1, 1, 2, 3, 4]

    fit_results_mods = OrderedDict()

    c1 = rt.TCanvas('c0', '', 800, 600)
    for mod in mods:
        f = rt.TFile(inFile, 'r')
        if mod < 0:
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
            print f
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

            fit_results_mods[mod]['results'][i].Print('V')

            fit_results_mods[mod][histName].Draw('A*')

            fit_results_mods[mod]['funcs'][i].Draw('same')

            c1.SaveAs('d_' + histName + '_' + str(i) + '.png')
        c1.SaveAs('d_' + histName + '_' + str(i) + '.root')
        rt.gStyle.SetOptFit(0)
        fit_results_mods[mod][histName + '_all'].Draw('A*')
        for i, f in enumerate(fit_results_mods[mod]['funcs']):
            f.Draw('same')
        c1.SaveAs('d_' + histName + '_quad_all' + '.png')
        c1.SaveAs('d_' + histName + '_quad_all' + '.root')
    return fit_results_mods



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


if __name__ == "__main__":
    quad_list = [[-2.7, -2.1, 1.42, 7.55, 12.79, 0.287],
                 [-1.6, -1.01, 1.465, -1.387, 3.182, -1.798],
                 [-0.6, 0.0, 1.452, 0.6599, 2.931, -0.053],
                 [0.5, 1.15, 1.4, 7.073, 11.65, 3.29],
                 [1.53, 2.10, 1.34, -5.985, 9.299, -0.4137],
                 [2.56, 3.2, 1.34, 1.633, 3.098, 3.468]]
    mods = [-4, -3, -2, -1, 1, 2, 3, 4]



    #fit_func_root_quad('./files/PixelTree_out_test.root', 'Graph', quad_list)
    fit_func_root_quad_mods('./files/PixelTree_out_31218.root', 'gr_r_phi_{}', quad_list)
    #fit_func_root_2d('./files/PixelTree_out_test.root', 'hClustOccIn_clean')
    #    fitFuncRoot('clus_ana_317650_10Jul18.root', 'hclumultxPhi1', 0.5, 12, -0.5, 2.1)
    # fitFuncRoot('clus_ana_317650_10Jul18.root', 'hpixcharPhi1', -0.3, 13, -5, 14)
    #fit_func_root_sine('./files/PixelTree_out_test.root', 'hClustOccIn_clean', -0.4, -1.3, 33., 3.)
    plot_circle(-0.229685, -1.28837, 33.1579, 1.53185, 'circle_full_new')
    #plot_circle_func(-0.085, -1.288, 33.16, 0.567354, 'circle_full_test')
    plot_circle_func(-0.085, -1.288, 33.16, 1.53186, 'circle_full_test_v2')
