#!/usr/bin/env python

import ROOT as rt
import numpy as np
from pixelMapping_cfi import *
import time
from array import array
from collections import OrderedDict


def read_clusters(input_files):
    tChain = rt.TChain('pixelTree')
    for f in input_files:
        tChain.Add(f)

    of = rt.TFile(input_files[0].split('_')[0] + '_out_61218.root', 'recreate')

    nEvents = tChain.GetEntries()

    r_phi_array = []

    clust_per_mod = OrderedDict()
    mods = [-4, -3, -2, -1, 1, 2, 3, 4]

    for mod in mods:
        clust_per_mod[mod] = {}
        if mod < 0:
            clust_per_mod[mod]['TH2F'] = rt.TH2F('h_clust_r_phi_neg_'+str(abs(mod)), 'cluster r vs phi', 1000, -rt.TMath.Pi(), rt.TMath.Pi(), 1000, 2.0, 3.5)
            clust_per_mod[mod]['TH1F'] = rt.TH1F('h_clust_phi_neg_'+str(abs(mod)), 'cluster phi', 1000, -rt.TMath.Pi(), rt.TMath.Pi())
            clust_per_mod[mod]['TGraph'] = None
            clust_per_mod[mod]['r_phi_array'] = []
        else:
            clust_per_mod[mod]['TH2F'] = rt.TH2F('h_clust_r_phi_'+str(abs(mod)), 'cluster r vs phi', 1000, -rt.TMath.Pi(), rt.TMath.Pi(), 1000, 2.0, 3.5)
            clust_per_mod[mod]['TH1F'] = rt.TH1F('h_clust_phi_'+str(abs(mod)), 'cluster phi', 1000, -rt.TMath.Pi(), rt.TMath.Pi())
            clust_per_mod[mod]['TGraph'] = None
            clust_per_mod[mod]['r_phi_array'] = []

    clustOccOut_init_prehot = rt.TH1F('hClustOccOut_init_prehot', 'Occ vs phi', 100, -3.5, 3.5)
    clustOccIn_init_prehot = rt.TH1F('hClustOccIn_init_prehot', 'Occ vs phi', 100, -3.5, 3.5)
    clustGXY = rt.TH2F('hClustGXGY', 'GY vs GX', 1000, -5, 5, 1000, -5, 5)
    clustIn = rt.TH2F('hClustIn', 'GY vs GX', 100, -5, 5, 100, -5, 5)
    clustOut = rt.TH2F('hClustOut', 'GY vs GX', 100, -5, 5, 100, -5, 5)
    clustOccOut_init = rt.TH1F('hClustOccOut_init', 'Occ vs phi', 100, -3.5, 3.5)
    clustOccIn_init = rt.TH1F('hClustOccIn_init', 'Occ vs phi', 100, -3.5, 3.5)
    clustOccOut = rt.TH1F('hClustOccOut', 'Occ vs phi', 100, -3.5, 3.5)
    clustOccIn = rt.TH1F('hClustOccIn', 'Occ vs phi', 1000, -rt.TMath.Pi(), rt.TMath.Pi())
    clustOccIn_clean = rt.TH1F('hClustOccIn_clean', 'Occ vs phi', 1000, -rt.TMath.Pi(), rt.TMath.Pi())
    clustOccIn_prescale = rt.TH1F('hClustOccIn_prescale', 'Occ vs phi', 1000, -rt.TMath.Pi(), rt.TMath.Pi())
    clustOccIn_plus = rt.TH1F('hClustOccIn_1', 'Occ vs phi', 1000, -rt.TMath.Pi(), rt.TMath.Pi())
    clustOccIn_plus_clean = rt.TH1F('hClustOccIn_1_cleaned', 'Occ vs phi', 1000, -rt.TMath.Pi(), rt.TMath.Pi())

    clustR_init = rt.TH1F('hClustR_init', 'Clusters R', 100, 2.5, 3.0)
    clustR = rt.TH1F('hClustR', 'Clusters R', 100, 2.5, 3.0)
    clustRPhi = rt.TH2F('hClustRPhi', 'Clusters R vs phi', 1000, -3.5, 3.5, 1000, 2.0, 3.5)
    clustRPhi_outer = rt.TH2F('hClustRPhi_outer', 'Clusters R vs phi', 1000, -3.5, 3.5, 1000, 2.0, 4.0)
    clustRZ = rt.TH2F('hClustRZ', 'Clusters R vs Z', 1000, -30., 30., 1000, 0., 20.)
    clustRPhiZ = rt.TH3F('hClustRPhiZ', 'Clusters R vs phi vs z', 500, -3.5, 3.5, 500, 2.0, 4.0, 500, -30, 30)
    clustPerRoc = rt.TH1F('hClustPerRoc', 'Clusters per roc pair along ladders', 10, -4, 4)

    nClusTot = 0.

    # load numpy array that contains list of hot pixels
    hot_pixels = np.load('hot_pixels_5_61218.npy')

    # turn the numpy array into a python list for ease of use in future conditional statement
    # needs to be updated at some point so that it uses a numpy array instead of list
    list_hot_pixels = list(hot_pixels)
    for i, pix in enumerate(list_hot_pixels):
        list_hot_pixels[i] = list(pix)

    for iev in xrange(nEvents):
        if iev % 1 == 0:
            print('Event: ', iev, '/ ', nEvents)
        tChain.GetEntry(iev)
        nClus = tChain.ClN
        nClusTot += nClus
        if iev == 0:
            print(nClus)
        for iCl in xrange(nClus):

            layer = tChain.ClLayer[iCl]

            if layer != 1: continue

            clusGX = tChain.ClGx[iCl]
            clusGY = tChain.ClGy[iCl]
            clusGZ = tChain.ClGz[iCl]
            clR = np.sqrt(clusGX ** 2 + clusGY ** 2)
            try:
                clPhi = np.arctan(clusGX / clusGY)
            except ZeroDivisionError:
                clusGY += 0.00001
                clPhi = np.arctan(clusGX / clusGY)
            clCol = tChain.ClCol[iCl]
            clRow = tChain.ClRow[iCl]
            clLadder = tChain.ClLadder[iCl]
            clMod = tChain.ClModule[iCl]

            if clLadder in [-6, -4, -2, 1, 3, 5]:
                outer = True
                inner = False
                clustOut.Fill(clusGY, clusGX)
            elif clLadder in [-5, -3, -1, 2, 4, 6]:
                inner = True
                outer = False
                clustIn.Fill(clusGY, clusGX)

            if clusGX < 0 < clusGY:
                clPhi += np.pi
            elif clusGX < 0 and clusGY < 0:
                clPhi -= np.pi

            if clR < 2.76 and -1.632 < clPhi < -1.566: continue
            if clR < 2.73 and -0.0003 < clPhi < 0.03: continue

            if 8 < clCol < 46:
                notEdgeCol = True
                clustPerRoc.Fill(-4)
            elif 60 < clCol < 96:
                notEdgeCol = True
                clustPerRoc.Fill(-3)
            elif 112 < clCol < 148:
                notEdgeCol = True
                clustPerRoc.Fill(-2)
            elif 164 < clCol < 200:
                notEdgeCol = True
                clustPerRoc.Fill(-1)
            elif 216 < clCol < 252:
                notEdgeCol = True
                clustPerRoc.Fill(1)
            elif 268 < clCol < 304:
                notEdgeCol = True
                clustPerRoc.Fill(2)
            elif 320 < clCol < 356:
                notEdgeCol = True
                clustPerRoc.Fill(3)
            elif 372 < clCol < 408:
                notEdgeCol = True
                clustPerRoc.Fill(4)
            else:
                notEdgeCol = False

            if 8 < clRow < 72:
                notEdgeRow_8 = True
            elif 88 < clRow < 152:
                notEdgeRow_8 = True
            else:
                notEdgeRow_8 = False

            clustGXY.Fill(clusGY, clusGX)

            if inner:
                clustOccIn_init_prehot.Fill(clPhi)
            elif outer:
                clustOccOut_init_prehot.Fill(clPhi)

            # check to see if pixel exists in the list of hot pixels, the inefficient part of the code
            if [int(clRow), int(clCol), clMod, clLadder] in list_hot_pixels:
                continue

            clustRPhiZ.Fill(clPhi, clR, clusGZ)

            if clLadder == 6:
                clustRZ.Fill(clusGZ, clR)
            if inner:
                clustOccIn_init.Fill(clPhi)
                clustR_init.Fill(clR)

                if notEdgeCol and notEdgeRow_8:
                    clust_per_mod[clMod]['TH2F'].Fill(clPhi, clR)
                    clust_per_mod[clMod]['TH1F'].Fill(clPhi)
                    clust_per_mod[clMod]['r_phi_array'].append((clPhi, clR))

                    clustOccIn.Fill(clPhi)
                    clustOccIn_prescale.Fill(clPhi)
                    clustRPhi.Fill(clPhi, clR)
                    r_phi_array.append((clPhi, clR))

                    clustR.Fill(clR)
                    clustOccIn_plus.Fill(clPhi)

            elif outer:
                clustOccOut_init.Fill(clPhi)
                if notEdgeCol and notEdgeRow_8:
                    clustOccOut.Fill(clPhi)
                    clustRPhi_outer.Fill(clPhi, clR)

    clustOccIn_plus.Scale(1. / nEvents)
    clustOccIn.Scale(1. / nEvents)
    clustOccIn_prescale.Scale(1. / nEvents)
    r_phi_map = []

    for ibin in xrange(clustRPhi.GetNbinsX()):
        r_values = []
        for jbin in xrange(clustRPhi.GetNbinsY()):
            if clustRPhi.GetBinContent(ibin,jbin) > 0.:
                r_values.append(clustRPhi.GetYaxis().GetBinCenter(jbin))
        if not r_values:
            r_mean = 0
        else:
            r_mean = sum(r_values) / float(len(r_values))
        r_phi_map.append(r_mean)

    clustOccIn_plus_scaled = clustOccIn_plus.Clone("hclustOccIn_1_scaled")
    clustOccIn_scaled = clustOccIn.Clone("hclustOccIn_scaled")

    for ibin in xrange(clustOccIn.GetNbinsX()):
        occ = clustOccIn.GetBinContent(ibin)
        occ_plus = clustOccIn_plus.GetBinContent(ibin)
        occ_scaled = occ * (r_phi_map[ibin]**2/2.776**2)
        occ_plus_scaled = occ_plus * (r_phi_map[ibin]**2/2.776**2)

        if occ_scaled > 0.:
            clustOccIn_scaled.SetBinContent(ibin, occ_scaled)
        if occ_plus_scaled > 0.:
            clustOccIn_plus_scaled.SetBinContent(ibin, occ_plus_scaled)

    for i in xrange(clustOccIn_plus.GetNbinsX()):
        occBin = clustOccIn_scaled.GetBinContent(i)
        phiBin = clustOccIn_scaled.GetBinCenter(i)
        plusBin = clustOccIn_plus_scaled.GetBinContent(i)

        for phi1, phi2, occ in bin_edge_map:
            if phi1 < phiBin < phi2 and occBin > occ:
                clustOccIn_clean.SetBinContent(i, occBin)
                clustOccIn_plus_clean.SetBinContent(i, plusBin)
    r_vals = array('d')
    phi_vals = array('d')
    for phi, r in r_phi_array:
        r_vals.append(r)
        phi_vals.append(phi)
    n_vals = len(r_phi_array)
    g_r_phi = rt.TGraph(n_vals, phi_vals, r_vals)

    for mod in clust_per_mod.keys():
        clust_per_mod[mod]['TH2F'].Scale(1. / nEvents)
        clust_per_mod[mod]['TH1F'].Scale(1. / nEvents)
        r_vals = array('d')
        phi_vals = array('d')
        for phi, r in clust_per_mod[mod]['r_phi_array']:
            r_vals.append(r)
            phi_vals.append(phi)
        n_vals = len(clust_per_mod[mod]['r_phi_array'])
        clust_per_mod[mod]['TGraph'] = rt.TGraph(n_vals, phi_vals, r_vals)
        if mod < 0:
            clust_per_mod[mod]['TGraph'].SetName('gr_r_phi_neg_'+str(mod))
        else:
            clust_per_mod[mod]['TGraph'].SetName('gr_r_phi_'+str(mod))
        clust_per_mod[mod]['TGraph'].Write()
    of.Write()
    g_r_phi.Write()
    of.Close()


if __name__ == "__main__":
    from inputFiles_cfi import *
    t_start = time.time()
    read_clusters(input_files_2)
    t_stop = time.time()
    print t_stop - t_start

