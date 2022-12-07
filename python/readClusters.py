#!/usr/bin/env python

# readClusters.py

import sys
import os
import time
import ROOT as rt
import numpy as np
from pixelMapping_cfi import *
from tools import getChain, get_file_list, get_eos_file_list
from array import array
from collections import OrderedDict

def make_cluster_map(input_files_, num_files, doClusterSizeCut, doNumberClustersCut):
    """
    function to make a cluster occupancy map of layer 1 pixels, with 1:1 pixel, index correspondence
    :param input_files: list of PixelTree ntuples to be processed
    :return: numpy file with cluster occupancy
    """
    debug = False
    chain = getChain(input_files_, num_files)

    # creating empty array to be used as map
    occ_array = np.full((1920, 3328), None)

    n_events = chain.GetEntries()
    max_events = -1

    print("Number of events: {0}".format(n_events))

    info_dict = {}
    phi_dict  = {}

    for iev, event in enumerate(chain):
        if max_events > 0 and iev > max_events:
            break
        if iev % 10 == 0:
            print('Event', iev, ' / ', n_events)
            sys.stdout.flush()
        
        n_cl = event.ClN
        # cut on number of clusters per event
        if doNumberClustersCut:
            if n_cl < 2000:
                continue
        #print("event: {0}, number of clusters: {1}".format(iev, n_cl))

        for icl in range(n_cl):

            detid_cl = event.ClDetId[icl]

            # skipping detids outside of layer 1
            if not(303000000 < detid_cl < 304000000):
                continue  # layer 1

            clus_size   = chain.ClSize[icl]
            clus_size_x = chain.ClSizeX[icl]
            clus_size_y = chain.ClSizeY[icl]
            
            # cut on cluster size
            if doClusterSizeCut:
                if clus_size < 2:
                    continue

            # cut on number of tracks per cluster
            
            # select on track clusters
            #if chain.ClTkN[icl] < 1: continue 
            
            # select off track clusters
            #if chain.ClTkN[icl] >= 1: continue 

            row_cl    = int(event.ClRow[icl])
            col_cl    = int(event.ClCol[icl])
            ladder_cl = int(event.ClLadder[icl])

            # remove if not an inner ladder
            if ladder_cl not in [-5, -3, -1, 2, 4, 6]: continue

            charge = event.ClCharge[icl]
            #if charge > 1000: continue
            gx_cl = event.ClGx[icl]
            gy_cl = event.ClGy[icl]
            gr_cl = np.sqrt(gx_cl**2 + gy_cl**2) 
            gz_cl = event.ClGz[icl]
            try:
                phi_cl = np.arctan2(gy_cl, gx_cl)
            except ZeroDivisionError:
                print("ERROR: ZeroDivisionError")
                gy_cl += 0.00001
                phi_cl = np.arctan2(gy_cl, gx_cl)

            # converting detids to a ladder index and module index
            ladder_index = int(ladder_map(detid_cl))
            module_index = int(module_map(detid_cl, ladder_index))

            array_row = row_cl + (160*ladder_index)
            array_col = col_cl + (416*module_index)
 
            if not occ_array[array_row][array_col]:
                occ_array[array_row][array_col] = [0, gr_cl, phi_cl, gz_cl]
 
            edge_row, edge_col = edge_pix_map(row_cl, col_cl)
            if edge_row or edge_col: continue 

            occ_array[row_cl+(160*ladder_index)][col_cl+(416*module_index)][0] += 1

            if debug:
                info = "Event {0}: ladder = {1}, module = {2}, cluster {3}, x = {4:.3f}, y = {5:.3f}, z = {6:.3f}, phi = {7:.3f}".format(iev, ladder_index, module_index, icl, gx_cl, gy_cl, gz_cl, phi_cl)
                
                # record info
                if ladder_index in info_dict:
                    info_dict[ladder_index].append(info)
                else:
                    info_dict[ladder_index] = [info]
                
                # record phi
                if ladder_index in phi_dict:
                    phi_dict[ladder_index].append(phi_cl)
                else:
                    phi_dict[ladder_index] = [phi_cl]

    if debug:
        # print(phi)
        for ladder_index in phi_dict:
            phi_list = phi_dict[ladder_index]
            phi_min  = np.min(phi_list)
            phi_max  = np.max(phi_list)
            phi_mean = np.mean(phi_list)
            phi_std  = np.std(phi_list)
            print("ladder = {0}: phi range = [{1:.3f}, {2:.3f}], phi_mean = {3:.3f} +/- {4:.3f}".format(ladder_index, phi_min, phi_max, phi_mean, phi_std))
            # split ladder 3 into phi < 0 and phi > 0, as ladder three crosses phi = -pi = pi
            if ladder_index == 3:
                phi_neg_list = []
                phi_pos_list = []
                for phi in phi_list:
                    if phi < 0.0:
                        phi_neg_list.append(phi)
                    else:
                        phi_pos_list.append(phi)
                # phi < 0
                phi_neg_min  = np.min(phi_neg_list)
                phi_neg_max  = np.max(phi_neg_list)
                phi_neg_mean = np.mean(phi_neg_list)
                phi_neg_std  = np.std(phi_neg_list)
                # phi >= 0
                phi_pos_min  = np.min(phi_pos_list)
                phi_pos_max  = np.max(phi_pos_list)
                phi_pos_mean = np.mean(phi_pos_list)
                phi_pos_std  = np.std(phi_pos_list)
                print("ladder = {0}: require phi <  0: phi range = [{1:.3f}, {2:.3f}], phi_mean = {3:.3f} +/- {4:.3f}".format(ladder_index, phi_neg_min, phi_neg_max, phi_neg_mean, phi_neg_std))
                print("ladder = {0}: require phi >= 0: phi range = [{1:.3f}, {2:.3f}], phi_mean = {3:.3f} +/- {4:.3f}".format(ladder_index, phi_pos_min, phi_pos_max, phi_pos_mean, phi_pos_std))

    return occ_array

def read_clusters(input_files, f_name_):
    tChain = rt.TChain('pixelTree')
    for f in input_files:
        tChain.Add(f)

    of = rt.TFile('./output_files/'+f_name_+'.root', 'recreate')

    nEvents = tChain.GetEntries()

    r_phi_array = []

    clust_per_mod = OrderedDict()
    clust_per_mod_outer = OrderedDict()
    clust_per_mod_pair = OrderedDict()
    mods = [-4, -3, -2, -1, 1, 2, 3, 4]
    mods_pairs = ['-4-3', '-3-2', '-2-1', '-11', '12', '23', '34', '-44', '-33', '-22']

    for mod in mods:
        clust_per_mod[mod] = {}
        clust_per_mod_outer[mod] = {}
        if mod < 0:
            clust_per_mod[mod]['TH2F'] = rt.TH2F('h_clust_r_phi_neg_'+str(abs(mod)), 'cluster r vs phi', 1000, -rt.TMath.Pi(), rt.TMath.Pi(), 1000, 2.0, 3.5)
            clust_per_mod[mod]['TH1F'] = rt.TH1F('h_clust_phi_neg_'+str(abs(mod)), 'cluster phi', 5000, -rt.TMath.Pi(), rt.TMath.Pi())
            clust_per_mod[mod]['quad'] = rt.TH1F('h_clust_missing_neg_'+str(abs(mod)), 'cluster quadrant', 4, 1, 5)
            clust_per_mod[mod]['TGraph'] = None
            clust_per_mod[mod]['r_phi_array'] = []
            clust_per_mod_outer[mod]['TH2F'] = rt.TH2F('h_clust_outer_r_phi_neg_'+str(abs(mod)), 'cluster r vs phi', 1000, -rt.TMath.Pi(), rt.TMath.Pi(), 1000, 2.0, 3.5)
            clust_per_mod_outer[mod]['TH1F'] = rt.TH1F('h_clust_outer_phi_neg_'+str(abs(mod)), 'cluster phi', 5000, -rt.TMath.Pi(), rt.TMath.Pi())
            clust_per_mod_outer[mod]['quad'] = rt.TH1F('h_clust_outer_missing_neg_'+str(abs(mod)), 'cluster quadrant', 4, 1, 5)
            clust_per_mod_outer[mod]['TGraph'] = None
            clust_per_mod_outer[mod]['r_phi_array'] = []
        else:
            clust_per_mod[mod]['TH2F'] = rt.TH2F('h_clust_r_phi_'+str(abs(mod)), 'cluster r vs phi', 10000, -rt.TMath.Pi(), rt.TMath.Pi(), 1000, 2.0, 3.5)
            clust_per_mod[mod]['TH1F'] = rt.TH1F('h_clust_phi_'+str(abs(mod)), 'cluster phi', 5000, -rt.TMath.Pi(), rt.TMath.Pi())
            clust_per_mod[mod]['quad'] = rt.TH1F('h_clust_missing_'+str(abs(mod)), 'cluster quadrant', 4, 1, 5)
            clust_per_mod[mod]['TGraph'] = None
            clust_per_mod[mod]['r_phi_array'] = []
            clust_per_mod_outer[mod]['TH2F'] = rt.TH2F('h_clust_outer_r_phi_'+str(abs(mod)), 'cluster r vs phi', 1000, -rt.TMath.Pi(), rt.TMath.Pi(), 1000, 2.0, 3.5)
            clust_per_mod_outer[mod]['TH1F'] = rt.TH1F('h_clust_outer_phi_'+str(abs(mod)), 'cluster phi', 5000, -rt.TMath.Pi(), rt.TMath.Pi())
            clust_per_mod_outer[mod]['quad'] = rt.TH1F('h_clust_outer_missing_'+str(abs(mod)), 'cluster quadrant', 4, 1, 5)
            clust_per_mod_outer[mod]['TGraph'] = None
            clust_per_mod_outer[mod]['r_phi_array'] = []

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
    beamspot_x = rt.TH1F('hBsX', 'beam spot location, x direction', 50, -1, 1)
    beamspot_y = rt.TH1F('hBsY', 'beam spot location, y direction', 50, -1, 1)
    beamspot_width_x = rt.TH1F('hwidthBsX', 'beam spot width, x direction', 50, -1, 1)
    beamspot_width_y = rt.TH1F('hwidthBsY', 'beam spot width, y direction', 50, -1, 1)

    nClusTot = 0.

    # load numpy array that contains list of hot pixels
    # hot_pixels = np.load('hot_pixels_rereco_5_29119.npy')
    # hot_pixels = np.load('hot_pixels_design_3_10419.npy')
    # hot_pixels = np.load('hot_pixels_realistic_3_21519.npy')
    # hot_rocs = np.load('hot_rocs_realistic_3_10419.npy')
    # hot_pixels = np.load('hot_pixels_express_5_29119.npy')

    # turn the numpy array into a python list for ease of use in future conditional statement
    # needs to be updated at some point so that it uses a numpy array instead of list
    # list_hot_pixels = list(hot_pixels)
    # for i, pix in enumerate(list_hot_pixels):
    #    list_hot_pixels[i] = list(pix)
    # hot_pixels = [list(pix) for pix in hot_pixels]
    # hot_rocs = [list(roc) for roc in hot_rocs]
    for iev in range(nEvents):
        if iev % 10 == 0:
            print('Event: ', iev, '/ ', nEvents)
            sys.stdout.flush()
        tChain.GetEntry(iev)
        nClus = tChain.ClN
        nClusTot += nClus
        if iev == 0:
            print(nClus)
            sys.stdout.flush()
        beamspot_x.Fill(tChain.BsX)
        beamspot_y.Fill(tChain.BsY)
        beamspot_width_x.Fill(tChain.Bs_widthX)
        beamspot_width_y.Fill(tChain.Bs_widthY)
        for iCl in range(nClus):

            layer = tChain.ClLayer[iCl]

            if layer != 1: continue
            clus_size_x = tChain.ClSizeX[iCl]
            clus_size_y = tChain.ClSizeY[iCl]
            clus_charge = tChain.ClCharge[iCl]
            if clus_charge < 0: continue
            #if clus_charge < 50: continue
            # if clus_size_x != 2: continue
            # if clus_size_y != 2: continue
            if clus_size_y + clus_size_x < 3: continue
            # if tChain.ClTkN[iCl] > 0: continue 

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
            # if [int(clRow), int(clCol), clMod, clLadder] in list_hot_pixels:
            #     continue

            clustRPhiZ.Fill(clPhi, clR, clusGZ)

            if clLadder == 6:
                clustRZ.Fill(clusGZ, clR)
            if inner:
                clustOccIn_init.Fill(clPhi)
                clustR_init.Fill(clR)

                if notEdgeCol and notEdgeRow_8:
                    # check to see if pixel exists in the list of hot pixels, the inefficient part of the code
                    # roc_row, roc_col = roc_map(int(clRow), int(clCol))
                    # try:
                    #     hot_rocs.index([int(roc_row), int(roc_col), clMod, clLadder])
                    #     continue
                    # except ValueError:
                    #     pass
                    # try:
                    #     hot_pixels.index([int(clRow), int(clCol), clMod, clLadder])
                    #     continue
                    # except ValueError:
                    #     pass
                    # for mod in mods_pairs:
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
                    roc_row, roc_col = roc_map(int(clRow), int(clCol))
                    # try:
                    #     hot_rocs.index([int(roc_row), int(roc_col), clMod, clLadder])
                    #     continue
                    # except ValueError:
                    #     pass
                    # try:
                    #     hot_pixels.index([int(clRow), int(clCol), clMod, clLadder])
                    #     continue
                    # except ValueError:
                    #     pass
                    clustOccOut.Fill(clPhi)
                    clustRPhi_outer.Fill(clPhi, clR)
                    clust_per_mod_outer[clMod]['TH2F'].Fill(clPhi, clR)
                    clust_per_mod_outer[clMod]['TH1F'].Fill(clPhi)
                    clust_per_mod_outer[clMod]['r_phi_array'].append((clPhi, clR))

    # clustOccIn_plus.Scale(1. / nEvents)
    # clustOccIn.Scale(1. / nEvents)
    # clustOccIn_prescale.Scale(1. / nEvents)
    # r_phi_map = []

    # for ibin in range(clustRPhi.GetNbinsX()):
    #     r_values = []
    #     for jbin in range(clustRPhi.GetNbinsY()):
    #         if clustRPhi.GetBinContent(ibin,jbin) > 0.:
    #             r_values.append(clustRPhi.GetYaxis().GetBinCenter(jbin))
    #     if not r_values:
    #         r_mean = 0
    #     else:
    #         r_mean = sum(r_values) / float(len(r_values))
    #     r_phi_map.append(r_mean)

    # clustOccIn_plus_scaled = clustOccIn_plus.Clone("hclustOccIn_1_scaled")
    # clustOccIn_scaled = clustOccIn.Clone("hclustOccIn_scaled")

    # for ibin in range(clustOccIn.GetNbinsX()):
    #     occ = clustOccIn.GetBinContent(ibin)
    #     occ_plus = clustOccIn_plus.GetBinContent(ibin)
    #     occ_scaled = occ * (r_phi_map[ibin]**2/2.776**2)
    #     occ_plus_scaled = occ_plus * (r_phi_map[ibin]**2/2.776**2)

    #     if occ_scaled > 0.:
    #         clustOccIn_scaled.SetBinContent(ibin, occ_scaled)
    #     if occ_plus_scaled > 0.:
    #         clustOccIn_plus_scaled.SetBinContent(ibin, occ_plus_scaled)

    # for i in range(clustOccIn_plus.GetNbinsX()):
    #     occBin = clustOccIn_scaled.GetBinContent(i)
    #     phiBin = clustOccIn_scaled.GetBinCenter(i)
    #     plusBin = clustOccIn_plus_scaled.GetBinContent(i)

    #     for phi1, phi2, occ in bin_edge_map:
    #         if phi1 < phiBin < phi2 and occBin > occ:
    #             clustOccIn_clean.SetBinContent(i, occBin)
    #             clustOccIn_plus_clean.SetBinContent(i, plusBin)
    # r_vals = array('d')
    # phi_vals = array('d')
    # for phi, r in r_phi_array:
    #     r_vals.append(r)
    #     phi_vals.append(phi)
    # n_vals = len(r_phi_array)
    # g_r_phi = rt.TGraph(n_vals, phi_vals, r_vals)

    # for mod in clust_per_mod.keys():
    #     # clust_per_mod[mod]['TH2F'].Scale(1. / nEvents)
    #     # clust_per_mod[mod]['TH1F'].Scale(1. / nEvents)
    #     r_vals = array('d')
    #     phi_vals = array('d')
    #     for phi, r in clust_per_mod[mod]['r_phi_array']:
    #         r_vals.append(r)
    #         phi_vals.append(phi)
    #     n_vals = len(clust_per_mod[mod]['r_phi_array'])
    #     clust_per_mod[mod]['TGraph'] = rt.TGraph(n_vals, phi_vals, r_vals)
    #     if mod < 0:
    #         clust_per_mod[mod]['TGraph'].SetName('gr_r_phi_neg_'+str(mod))
    #     else:
    #         clust_per_mod[mod]['TGraph'].SetName('gr_r_phi_'+str(mod))
    #     clust_per_mod[mod]['TGraph'].Write()
    #     for phi, r in clust_per_mod_outer[mod]['r_phi_array']:
    #         r_vals.append(r)
    #         phi_vals.append(phi)
    #     n_vals = len(clust_per_mod_outer[mod]['r_phi_array'])
    #     clust_per_mod_outer[mod]['TGraph'] = rt.TGraph(n_vals, phi_vals, r_vals)
    #     if mod < 0:
    #         clust_per_mod_outer[mod]['TGraph'].SetName('gr_r_phi_neg_'+str(mod))
    #     else:
    #         clust_per_mod_outer[mod]['TGraph'].SetName('gr_r_phi_'+str(mod))
    #     clust_per_mod_outer[mod]['TGraph'].Write()

    # for pix in hot_pixels:
    #     row = pix[0]
    #     col = pix[1]

    #     edge_row, edge_col = edge_pix_map(row, col)

    #     if edge_row or edge_col:
    #         continue
    #     mod = pix[2]
    #     lad = pix[3]

    #     if lad == -6 or lad == 1:
    #         continue

    #     quad = ladder_quad_map(lad)

    #     inner_ladders = np.array([-5, -3, -1, 2, 4, 6])
    #     outer_ladders = np.array([-6, -4, -2, 1, 3, 5])
    #     if np.any(inner_ladders[inner_ladders==lad]):
    #         clust_per_mod[mod]['quad'].Fill(quad)
    #     elif np.any(outer_ladders[outer_ladders==lad]):
    #         clust_per_mod_outer[mod]['quad'].Fill(quad)

    # for roc in hot_rocs:
    #     mod = roc[2]
    #     lad = roc[3]

    #     if lad == -6 or lad == 1:
    #         continue

    #     quad = ladder_quad_map(str(lad))

    #     inner_ladders = np.array([-5, -3, -1, 2, 4, 6])
    #     outer_ladders = np.array([-6, -4, -2, 1, 3, 5])
    #     n_pixels = 2304
    #     if np.any(inner_ladders[inner_ladders==lad]):
    #         for ipix in range(n_pixels):
    #             clust_per_mod[mod]['quad'].Fill(quad)
    #     elif np.any(outer_ladders[outer_ladders==lad]):
    #         for ipix in range(n_pixels):
    #             clust_per_mod_outer[mod]['quad'].Fill(quad)

    of.Write()
    #g_r_phi.Write()
    of.Close()


def run(directory, output_file, message, num_files, doClusterSizeCut, doNumberClustersCut):
    print(message)
    file_list   = get_file_list(directory)
    occ_map     = make_cluster_map(file_list, num_files, doClusterSizeCut, doNumberClustersCut)
    np.save(output_file, occ_map)

def runEOS(directory, output_file, message, num_files, doClusterSizeCut, doNumberClustersCut):
    print(message)
    file_list   = get_eos_file_list(directory)
    occ_map     = make_cluster_map(file_list, num_files, doClusterSizeCut, doNumberClustersCut)
    np.save(output_file, occ_map)

def runSingleMuon2017B():
    # Single Muon 2017B v1
    # 20 files: 125008 events
    #directory   = '/mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/user/caleb/PixelTrees/SingleMuon/crab_PixelTree_SingleMuon_2017B_RAW_v1/210810_160009/0000' 
    
    # Single Muon 2017B v2
    # 20 files: 124586 events
    # 23 files: 142531 events
    directory   = '/mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/user/caleb/PixelTrees/SingleMuon/crab_PixelTree_SingleMuon_2017B_RAW_v2/210909_000356/0000' 
    
    #output_file = 'SingleMuon_2017B_MoreEvents_AllClusters.npy'
    #output_file = 'SingleMuon_2017B_MoreEvents_ClusterSize2_AllClusters.npy'
    
    output_file = 'SingleMuon_2017B_Legacy_MoreEvents_ClusterSize2_NumberClusters2000_AllClusters.npy'
    
    message     = 'Running over SingleMuon 2017B PixelTrees.'
    num_files   = 23
    doClusterSizeCut = True
    doNumberClustersCut = True
    
    run(directory, output_file, message, num_files, doClusterSizeCut, doNumberClustersCut)

def runZeroBias2017B():
    
    # Zero Bias 2017B v1
    #   8 files:   68332 events
    #  15 files:  136644 events
    # 150 files: 1277957 events 
    #directory   = '/mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/user/caleb/PixelTrees/ZeroBias/crab_PixelTree_ZeroBias_2017B_RAW_v2/210824_165140/0000'
    
    # Zero Bias 2017B v2
    # 10 files:  94184 events
    # 15 files: 141746 events
    # 20 files: 183883 events
    directory   = '/mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/user/caleb/PixelTrees/ZeroBias/crab_PixelTree_ZeroBias_2017B_RAW_v4/210909_195432/0000'
    
    #output_file = 'ZeroBias_2017B_MoreEvents_AllClusters.npy' 
    #output_file = 'ZeroBias_2017B_MoreEvents_ClusterSize2_AllClusters.npy' 
    #output_file = 'ZeroBias_2017B_MoreEvents2_AllClusters.npy' 
    #output_file = 'ZeroBias_2017B_MoreEvents2_ClusterSize2_AllClusters.npy' 
    
    #output_file = 'ZeroBias_2017B_AllClusters.npy' 
    #output_file = 'ZeroBias_2017B_ClusterSize2_AllClusters.npy' 
    #output_file = 'ZeroBias_2017B_ClusterSize2_NumberClusters2000_AllClusters.npy' 
    
    output_file = 'ZeroBias_2017B_Legacy_MoreEvents_ClusterSize2_NumberClusters2000_AllClusters.npy' 
    
    message     = 'Running over ZeroBias 2017B PixelTrees.'
    num_files   = 15
    doClusterSizeCut = True
    doNumberClustersCut = True
    
    run(directory, output_file, message, num_files, doClusterSizeCut, doNumberClustersCut)

def runSingleMuon2018C():
    
    # Single Muon
    # 7041 events in first file, about 2.6 hours run time
    #  1 file:    7041 events
    # 12 files:  70667 events
    # 20 files: 115126 events
    directory   = '/mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/user/caleb/PixelTrees/SingleMuon/crab_PixelTree_SingleMuon_2018C_RAW_Run319337_v1/210403_235502/0000'
    output_file = 'SingleMuon_MoreEvents_AllClusters.npy'
    #output_file = 'SingleMuon_MoreEvents_ClusterSize2_AllClusters.npy'
    message     = 'Running over SingleMuon PixelTrees.'
    num_files   = 20
    doClusterSizeCut = False
    doNumberClustersCut = False

    run(directory, output_file, message, num_files, doClusterSizeCut, doNumberClustersCut)

def runZeroBias2018C():
    
    # Zero Bias
    # 1354 events in first file, about 0.5 hours run time
    #   1 file:    1354 events
    #  2 files:   69598 events
    #  3 files:  137809 events
    # 20 files: 1322727 events
    directory   = '/mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/user/caleb/PixelTrees/ZeroBias/crab_PixelTree_ZeroBias_2018C_RAW_AllRuns_v1/210405_171418/0000'
    #output_file = 'ZeroBias_MoreEvents2_AllClusters.npy' 
    output_file = 'ZeroBias_MoreEvents2_ClusterSize2_AllClusters.npy' 
    message     = 'Running over ZeroBias PixelTrees.'
    num_files   = 20
    doClusterSizeCut = True
    doNumberClustersCut = False
    
    run(directory, output_file, message, num_files, doClusterSizeCut, doNumberClustersCut)

def runMinBias2017B():
    
    # Min Bias
    directory   = '/mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/user/caleb/PixelTrees/MinimumBias/crab_PixelTree_MinBias_2017B_RAW_v3/210715_152921/0000'
    output_file = 'MinBias_2017B_MoreEvents_AllClusters.npy' 
    #output_file = 'MinBias_2017B_OnTrack.npy' 
    #output_file = 'MinBias_2017B_OffTrack.npy' 
    message     = 'Running over MinBias PixelTrees.'
    num_files   = 30
    doClusterSizeCut = False
    doNumberClustersCut = False
    
    run(directory, output_file, message, num_files, doClusterSizeCut, doNumberClustersCut)

def runMinBias2018C():
    
    # Min Bias
    #   1 file:   41084 events
    #  2 files:   82206 events
    #  3 files:  123718 events
    # 30 files: 1227505 events
    directory   = '/mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/user/caleb/PixelTrees/MinimumBias2/crab_PixelTree_MinBias_2018C_RAW_v2/210707_165008/0000'
    #output_file = 'MinBias_2018C_MoreEvents2_AllClusters.npy' 
    output_file = 'MinBias_2018C_MoreEvents2_ClusterSize2_AllClusters.npy' 
    #output_file = 'MinBias_2018C_AllClusters.npy' 
    #output_file = 'MinBias_2018C_OnTrack.npy' 
    #output_file = 'MinBias_2018C_OffTrack.npy' 
    message     = 'Running over MinBias PixelTrees.'
    num_files   = 30
    doClusterSizeCut = True
    doNumberClustersCut = False
    
    run(directory, output_file, message, num_files, doClusterSizeCut, doNumberClustersCut)

def runExpressData2021():
    directory   = '/store/user/lpcsusylep/PixelTrees/ExpressPhysics/crab_PixelTree_Express_2021_Run346512_v1/211110_190014/0000'
    #output_file = 'ExpressPhysics_2021_Run346512_v1.npy'
    output_file = 'ExpressPhysics_2021_Run346512_v2.npy'
    message     = 'Running over Express Data from 2021.'
    #num_files   = 20
    num_files   = -1
    doClusterSizeCut    = True
    doNumberClustersCut = False
    runEOS(directory, output_file, message, num_files, doClusterSizeCut, doNumberClustersCut)

def runZeroBias2022F():
    directory   = '/store/user/lpcsusylep/PixelTrees/ZeroBias/crab_PixelTree_ZeroBias_2022F_RAW_v1_Run362154_v2/221206_141708/0000'
    #output_file = 'data/ZeroBias_2022F_nFiles1_NoCuts.npy'
    output_file = 'data/ZeroBias_2022F_nFiles1_ClustSize2_nClust2000.npy'
    message     = 'Running over Zero Bias 2022F pixel trees.'
    num_files   = 1
    doClusterSizeCut    = True
    doNumberClustersCut = True
    runEOS(directory, output_file, message, num_files, doClusterSizeCut, doNumberClustersCut)

def runTTBarPU():
    
    # TTBar with pileup
    directory   = '/mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/user/eschmitz/PixelTrees/RelValTTbar_13TeV/crab_RelValTTbar_13TeVRAW_eschmitzcrab_design_0p2/190930_182712/0000'
    #output_file = 'TTBar_pileup_0p2_ClusterSize2_AllClusters.npy' 
    output_file = 'TTBar_pileup_0p2_AllClusters.npy' 
    #output_file = 'TTBar_pileup_0p2_OnTrack.npy' 
    #output_file = 'TTBar_pileup_0p2_OffTrack.npy' 
    message     = 'Running over TTBar PU PixelTrees.'
    num_files   = 99
    doClusterSizeCut = False
    doNumberClustersCut = False
    
    run(directory, output_file, message, num_files, doClusterSizeCut, doNumberClustersCut)

def runTTBar():
    
    # TTBar without pileup
    directory   = '/mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/user/eschmitz/PixelTrees/RelValTTbar_13TeV/crab_RelValTTbar_13TeVdesign_0p2_GEN_SIM/190819_222136/0000'
    #output_file = 'TTBar_0p2_ClusterSize2_AllClusters.npy' 
    output_file = 'TTBar_0p2_AllClusters.npy' 
    #output_file = 'TTBar_0p2_OnTrack.npy' 
    #output_file = 'TTBar_0p2_OffTrack.npy' 
    message     = 'Running over TTBar PixelTrees.'
    num_files   = 98
    doClusterSizeCut = False
    doNumberClustersCut = False
    
    run(directory, output_file, message, num_files, doClusterSizeCut, doNumberClustersCut)
    
def main():
    # main
    
    t_start = time.time()
    
    #runZeroBias2017B()
    #runSingleMuon2017B()
    #runMinBias2017B()
    
    #runSingleMuon2018C()
    #runMinBias2018C()
    #runZeroBias2018C()
    
    #runTTBarPU()
    #runTTBar()

    #runExpressData2021()

    runZeroBias2022F()
    
    t_stop = time.time()
    print("run time (sec): {0}".format(t_stop - t_start))

if __name__ == "__main__":
    main()

