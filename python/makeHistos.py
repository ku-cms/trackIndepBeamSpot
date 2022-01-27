# makeHistos.py

import ROOT as rt
import root_numpy as rnp
import numpy as np
import tools
import string
import os
import csv
import scipy.optimize as so
import scipy.stats as ss
import pandas as pd

alpha_low = string.ascii_lowercase

# get ladder for given phi based on phi binning 
def getLadder(phi):
    if phi > np.pi or phi < -np.pi:
        print("ERROR: phi = {0} is outside of [-pi, pi]".format(phi))
        return -999
    phi_bin_edges = [-np.pi, -2.30, -1.75, -1.25, -0.75, -0.25, 0.25, 0.75, 1.25, 1.75, 2.30, 2.90, np.pi]
    for i in range(len(phi_bin_edges)):
        if np.isnan(phi):
            continue
        if phi >= phi_bin_edges[i] and phi < phi_bin_edges[i+1]:
            return i
    #print("No valid ladder for phi = {0}".format(phi))
    return -999

# TODO: Fix cases of repeats in ladder nums (e.g. two 11s)
# TODO: Check cases where missing ladders are on left and right... 
#       Assigned ladders may need to shift left or right

# get ladder indices given phi values
def getLadderNums(phi_vals):
    ladder_nums = []
    ladder      = -1
    last_ladder = -1
    for i, phi in enumerate(phi_vals):
        if np.isnan(phi):
            ladder = -1
        else:
            ladder = getLadder(phi)
        # each ladder index can only be used once
        if ladder >= 0 and ladder == last_ladder:
            ladder += 1
        # max ladder is 11
        if ladder > 11:
            ladder = 11
        last_ladder = ladder
        ladder_nums.append(ladder)
    return ladder_nums

# get occupancy per ladder number
def getLadderOccupancy(ladder_nums, occ_phi_ring):
    len_nums    = len(ladder_nums)
    ladder_occ  = np.zeros(len_nums)
    for i, num in enumerate(ladder_nums):
        if num >= 0: 
            ladder_occ[num] = occ_phi_ring[i]
    #print("ladder_nums: {0}".format(ladder_nums))
    #print("ladder_occ: {0}".format(ladder_occ))
    return ladder_occ

def read_file(input_file_):
    return np.load(input_file_, allow_pickle=True, encoding='latin1')

def remake_arrays(input_arr_, root_output_name, csv_output_name):
    useWeightedAve = False
    fixPhi         = True

    print("Running to create output file: {0}".format(root_output_name))
    
    w_r_bins = 0.01

    # need z-binning corresponding to 1 roc
    w_z_bins = 52  # # of pixels in a roc

    # need phi-binning corresponding to 1 roc (maybe 2?)
    w_phi_bins = 80

    n_z_bins = int(3328 / w_z_bins)  # 3328 is number of pixels in a ladder row

    n_phi_bins = int(960 / w_phi_bins)  # 1440 is number of pixels around phi for all ladders, 960 for inner ladders

    inner_array = np.array([row for row in input_arr_ if not np.all(row==None)])
    cleaned_array = np.array([[x if x is not None else [0, np.nan, np.nan, np.nan] for x in row]
                              for row in inner_array])
    #r_min = np.nanmin(cleaned_array[:, :, 1])
    #r_max = np.nanmax(cleaned_array[:, :, 1])

    # separate pixels into groups corresponding to rocs in phi and z
    array_by_rocs = np.array([cleaned_array[j*w_phi_bins:(j+1)*w_phi_bins, i*w_z_bins:(i+1)*w_z_bins] for i in range(n_z_bins) for j in range(n_phi_bins)])

    roc_index = range(0, n_z_bins*n_phi_bins)

    # minus - 0-383
    # plus - 384-767
    occ = []
    x_array = []
    y_array = []
    z_array = []
    z_err_array = []
    phi_array = []
    phi_err_array = []
    r_array = []
    r_err_array = []

    z_hl = []
    z_err_hl = []
    phi_hl = []
    r_hl = []
    r_err_hl = []
    occ_hl = []

    occ_ring = []
    r_ring = []
    phi_ring = []
    z_ring = []

    n_ladders = 12
    n_rings   = 64

    for x in range(n_ladders):
        occ_hl.append([])
        r_hl.append([])
        r_err_hl.append([])
        phi_hl.append([])
        z_hl.append([])
        z_err_hl.append([])
       
    for x in range(n_rings):
        occ_ring.append([])
        r_ring.append([])
        phi_ring.append([])
        z_ring.append([])

    i_ladder = 0
    i_ring = 0

    # section off rocs into roc ladders (12 'ladders'), each true ladder is split in half 6 * 2 = 12
    for roc in roc_index:

        occ_tmp = np.concatenate(array_by_rocs[roc, :, :, 0])
        r       = np.concatenate(array_by_rocs[roc, :, :, 1])
        phi     = np.concatenate(array_by_rocs[roc, :, :, 2])
        z       = np.concatenate(array_by_rocs[roc, :, :, 3])

        x       = r[~np.isnan(z)] * np.cos(phi[~np.isnan(z)])
        y       = r[~np.isnan(z)] * np.sin(phi[~np.isnan(z)])
        r       = r[~np.isnan(z)]
        phi     = phi[~np.isnan(z)]
        occ_tmp = occ_tmp[~np.isnan(z)]
        z       = z[~np.isnan(z)]

        r = np.sqrt(x**2 + y**2 + z**2)

        # WARNING: phi = -pi = pi issue for roc % 12 == 3
        # - roc 3 in ladder 3 (roc % 12 == 3) crosses phi = -pi = +pi
        # - do not take a normal average
        # - split into phi < 0 and phi >= 0
        # - for phi < 0, find the different from -pi, then add this to +pi for the average
        # - make sure the final average is within [-pi, +pi]... if avg > pi, then it should be set to -pi < new_avg < 0
        if fixPhi and (roc % 12 == 3):
            phi_neg = [val for val in phi if val <  0.0]
            phi_pos = [val for val in phi if val >= 0.0]
            # for -pi < phi < 0, find absolute value of different from -pi, and add to +pi
            phi_neg_fixed = [np.pi + abs(-np.pi - val) for val in phi if val < 0.0]
            phi_fixed = phi_pos + phi_neg_fixed
            #print("roc {0}: z_avg = {1:.3f}, phi_avg = {2:.3f}, phi_fixed_avg = {3:.3f}".format(roc, np.average(z), np.average(phi), np.average(phi_fixed)))
        else:
            phi_fixed = phi

        #print("roc {0}: z_avg = {1:.3f}, phi_avg = {2:.3f}, phi_fixed_avg = {3:.3f}".format(roc, np.average(z), np.average(phi), np.average(phi_fixed)))
        
        occ.append(np.sum(occ_tmp))
        
        if useWeightedAve:
            # use weights
            x_array.append(np.average(x,            weights=occ_tmp))
            y_array.append(np.average(y,            weights=occ_tmp))
            z_array.append(np.average(z,            weights=occ_tmp))
            z_err_array.append(np.std(z))
            phi_array.append(np.average(phi_fixed,  weights=occ_tmp))
            phi_err_array.append(np.std(phi_fixed,  weights=occ_tmp))
            r_array.append(np.average(r,            weights=occ_tmp))
            r_err_array.append(np.std(r))

            occ_hl[i_ladder].append(np.sum(occ_tmp))
            r_hl[i_ladder].append(np.average(r,             weights=occ_tmp))
            r_err_hl[i_ladder].append(np.std(r))
            phi_hl[i_ladder].append(np.average(phi_fixed,   weights=occ_tmp))
            z_hl[i_ladder].append(np.average(z,             weights=occ_tmp))
            z_err_hl[i_ladder].append(np.std(z))

            occ_ring[i_ring].append(np.sum(occ_tmp))
            r_ring[i_ring].append(np.average(r,             weights=occ_tmp))
            phi_ring[i_ring].append(np.average(phi_fixed,   weights=occ_tmp))
            z_ring[i_ring].append(np.average(z,             weights=occ_tmp))
        
        else:
            # do not use weights
            x_array.append(np.average(x))
            y_array.append(np.average(y))
            z_array.append(np.average(z))
            z_err_array.append(np.std(z))
            phi_array.append(np.average(phi_fixed))
            phi_err_array.append(np.std(phi_fixed))
            r_array.append(np.average(r))
            r_err_array.append(np.std(r))

            occ_hl[i_ladder].append(np.sum(occ_tmp))
            r_hl[i_ladder].append(np.average(r))
            r_err_hl[i_ladder].append(np.std(r))
            phi_hl[i_ladder].append(np.average(phi_fixed))
            z_hl[i_ladder].append(np.average(z))
            z_err_hl[i_ladder].append(np.std(z))

            occ_ring[i_ring].append(np.sum(occ_tmp))
            r_ring[i_ring].append(np.average(r))
            phi_ring[i_ring].append(np.average(phi_fixed))
            z_ring[i_ring].append(np.average(z))
        
        i_ladder += 1
        if i_ladder == n_ladders:
            i_ring += 1
            i_ladder = 0

    occ = np.array(occ)
    x_array = np.array(x_array)
    y_array = np.array(y_array)
    z_array = np.array(z_array)
    z_err_array = np.array(z_err_array)
    phi_array = np.array(phi_array)
    phi_err_array = np.array(phi_err_array)
    r_array = np.array(r_array)
    r_err_array = np.array(r_err_array)

    phi_sort = np.argsort(phi_array)
    z_sort = np.argsort(z_array)

    occ = occ[z_sort]
    x_array = x_array[z_sort]
    y_array = y_array[z_sort]
    z_array = z_array[z_sort]
    z_err_array = z_err_array[z_sort]
    phi_array = phi_array[z_sort]
    phi_err_array = phi_err_array[z_sort]
    r_array = r_array[z_sort]
    r_err_array = r_err_array[z_sort]

    # removing rocs
    remove_z     =  (z_array > -25 ) * (z_array < 25)
    remove_blips =  (z_array < -21) + (z_array > -20)
    remove_blips *= (z_array < -14.5) + (z_array > -13.5)
    remove_blips *= (z_array < -7.5) + (z_array > -6.5)
    remove_blips *= (z_array < 5.75) + (z_array > 6.5)
    remove_blips *= (z_array < 12.5) + (z_array > 13.5)
    remove_blips *= (z_array < 19) + (z_array > 20)

    remove_hl =  (phi_array < -2.3) + (phi_array > -2)
    remove_hl *= (phi_array < -1.3) + (phi_array > -1)
    remove_hl *= (phi_array < -0.25) + (phi_array > 0)
    remove_hl *= (phi_array < 0.8) + (phi_array > 1.1)
    remove_hl *= (phi_array < 1.8) + (phi_array > 2.2)
    remove_hl *= (phi_array < 2.4) + (phi_array > 2.7)

    remove_hl_low =  (phi_array > -2.3) * (phi_array < -2)
    remove_hl_low *= (phi_array > -1.3) * (phi_array < -1)
    remove_hl_low *= (phi_array > -0.25) * (phi_array < 0)
    remove_hl_low *= (phi_array > 0.8) * (phi_array < 1.1)
    remove_hl_low *= (phi_array > 1.8) * (phi_array < 2.2)
    remove_hl_low *= (phi_array > 2.4) * (phi_array < 2.7)

    #occ = occ[remove_z*remove_blips*remove_hl_low]
    #x_array = x_array[remove_z*remove_blips*remove_hl_low]
    #y_array = y_array[remove_z*remove_blips*remove_hl_low]
    #z_array = z_array[remove_z*remove_blips*remove_hl_low]
    #phi_array = phi_array[remove_z*remove_blips*remove_hl_low]
    #r_array = r_array[remove_z*remove_blips*remove_hl_low]
    #occ = occ[remove_z*remove_blips*remove_hl]
    #x_array = x_array[remove_z*remove_blips*remove_hl]
    #y_array = y_array[remove_z*remove_blips*remove_hl]
    #z_array = z_array[remove_z*remove_blips*remove_hl]
    #phi_array = phi_array[remove_z*remove_blips*remove_hl]
    #r_array = r_array[remove_z*remove_blips*remove_hl]

    #occ = occ[remove_z*remove_blips]
    #x_array = x_array[remove_z*remove_blips]
    #y_array = y_array[remove_z*remove_blips]
    #z_array = z_array[remove_z*remove_blips]
    #phi_array = phi_array[remove_z*remove_blips]
    #r_array = r_array[remove_z*remove_blips]
    
    # remove z
    occ             = occ[              remove_z*remove_blips]
    x_array         = x_array[          remove_z*remove_blips]
    y_array         = y_array[          remove_z*remove_blips]
    z_array         = z_array[          remove_z*remove_blips]
    phi_array       = phi_array[        remove_z*remove_blips]
    r_array         = r_array[          remove_z*remove_blips]
    r_err_array     = r_err_array[      remove_z*remove_blips]
    z_err_array     = z_err_array[      remove_z*remove_blips]
    phi_err_array   = phi_err_array[    remove_z*remove_blips]
    
    # remove phi: should be after removing z
    remove_phi = (phi_array > -np.pi) * (phi_array < np.pi)
    
    occ             = occ[              remove_phi]
    x_array         = x_array[          remove_phi]
    y_array         = y_array[          remove_phi]
    z_array         = z_array[          remove_phi]
    phi_array       = phi_array[        remove_phi]
    r_array         = r_array[          remove_phi]
    r_err_array     = r_err_array[      remove_phi]
    z_err_array     = z_err_array[      remove_phi]
    phi_err_array   = phi_err_array[    remove_phi]

    phi_sort = np.argsort(phi_array)
    z_sort   = np.argsort(z_array)
    r_sort   = np.argsort(r_array)

    occ_r = occ[r_sort]
    x_r_array = x_array[r_sort]
    y_r_array = y_array[r_sort]
    z_r_array = z_array[r_sort]
    z_err_r_array = z_err_array[r_sort]
    phi_r_array = phi_array[r_sort]
    r_r_array = r_array[r_sort]
    r_err_r_array = r_err_array[r_sort]

    occ_z = occ[z_sort]
    x_z_array = x_array[z_sort]
    y_z_array = y_array[z_sort]
    z_z_array = z_array[z_sort]
    z_err_z_array = z_err_array[z_sort]
    phi_z_array = phi_array[z_sort]
    r_z_array = r_array[z_sort]
    r_err_z_array = r_err_array[z_sort]

    occ_phi = occ[phi_sort]
    x_phi_array = x_array[phi_sort]
    y_phi_array = y_array[phi_sort]
    z_phi_array = z_array[phi_sort]
    z_err_phi_array = z_err_array[phi_sort]
    phi_phi_array = phi_array[phi_sort]
    phi_err_phi_array = phi_err_array[phi_sort]
    r_phi_array = r_array[phi_sort]
    r_err_phi_array = r_err_array[phi_sort]

    ####### begin condensing for projections ###############
    z_condense = []
    z_err_condense = []
    occ_z_condense = []
    n_half_ladders = 12
    for iz, z in enumerate(z_z_array):
        if iz % n_half_ladders == 0:
            z_section = iz//n_half_ladders
            z_condense.append(np.average(z_z_array[iz:n_half_ladders*(z_section+1)], weights=occ_z[iz:n_half_ladders*(z_section+1)]))
            z_err_condense.append(np.sqrt(np.sum(z_err_z_array[iz:n_half_ladders*(z_section+1)]**2)))
            occ_z_condense.append(np.sum(occ_z[iz:n_half_ladders*(z_section+1)]))
        else:
            continue
    phi_condense = []
    phi_err_condense = []
    occ_phi_condense = []
    n_z_sections = 64
    for iphi, phi in enumerate(phi_phi_array):
        if iphi % n_z_sections == 0:
            phi_section = iphi//n_z_sections
            phi_condense.append(np.average(phi_phi_array[iphi:n_z_sections*(phi_section+1)], weights=occ_phi[iphi:n_z_sections*(phi_section+1)]))
            phi_err_condense.append(np.sqrt(np.sum(phi_err_phi_array[iphi:n_z_sections*(phi_section+1)]**2)))
            occ_phi_condense.append(np.sum(occ_phi[iphi:n_z_sections*(phi_section+1)]))
        else:
            continue
    r_condense = []
    r_err_condense = []
    occ_r_condense = []
    for ir, r in enumerate(r_r_array):
        if ir % n_half_ladders == 0:
            r_section = ir//n_half_ladders
            r_condense.append(np.average(r_r_array[ir:n_half_ladders*(r_section+1)], weights=occ_r[ir:n_half_ladders*(r_section+1)]))
            r_err_condense.append(np.sqrt(np.sum(r_err_r_array[ir:n_half_ladders*(r_section+1)]**2)))
            occ_r_condense.append(np.sum(occ_r[ir:n_half_ladders*(r_section+1)]))
        else:
            continue
    r_condense_comb = []
    r_err_condense_comb = []
    occ_r_condense_comb = []
    for ir, r in enumerate(r_condense):
        if ir % 2 == 0:
            if ir+1 < len(r_condense):
                r_condense_comb.append((r+r_condense[ir+1])/2)
                r_err_condense_comb.append(np.sqrt(r_err_condense[ir]**2+r_err_condense[ir+1]**2))
                occ_r_condense_comb.append(np.sum([occ_r_condense[ir]+occ_r_condense[ir+1]]))
            else:
                r_condense_comb.append(r)
                r_err_condense_comb.append(r_err_condense[ir])
                occ_r_condense_comb.append(occ_r_condense[ir])
        else:
            continue
    ########### end condensing for projections ################

    gr_r_hl = []
    gr_z_hl = []
    phi_avg_hl = np.array([np.mean(phi) for phi in phi_hl])
    avg_phi_sort = np.argsort(phi_avg_hl)
    r_hl = np.array(r_hl)
    r_err_hl = np.array(r_err_hl)
    z_hl = np.array(z_hl)
    z_err_hl = np.array(z_err_hl)
    occ_hl = np.array(occ_hl)
    r_hl = r_hl[avg_phi_sort]
    r_err_hl = r_err_hl[avg_phi_sort]
    z_hl = z_hl[avg_phi_sort]
    z_err_hl = z_err_hl[avg_phi_sort]
    occ_hl = occ_hl[avg_phi_sort]

    for hl in range(n_ladders):
        r_hl[hl] = np.array(r_hl[hl])
        r_err_hl[hl] = np.array(r_err_hl[hl])
        z_hl[hl] = np.array(z_hl[hl])
        z_err_hl[hl] = np.array(z_err_hl[hl])
        occ_z_hl = np.array(occ_hl[hl])
        occ_r_hl = np.array(occ_hl[hl])

        remove_z_hl = (z_hl[hl] > -12.5) * (z_hl[hl] < 12.5)
        remove_z_hl += (z_hl[hl] < -12.5) + (z_hl[hl] > 12.5)
        remove_blips_hl = (z_hl[hl] < -21) + (z_hl[hl] > -20)
        remove_blips_hl *= (z_hl[hl] < -14.5) + (z_hl[hl] > -13.5)
        remove_blips_hl *= (z_hl[hl] < -7.5) + (z_hl[hl] > -6.5)
        remove_blips_hl *= (z_hl[hl] < 5.75) + (z_hl[hl] > 6.5)
        remove_blips_hl *= (z_hl[hl] < 12.5) + (z_hl[hl] > 13.5)
        remove_blips_hl *= (z_hl[hl] < 19) + (z_hl[hl] > 20)
        z_new_hl = z_hl[hl][remove_z_hl*remove_blips_hl]
        r_new_hl = r_hl[hl][remove_z_hl*remove_blips_hl]
        r_err_new_hl = r_err_hl[hl][remove_z_hl*remove_blips_hl]
        z_err_new_hl = z_err_hl[hl][remove_z_hl*remove_blips_hl]
        occ_z_hl = occ_z_hl[remove_z_hl*remove_blips_hl]
        occ_r_hl = occ_r_hl[remove_z_hl*remove_blips_hl]

        r_sort = np.argsort(r_new_hl)
        z_sort = np.argsort(z_new_hl)
        r_new_hl = r_new_hl[r_sort]
        r_err_new_hl = r_new_hl[r_sort]
        occ_r_hl = occ_r_hl[r_sort]

        
        r_condense_hl = []
        r_err_condense_hl = []
        occ_r_condense_hl = []
        #print("length of r_new_hl: {0}".format(len(r_new_hl)))
        for i, r in enumerate(r_new_hl):
            if i % 2 == 0:
                #print("i = {0}, r = {1}".format(i, r))
                # check i + 1 is still valid for the list
                if i + 1 < len(r_new_hl):
                    r_condense_hl.append(np.mean([r_new_hl[i], r_new_hl[i+1]]))
                    r_err_condense_hl.append(np.sqrt(r_err_new_hl[i]**2 + r_new_hl[i+1]**2))
                    occ_r_condense_hl.append(np.sum([occ_r_hl[i]+occ_r_hl[i+1]]))
                # if i + 1 is not valid, don't use i + 1
                else:
                    r_condense_hl.append(r_new_hl[i])
                    r_err_condense_hl.append(r_err_new_hl[i])
                    occ_r_condense_hl.append(occ_r_hl[i])
            else:
                continue

        z_new_hl = z_new_hl[z_sort]
        z_err_new_hl = z_err_new_hl[z_sort]
        occ_z_hl = occ_z_hl[z_sort]

        gr_r_hl.append(rt.TGraph())
        rnp.fill_graph(gr_r_hl[hl], np.swapaxes([r_condense_hl, occ_r_condense_hl], 0, 1))
        gr_r_hl[hl].SetName("gr_r_occ_hl_"+str(hl))

        gr_z_hl.append(rt.TGraph())
        rnp.fill_graph(gr_z_hl[hl], np.swapaxes([z_new_hl, occ_z_hl], 0, 1))
        gr_z_hl[hl].SetName("gr_z_occ_hl_" + str(hl))

    # Remove NANs in 2D array
    for i in range(len(z_ring)):
        z_i = np.array(z_ring[i])
        mask = ~np.isnan(z_i)
        z_i_new = z_i[mask]
        z_ring[i] = z_i_new
    
    gr_phi_ring = []
    gr_phi_ring_subtracted = []
    gr_phi_ring_postcut = []
    ring_array = []
    phi_per_ring_arrary = []
    z_avg_ring = np.array([np.mean(z) for z in z_ring])
    avg_z_sort = np.argsort(z_avg_ring)
    
    #print("length z_ring: {0}".format(len(z_ring)))
    #print("length z_avg_ring: {0}".format(len(z_avg_ring)))
    #print("length avg_z_sort: {0}".format(len(avg_z_sort)))
    #for i in range(len(z_ring)):
    #    print("length of z[{0}]: {1}".format(i, len(z_ring[i])))
    #print("z_ring: {0}".format(z_ring))
    #print("z_avg_ring: {0}".format(z_avg_ring))
    #print("avg_z_sort: {0}".format(avg_z_sort))
    
    occ_ring = np.array(occ_ring)
    phi_ring = np.array(phi_ring)
    phi_ring = phi_ring[avg_z_sort]
    occ_ring = occ_ring[avg_z_sort]
    z_avg_ring_sorted = z_avg_ring[avg_z_sort]
        
    #print("length occ_ring: {0}".format(len(occ_ring)))
    #print("length phi_ring: {0}".format(len(phi_ring)))
    
    # using fixed occupancy cut
    useFixedCut = False
    #min_occupancy = 1
    min_occupancy = 20000
    onlyGoodRings = True
    phi_ring_sum = np.zeros(12)
    occ_phi_ring_subtracted_sum = np.zeros(12)
    num_good_rings = 0

    # 2D histograms
    h2d_occupancy       = rt.TH2F("h2d_occupancy",      "h2d_occupancy",      64, 0.0, 64.0, 12, 0.0, 12.0)
    h2d_occupancy_cut   = rt.TH2F("h2d_occupancy_cut",  "h2d_occupancy_cut",  64, 0.0, 64.0, 12, 0.0, 12.0)
    h2d_occupancy_mask  = rt.TH2F("h2d_occupancy_mask", "h2d_occupancy_mask", 64, 0.0, 64.0, 12, 0.0, 12.0)
    
    # output to csv file
    output_column_titles = ["index", "ring", "ladder", "phi", "occupancy", "mask"]
    with open(csv_output_name, 'w', newline='') as output_csv:
        output_writer = csv.writer(output_csv)
        output_writer.writerow(output_column_titles)
        index = 0
        for ring in range(n_rings):
            # sort by phi for phi distribution

            #print("ring {0}: z = {1}".format(ring, z_avg_ring_sorted[ring]))
            ring_array.append(ring)
            phi_ring[ring]  = np.array(phi_ring[ring])
            occ_phi_ring    = np.array(occ_ring[ring])
            phi_sort        = np.argsort(phi_ring[ring])
            phi_ring[ring]  = phi_ring[ring][phi_sort]
            occ_phi_ring    = occ_phi_ring[phi_sort]
            
            n_vals  = len(occ_phi_ring)
            if n_vals != n_ladders:
                print("ERROR for ring = {0}: there are {1} values, expected {2} values".format(ring, n_vals, n_ladders))
            avg     = np.mean(occ_phi_ring)
            std_dev = np.std(occ_phi_ring)
            delta   = 0.30 * avg
            #print("ring {0}: n_vals = {1}, avg = {2:.2f}, std_dev = {3:.2f}".format(ring, n_vals, avg, std_dev))

            # --- cut on occupancy --- #
            # default:  use fixed value for cut
            # variable: based on percent difference from average
            if not useFixedCut:
                min_occupancy = avg - delta
            occupancy_cut = occ_phi_ring >= min_occupancy
            
            occ_phi_ring_postcut  = occ_phi_ring[ occupancy_cut ]
            phi_ring_postcut      = phi_ring[ring][ occupancy_cut ]
            length_before_cut     = len(occ_phi_ring)
            length_after_cut      = len(occ_phi_ring_postcut)
            phi_per_ring_arrary.append(length_after_cut)
            #print("Ring {0}: num. points: before cut: {1}, after cut: {2}".format(ring, length_before_cut, length_after_cut))
            #print("phi_ring[{0}]: {1}".format(ring, phi_ring[ring]))
            #print("occ_phi_ring: {0}".format(occ_phi_ring))
            #print("occ_phi_ring_postcut: {0}".format(occ_phi_ring_postcut))

            # subtract average
            occ_phi_ring_subtracted = occ_phi_ring - avg
            
            # skip rings if there are NANs
            num_nans = len(phi_ring[ring][np.isnan(phi_ring[ring])])
            if num_nans == 0:
                #print(" --- good ring: {0}".format(ring))
                num_good_rings += 1
                phi_ring_sum = phi_ring_sum + phi_ring[ring]
            
            if onlyGoodRings:
                if num_nans == 0:
                    occ_phi_ring_subtracted_sum = occ_phi_ring_subtracted_sum + occ_phi_ring_subtracted
            else:
                occ_phi_ring_subtracted_sum = occ_phi_ring_subtracted_sum + occ_phi_ring_subtracted

            gr_phi_ring.append(rt.TGraph())
            rnp.fill_graph(gr_phi_ring[ring], np.swapaxes([phi_ring[ring], occ_phi_ring], 0, 1))
            gr_phi_ring[ring].SetName("gr_phi_occ_ring_standard_"+str(ring))

            # gr_phi_ring_subtracted
            gr_phi_ring_subtracted.append(rt.TGraph())
            rnp.fill_graph(gr_phi_ring_subtracted[ring], np.swapaxes([phi_ring[ring], occ_phi_ring_subtracted], 0, 1))
            gr_phi_ring_subtracted[ring].SetName("gr_phi_occ_ring_subtracted_"+str(ring))

            # fill after cut
            gr_phi_ring_postcut.append(rt.TGraph())
            rnp.fill_graph(gr_phi_ring_postcut[ring], np.swapaxes([phi_ring_postcut, occ_phi_ring_postcut], 0, 1))
            gr_phi_ring_postcut[ring].SetName("gr_phi_occ_ring_postcut_"+str(ring))
        
            # write to csv file and fill 2D histogram
            # be careful about missing phi points and phi = NAN
            ladder = 0
            ladderIndex = 0
            ladderFromPhi = 0
            ladder_nums = getLadderNums(phi_ring[ring])
            ladder_occ  = getLadderOccupancy(ladder_nums, occ_phi_ring)
            while ladder < n_ladders:
                phi = phi_ring[ring][ladder]
                occupancy = ladder_occ[ladder]
                # --- cut on occupancy --- #
                occupancy_after_cut = 0
                mask = 0
                # occupancy passing min occupancy cut
                if occupancy >= min_occupancy:
                    occupancy_after_cut = occupancy
                # mask:
                # 0 for points that pass selection
                # 1 for points that fail selection
                else:
                    mask = 1
                # set histo content
                h2d_occupancy.SetBinContent(ring+1, ladder+1, occupancy)
                h2d_occupancy_cut.SetBinContent(ring+1, ladder+1, occupancy_after_cut)
                h2d_occupancy_mask.SetBinContent(ring+1, ladder+1, mask)
                # write to csv file
                output_row = [index, ring, ladder, phi, occupancy, mask]
                output_writer.writerow(output_row)
                
                ladder += 1
                index  += 1

    # number of phi points after cut
    ring_array = np.array(ring_array)
    phi_per_ring_arrary = np.array(phi_per_ring_arrary)
    gr_num_phi_per_ring = rt.TGraph()
    rnp.fill_graph(gr_num_phi_per_ring, np.swapaxes([ring_array, phi_per_ring_arrary], 0, 1))
    gr_num_phi_per_ring.SetName("gr_num_phi_per_ring")

    phi_ring_avg = phi_ring_sum / num_good_rings

    #print("number of good rings: {0}".format(num_good_rings))
    #print("phi_condense (length {0}): {1}".format(len(phi_condense), phi_condense))
    #print("occ_phi_condense (length {0}): {1}".format(len(occ_phi_condense), occ_phi_condense))
    #print("phi_ring_avg (length {0}): {1}".format(len(phi_ring_avg), phi_ring_avg))
    #print("occ_phi_ring_subtracted_sum (length {0}): {1}".format(len(occ_phi_ring_subtracted_sum), occ_phi_ring_subtracted_sum))

    gr_phi = rt.TGraph()
    rnp.fill_graph(gr_phi, np.swapaxes([phi_condense, occ_phi_condense], 0, 1))
    gr_phi.SetName("gr_phi_occ")
    
    gr_phi_subtracted_sum = rt.TGraph()
    rnp.fill_graph(gr_phi_subtracted_sum, np.swapaxes([phi_ring_avg, occ_phi_ring_subtracted_sum], 0, 1))
    gr_phi_subtracted_sum.SetName("gr_phi_subtracted_sum_occ")

    gr_z = rt.TGraph()
    rnp.fill_graph(gr_z, np.swapaxes([z_condense, occ_z_condense], 0, 1))
    gr_z.SetName("gr_z_occ")

    gr_r3d = rt.TGraph()
    rnp.fill_graph(gr_r3d, np.swapaxes([r_condense, occ_r_condense], 0, 1))
    gr_r3d.SetName("gr_r_occ")
    
    gr_r3d = rt.TGraph()
    rnp.fill_graph(gr_r3d, np.swapaxes([r_condense_comb, occ_r_condense_comb], 0, 1))
    gr_r3d.SetName("gr_r_occ_pm_comb")

    # output ROOT file
    file_out = rt.TFile(root_output_name, "RECREATE")
    h2d_occupancy.Write()
    h2d_occupancy_cut.Write()
    h2d_occupancy_mask.Write()
    gr_phi.Write()
    gr_phi_subtracted_sum.Write()
    gr_z.Write()
    gr_r3d.Write()
    gr_num_phi_per_ring.Write()
    for hl in range(n_ladders):
        gr_r_hl[hl].Write()
        gr_z_hl[hl].Write()
    for ring in range(n_rings):
        gr_phi_ring[ring].Write()
        gr_phi_ring_subtracted[ring].Write()
        gr_phi_ring_postcut[ring].Write()
    file_out.Close()


if __name__ == "__main__":
    output_dir = "output"
    
    tools.makeDir(output_dir)
    inputs_v1 = [
        "TTBar_AllClusters_zsmear",
        #"TTBar_OnTrack_zsmear",
        #"TTBar_OffTrack_zsmear",
        #"TTBar_AllClusters",
        #"TTBar_OnTrack",
        #"TTBar_OffTrack",
        #"TTBar_0p2_AllClusters",
        #"TTBar_0p2_OnTrack",
        #"TTBar_0p2_OffTrack",
        #"TTBar_pileup_0p2_AllClusters",
        #"TTBar_pileup_0p2_OnTrack",
        #"TTBar_pileup_0p2_OffTrack",
        #"SingleMuon_AllClusters",
        #"SingleMuon_OnTrack",
        #"SingleMuon_OffTrack",
        #"ZeroBias_AllClusters",
        #"ZeroBias_OnTrack",
        #"ZeroBias_OffTrack",
        #"MinBias_2017B_AllClusters",
        #"MinBias_2017B_OnTrack",
        #"MinBias_2017B_OffTrack",
        #"MinBias_2018C_AllClusters",
        #"MinBias_2018C_OffTrack",
        #"MinBias_2018C_OnTrack",
    ]
    inputs_v2 = [
        "ZeroBias_2017B_AllClusters",
        "ZeroBias_2017B_ClusterSize2_AllClusters",
        "ZeroBias_2017B_ClusterSize2_NumberClusters2000_AllClusters",
    ]
    inputs_v3 = [
        "ZeroBias_2017B_MoreEvents_AllClusters",
        "ZeroBias_2017B_MoreEvents_ClusterSize2_AllClusters",
        "SingleMuon_2017B_MoreEvents_AllClusters",
        "SingleMuon_2017B_MoreEvents_ClusterSize2_AllClusters",
    ]
    inputs_v4 = [
        "ZeroBias_2017B_Legacy_MoreEvents_ClusterSize2_NumberClusters2000_AllClusters",
        "SingleMuon_2017B_Legacy_MoreEvents_ClusterSize2_NumberClusters2000_AllClusters",
    ]
    inputs_v5 = [
        "TTBar_AllClusters_zsmear",
        "ZeroBias_2017B_Legacy_MoreEvents_ClusterSize2_NumberClusters2000_AllClusters",
        "SingleMuon_2017B_Legacy_MoreEvents_ClusterSize2_NumberClusters2000_AllClusters",
    ]

    for sample in inputs_v5:
        in_array            = read_file("data/{0}.npy".format(sample))
        root_output_name    = "{0}/{1}.root".format(output_dir, sample)
        csv_output_name     = "{0}/{1}.csv".format(output_dir, sample)
        remake_arrays(in_array, root_output_name, csv_output_name)
 
