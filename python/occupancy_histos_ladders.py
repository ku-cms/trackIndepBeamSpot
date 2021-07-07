import numpy as np
#import probfit as pf
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
import scipy.optimize as so
import scipy.stats as ss
import ROOT as rt
import root_numpy as rnp

#import iminuit as im
#from pprint import pprint
import string
import pandas as pd

# design_0_nosmear_ge_2pix.npy
alpha_low = string.ascii_lowercase


def read_file(input_file_):
    return np.load(input_file_, allow_pickle=True, encoding='latin1')


def remake_arrays(input_arr_):

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
    r_min = np.nanmin(cleaned_array[:, :, 1])
    r_max = np.nanmax(cleaned_array[:, :, 1])

    # separate pixels into groups corresponding to rocs in phi and z
    array_by_rocs = np.array([cleaned_array[j*w_phi_bins:(j+1)*w_phi_bins, i*w_z_bins:(i+1)*w_z_bins] for i in range(n_z_bins) for j in range(n_phi_bins)])

    #roc_index = [0, 1]
    roc_index = range(0, n_z_bins*n_phi_bins)

    # fig, axs = plt.subplots(8, 8, sharex=False, sharey=False, figsize=(160, 160), tight_layout=True) #all rocs and modules
    # fig, axs = plt.subplots(12, 2, sharex=True, sharey=True, figsize=(20, 20), tight_layout=False) # fraction of rocs and modules
    # fig, axs = plt.subplots(1, sharex=True, sharey=True, figsize=(20, 20), tight_layout=True) # fraction of rocs and modules
    #fig = plt.figure() # fraction of rocs and modules
    #axs = fig.add_subplot(111, projection='3d') # fraction of rocs and modules

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
    n_rings = 64

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
        r = np.concatenate(array_by_rocs[roc, :, :, 1])
        phi = np.concatenate(array_by_rocs[roc, :, :, 2])
        z = np.concatenate(array_by_rocs[roc, :, :, 3])
        z_avg = np.nanmean(z)

        x = r[~np.isnan(z)] * np.cos(phi[~np.isnan(z)])
        y = r[~np.isnan(z)] * np.sin(phi[~np.isnan(z)])
        r = r[~np.isnan(z)]
        phi = phi[~np.isnan(z)]
        occ_tmp = occ_tmp[~np.isnan(z)]
        z = z[~np.isnan(z)]

        r = np.sqrt(x**2 + y**2 + z**2)

        occ.append(np.sum(occ_tmp))
        
        #x_array.append(np.average(x, weights=occ_tmp))
        #y_array.append(np.average(y, weights=occ_tmp))
        #z_array.append(np.average(z, weights=occ_tmp))
        #z_err_array.append(np.std(z))
        #phi_array.append(np.average(phi, weights=occ_tmp))
        #phi_err_array.append(np.average(phi, weights=occ_tmp))
        #r_array.append(np.average(r, weights=occ_tmp))
        #r_err_array.append(np.std(r))

        #occ_hl[i_ladder].append(np.sum(occ_tmp))
        #r_hl[i_ladder].append(np.average(r, weights=occ_tmp))
        #r_err_hl[i_ladder].append(np.std(r))
        #phi_hl[i_ladder].append(np.average(phi, weights=occ_tmp))
        #z_hl[i_ladder].append(np.average(z, weights=occ_tmp))
        #z_err_hl[i_ladder].append(np.std(z))

        #occ_ring[i_ring].append(np.sum(occ_tmp))
        #r_ring[i_ring].append(np.average(r, weights=occ_tmp))
        #phi_ring[i_ring].append(np.average(phi, weights=occ_tmp))
        #z_ring[i_ring].append(np.average(z, weights=occ_tmp))
        
        # remove weights
        x_array.append(np.average(x))
        y_array.append(np.average(y))
        z_array.append(np.average(z))
        z_err_array.append(np.std(z))
        phi_array.append(np.average(phi))
        phi_err_array.append(np.average(phi))
        r_array.append(np.average(r))
        r_err_array.append(np.std(r))

        occ_hl[i_ladder].append(np.sum(occ_tmp))
        r_hl[i_ladder].append(np.average(r))
        r_err_hl[i_ladder].append(np.std(r))
        phi_hl[i_ladder].append(np.average(phi))
        z_hl[i_ladder].append(np.average(z))
        z_err_hl[i_ladder].append(np.std(z))

        occ_ring[i_ring].append(np.sum(occ_tmp))
        r_ring[i_ring].append(np.average(r))
        phi_ring[i_ring].append(np.average(phi))
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
    remove_z =  (z_array > -12.5) * (z_array < 12.5)
    remove_z += (z_array < -12.5) + (z_array > 12.5)
    remove_z *= (z_array > -25 ) * (z_array < 25)
    remove_blips = (z_array < -21) + (z_array > -20)
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
    remove_phi = (phi_array > 3.0)
    
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
            r_condense_comb.append((r+r_condense[ir+1])/2)
            r_err_condense_comb.append(np.sqrt(r_err_condense[ir]**2+r_err_condense[ir+1]**2))
            occ_r_condense_comb.append(np.mean([occ_r_condense[ir], occ_r_condense[ir+1]]))
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
        print("length of r_new_hl: {0}".format(len(r_new_hl)))
        for i, r in enumerate(r_new_hl):
            if i % 2 == 0:
                #r_condense_hl.append((r + r_new_hl[i+1]) / 2)
                #r_err_condense_hl.append(np.sqrt(r_err_new_hl[i]**2 + r_new_hl[i+1]**2))
                #occ_r_condense_hl.append(np.mean([occ_r_hl[i], occ_r_hl[i+1]]))
                print("i = {0}, r = {1}".format(i, r))
                # check i + 1 is still valid for the list
                if i + 1 < len(r_new_hl):
                    r_condense_hl.append(np.mean([r_new_hl[i], r_new_hl[i+1]]))
                    r_err_condense_hl.append(np.sqrt(r_err_new_hl[i]**2 + r_new_hl[i+1]**2))
                    occ_r_condense_hl.append(np.mean([occ_r_hl[i], occ_r_hl[i+1]]))
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

    gr_phi_ring = []
    z_avg_ring = np.array([np.mean(z) for z in z_ring])
    avg_z_sort = np.argsort(z_avg_ring)
    occ_ring = np.array(occ_ring)
    phi_ring = np.array(phi_ring)
    phi_ring = phi_ring[avg_z_sort]
    occ_ring = occ_ring[avg_z_sort]

    for ring in range(n_rings):
        phi_ring[ring] = np.array(phi_ring[ring])
        occ_phi_ring = np.array(occ_ring[ring])
        phi_sort = np.argsort(phi_ring[ring])
        phi_ring[ring] = phi_ring[ring][phi_sort]
        occ_phi_ring = occ_phi_ring[phi_sort]

        gr_phi_ring.append(rt.TGraph())
        rnp.fill_graph(gr_phi_ring[ring], np.swapaxes([phi_ring[ring], occ_phi_ring], 0, 1))
        gr_phi_ring[ring].SetName("gr_phi_occ_ring_"+str(ring))


    gr_phi = rt.TGraph()
    rnp.fill_graph(gr_phi, np.swapaxes([phi_condense, occ_phi_condense], 0, 1))
    gr_phi.SetName("gr_phi_occ")

    gr_z = rt.TGraph()
    rnp.fill_graph(gr_z, np.swapaxes([z_condense, occ_z_condense], 0, 1))
    gr_z.SetName("gr_z_occ")

    gr_r3d = rt.TGraph()
    rnp.fill_graph(gr_r3d, np.swapaxes([r_condense, occ_r_condense], 0, 1))
    gr_r3d.SetName("gr_r_occ")
    gr_r3d = rt.TGraph()
    rnp.fill_graph(gr_r3d, np.swapaxes([r_condense_comb, occ_r_condense_comb], 0, 1))
    gr_r3d.SetName("gr_r_occ_pm_comb")

    #file_out = rt.TFile("output_0_charge_l200_nosmear.root", "RECREATE")
    #file_out = rt.TFile("output_0p1_neg0p08_charge_l200_nosmear.root", "RECREATE")
    #file_out = rt.TFile("output_neg0p1_0p2_charge_l200_nosmear.root", "RECREATE")
    #file_out = rt.TFile("output_0p2_0p19_charge_l200_nosmear.root", "RECREATE")
    #file_out = rt.TFile("output_neg0p3_neg0p32_charge_l200_nosmear.root", "RECREATE")
    #file_out = rt.TFile("output_0_charge_ge200_nosmear.root", "RECREATE")
    #file_out = rt.TFile("output_0p1_neg0p08_charge_ge200_nosmear.root", "RECREATE")
    #file_out = rt.TFile("output_neg0p1_0p2_charge_ge200_nosmear.root", "RECREATE")
    #file_out = rt.TFile("output_0p2_0p19_charge_ge200_nosmear.root", "RECREATE")
    #file_out = rt.TFile("output_neg0p3_neg0p32_charge_ge200_nosmear.root", "RECREATE")
    
    file_out = rt.TFile("test.root", "RECREATE")

    gr_phi.Write()
    gr_z.Write()
    gr_r3d.Write()
    for hl in range(n_ladders):
        gr_r_hl[hl].Write()
        gr_z_hl[hl].Write()
    for ring in range(n_rings):
        gr_phi_ring[ring].Write()
    file_out.Close()


if __name__ == "__main__":
    #in_array = read_file("design_0_nosmear_ge_2pix.npy")
    #in_array = read_file("design_0_no_outer_ge_2pix_smear.npy")
    #in_array = read_file("design_1p0_no_outer_ge_2pix_smear.npy")
    #in_array = read_file("design_0p1_0p08_no_outer_ge_2pix_smear.npy")
    #in_array = read_file("design_0p1_neg0p08_no_outer_ge_2pix_smear.npy")
    #in_array = read_file("design_z10_no_outer_ge_2pix.npy")
    #in_array = read_file("design_zneg10_no_outer_ge_2pix.npy")
    #in_array = read_file("design_0_nosmear_no_outer_ge_2pix.npy")
    #in_array = read_file("design_0_no_outer_ge_2pix_nosmear_phifix.npy")
    #in_array = read_file("design_z10_0_no_outer_ge_2pix_nosmear_phifix.npy")
    #in_array = read_file("design_0p1_no_outer_ge_2pix_nosmear_phifix.npy")
    #in_array = read_file("design_0_no_outer_ge_2pix_smear_phifix.npy")
    #in_array = read_file("design_0p1_no_outer_ge_2pix_smear_phifix.npy")
    #in_array = read_file("design_0_no_outer_all_pix_nosmear_phifix.npy")
    #in_array = read_file("design_0p1_no_outer_all_pix_nosmear_phifix.npy")
    #inrarray = read_file("design_0p2_no_outer_all_pix_nosmear.npy")
    #in_array = read_file("design_0p3_no_outer_all_pix_nosmear.npy")
    #in_array = read_file("design_0p1_II_no_outer_all_pix_nosmear.npy")
    #in_array = read_file("design_0_no_outer_all_pix_nosmear_chargege200.npy")
    #in_array = read_file("design_0p1_no_outer_all_pix_nosmear_chargege200.npy")
    #in_array = read_file("design_0p1_II_no_outer_all_pix_nosmear_chargege200.npy")
    #in_array = read_file("design_0p2_no_outer_all_pix_nosmear_chargege200.npy")
    #in_array = read_file("design_0p3_no_outer_all_pix_nosmear_chargege200.npy")
    #in_array = read_file("design_0_no_outer_all_pix_nosmear_chargel200.npy")
    #in_array = read_file("design_0p1_no_outer_all_pix_nosmear_chargel200.npy")
    #in_array = read_file("design_0p1_II_no_outer_all_pix_nosmear_chargel200.npy")
    #in_array = read_file("design_0p2_no_outer_all_pix_nosmear_chargel200.npy")
    #in_array = read_file("design_0p3_no_outer_all_pix_nosmear_chargel200.npy")
    
    #in_array = read_file("data/TTBar_AllClusters_zsmear.npy")
    in_array = read_file("data/SingleMuon_AllClusters.npy")
    remake_arrays(in_array)

