import numpy as np
import scipy.optimize as so
import scipy.stats as ss
import ROOT as rt
import root_numpy as rnp
import string
import pandas as pd
import os

alpha_low = string.ascii_lowercase

# creates directory if it does not exist
def makeDir(dir_name):
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)

def read_file(input_file_):
    return np.load(input_file_, allow_pickle=True, encoding='latin1')

def remake_arrays(input_arr_, file_out_name):
    useWeightedAve = False
    fixPhi         = True
    
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

    file_out = rt.TFile(file_out_name, "RECREATE")

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
    output_dir = "output"
    
    makeDir(output_dir)
    '''
    in_array    = read_file("data/TTBar_AllClusters_zsmear.npy")
    output_name = "{0}/TTBar_AllClusters_zsmear.root".format(output_dir)
    remake_arrays(in_array, output_name)

    in_array    = read_file("data/TTBar_OnTrack_zsmear.npy")
    output_name = "{0}/TTBar_OnTrack_zsmear.root".format(output_dir)
    remake_arrays(in_array, output_name)

    in_array    = read_file("data/TTBar_OffTrack_zsmear.npy")
    output_name = "{0}/TTBar_OffTrack_zsmear.root".format(output_dir)
    remake_arrays(in_array, output_name)
    
    in_array    = read_file("data/TTBar_AllClusters.npy")
    output_name = "{0}/TTBar_AllClusters.root".format(output_dir)
    remake_arrays(in_array, output_name)

    in_array    = read_file("data/TTBar_OnTrack.npy")
    output_name = "{0}/TTBar_OnTrack.root".format(output_dir)
    remake_arrays(in_array, output_name)

    in_array    = read_file("data/TTBar_OffTrack.npy")
    output_name = "{0}/TTBar_OffTrack.root".format(output_dir)
    remake_arrays(in_array, output_name)
    
    in_array    = read_file("data/SingleMuon_AllClusters.npy")
    output_name = "{0}/SingleMuon_AllClusters.root".format(output_dir)
    remake_arrays(in_array, output_name)
    
    in_array    = read_file("data/SingleMuon_OnTrack.npy")
    output_name = "{0}/SingleMuon_OnTrack.root".format(output_dir)
    remake_arrays(in_array, output_name)

    in_array    = read_file("data/SingleMuon_OffTrack.npy")
    output_name = "{0}/SingleMuon_OffTrack.root".format(output_dir)
    remake_arrays(in_array, output_name)
    
    in_array    = read_file("data/ZeroBias_AllClusters.npy")
    output_name = "{0}/ZeroBias_AllClusters.root".format(output_dir)
    remake_arrays(in_array, output_name)

    in_array    = read_file("data/ZeroBias_OnTrack.npy")
    output_name = "{0}/ZeroBias_OnTrack.root".format(output_dir)
    remake_arrays(in_array, output_name)

    in_array    = read_file("data/ZeroBias_OffTrack.npy")
    output_name = "{0}/ZeroBias_OffTrack.root".format(output_dir)
    remake_arrays(in_array, output_name)
    
    in_array    = read_file("data/TTBar_0p2_AllClusters.npy")
    output_name = "{0}/TTBar_0p2_AllClusters.root".format(output_dir)
    remake_arrays(in_array, output_name)

    in_array    = read_file("data/TTBar_0p2_OnTrack.py.npy")
    output_name = "{0}/TTBar_0p2_OnTrack.root".format(output_dir)
    remake_arrays(in_array, output_name)

    in_array    = read_file("data/TTBar_0p2_OffTrack.py.npy")
    output_name = "{0}/TTBar_0p2_OffTrack.root".format(output_dir)
    remake_arrays(in_array, output_name)

    in_array    = read_file("data/TTBar_pileup_0p2_AllClusters.npy")
    output_name = "{0}/TTBar_pileup_0p2_AllClusters.root".format(output_dir)
    remake_arrays(in_array, output_name)

    in_array    = read_file("data/TTBar_pileup_0p2_OnTrack.py.npy")
    output_name = "{0}/TTBar_pileup_0p2_OnTrack.root".format(output_dir)
    remake_arrays(in_array, output_name)

    in_array    = read_file("data/TTBar_pileup_0p2_OffTrack.py.npy")
    output_name = "{0}/TTBar_pileup_0p2_OffTrack.root".format(output_dir)
    remake_arrays(in_array, output_name)
   
    in_array    = read_file("data/MinBias_2017B_AllClusters.npy")
    output_name = "{0}/MinBias_2017B_AllClusters.root".format(output_dir)
    remake_arrays(in_array, output_name)
    
    in_array    = read_file("data/MinBias_2017B_OnTrack.npy")
    output_name = "{0}/MinBias_2017B_OnTrack.root".format(output_dir)
    remake_arrays(in_array, output_name)

    in_array    = read_file("data/MinBias_2017B_OffTrack.npy")
    output_name = "{0}/MinBias_2017B_OffTrack.root".format(output_dir)
    remake_arrays(in_array, output_name)
    '''
    in_array    = read_file("data/MinBias_2018C_AllClusters.npy")
    output_name = "{0}/MinBias_2018C_AllClusters.root".format(output_dir)
    remake_arrays(in_array, output_name)

    in_array    = read_file("data/MinBias_2018C_OffTrack.npy")
    output_name = "{0}/MinBias_2018C_OffTrack.root".format(output_dir)
    remake_arrays(in_array, output_name)

    in_array    = read_file("data/MinBias_2018C_OnTrack.npy")
    output_name = "{0}/MinBias_2018C_OnTrack.root".format(output_dir)
    remake_arrays(in_array, output_name)
    
