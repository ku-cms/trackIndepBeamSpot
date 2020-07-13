# from __future__ import print_function
import numpy as np
import probfit as pf
import matplotlib.pyplot as plt
import scipy.optimize as so
import scipy.stats as ss

import iminuit as im
from pprint import pprint
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

    n_z_bins = int(3328 / w_z_bins)  # 3328 is number of pixels in a ladder row

    inner_array = np.array([row for row in input_arr_ if not np.all(row==None)])
    cleaned_array = np.array([[x if x is not None else [0, np.nan, np.nan, np.nan] for x in row]
                              for row in inner_array])
    r_min = np.nanmin(cleaned_array[:, :, 1])
    r_max = np.nanmax(cleaned_array[:, :, 1])

    n_r_bins = int((r_max - r_min) / w_r_bins)

    array_by_rocs = np.array([cleaned_array[:, i*w_z_bins:(i+1)*w_z_bins] for i in range(n_z_bins)])

    #roc_index = [0, 1]
    roc_index = range(0, n_z_bins)

    # fig, axs = plt.subplots(8, 8, sharex=False, sharey=False, figsize=(160, 160), tight_layout=True) #all rocs and modules
    fig, axs = plt.subplots(1, sharex=True, sharey=True, figsize=(20, 20), tight_layout=True) # fraction of rocs and modules
    x, y = 0, 0

    occ = []
    r_sph = []
    z_array = []
    phi_array = []

    for roc in roc_index:
        true_roc = roc_map(roc)
        # if 23 < roc < 40: continue
        # if true_roc == 7: continue
        # if true_roc == 15: continue
        # if true_roc == 23: continue
        # if true_roc == 39: continue
        # if true_roc == 47: continue
        # if true_roc == 55: continue

        occ_tmp = np.concatenate(array_by_rocs[roc, :, :, 0])
        r = np.concatenate(array_by_rocs[roc, :, :, 1])
        phi = np.concatenate(array_by_rocs[roc, :, :, 2])
        z = np.concatenate(array_by_rocs[roc, :, :, 3])
        z_avg = np.nanmean(z)

        r_sph_tmp = np.sqrt(r**2 + z**2)
        z_array.append(np.average(z[~np.isnan(r_sph_tmp)], weights=occ_tmp[~np.isnan(r_sph_tmp)]))
        phi_array.append(np.average(phi[~np.isnan(r_sph_tmp)], weights=occ_tmp[~np.isnan(r_sph_tmp)]))
        r_sph.append(np.average(r_sph_tmp[~np.isnan(r_sph_tmp)], weights=occ_tmp[~np.isnan(r_sph_tmp)]))
        occ.append(np.nansum(occ_tmp))

    z_array = np.array(z_array)
    phi_array = np.array(phi_array)
    r_sph = np.array(r_sph)
    occ = np.array(occ)
    sort = np.argsort(r_sph)
    z_array = z_array[sort]
    phi_array = phi_array[sort]
    r_sph = r_sph[sort]
    occ = occ[sort]


    # for roc selection
    # removing outermost and innermost rocs
    remove_z = ((z_array > 1) + (z_array < -1)) * (z_array > -25) * (z_array < 25)
    #remove_z = (z_array > -22) * (z_array < 22)
    #remove_z *= (z_array < 14) + (z_array > 20)
    #remove_z *= (z_array > -14) + (z_array < -20)
    # remove_z = ((z_array > 11.5) + z_array < 8.5)) * (z_array > -25) * (z_array < 25)
    # remove_z = ((z_array > -8.5) + (z_array < -11.5)) * (z_array > -25) * (z_array < 25)
    remove_blips = (z_array < -21) + (z_array > -20)
    remove_blips *= (z_array < -14.5) + (z_array > -13.5)
    remove_blips *= (z_array < -7.5) + (z_array > -6.5)
    remove_blips *= (z_array < 5.75) + (z_array > 6.5)
    remove_blips *= (z_array < 12.5) + (z_array > 13.5)
    remove_blips *= (z_array < 19) + (z_array > 20)

    remove_hl = (phi_array < -2.3) + (phi_array > -2)
    remove_hl *= (phi_array < -1.3) + (phi_array > -1)
    remove_hl *= (phi_array < -0.25) + (phi_array > 0)
    remove_hl *= (phi_array < 0.8) + (phi_array > 1.1)
    remove_hl *= (phi_array < 1.8) + (phi_array > 2.2)
    remove_hl *= (phi_array < 2.4) + (phi_array > 2.7)
    #remove_hl = (phi_array > -2.3) * (phi_array < -2)
    #remove_hl += (phi_array > -1.3) * (phi_array < -1)
    #remove_hl += (phi_array > -0.25) * (phi_array < 0)
    #remove_hl += (phi_array > 0.8) * (phi_array < 1.1)
    #remove_hl += (phi_array > 1.8) * (phi_array < 2.2)
    #remove_hl += (phi_array > 2.4) * (phi_array < 2.7)

    #z_array = z_array[remove_z*remove_blips*remove_hl]
    #r_sph = r_sph[remove_z*remove_blips*remove_hl]
    #occ = occ[remove_z*remove_blips*remove_hl]

    #z_array = z_array[remove_z*remove_blips]
    #r_sph = r_sph[remove_z*remove_blips]
    #occ = occ[remove_z*remove_blips]
    # removing rocs at the edge of the modules
    #mask = np.array([1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1], dtype=bool)
    #r_sph = r_sph[mask]
    #occ = occ[mask]
    r_sph_condense = []
    occ_condense = []
    for ir, r in enumerate(r_sph):
        if ir % 2 == 0:
            r_sph_condense.append((r+r_sph[ir+1])/2)
            occ_condense.append(occ[ir]+occ[ir+1])
            #r_sph_condense.append(r)
            #occ_condense.append(occ[ir])
        else:
            continue
    occ_z = occ
    r_sph = np.array(r_sph_condense)
    occ = np.array(occ_condense)

    #axs.plot(r_sph, occ, 'b*', label='all z')
    axs.plot(z_array, occ_z, 'b*', label='all z')

    #bin_centers = np.array([(bins[i] + bins[i+1]) / 2. for i, b in enumerate(bins) if b != bins[-1]])

    # bin_chi2 = pf.costfunc.BinnedChi2(func, np.concatenate(r), weights=np.concatenate(occ),
    #                                  bins=n_r_bins, bound=(r_min_roc, r_max_roc))
    #def chi2(a, b, c, d, e, f):
    #    func_array = (func(r_sph, z_array, a, b, c, d, e, f) - occ)**2
    #    return np.sum(func_array)

    chi2 = pf.Chi2Regression(func, r_sph, occ)

    minuit = im.Minuit(chi2, a=1, b=1, c=1, d=1, e=1, f=1,
                       error_a=0.01, error_b=0.01, error_c=1, error_d=0.01, error_e=0.01, error_f=0.01,
                       fix_b=False,
                       limit_a=(None, None), limit_b=(None, None), limit_c=(None, None), limit_d=(None, None),
                       limit_e=(None, None), limit_f=(None, None),
                       errordef=1)
    #minuit = im.Minuit(chi2, a=2300000, b=2, c=20000,
    #                         error_a=1, error_b=0.01, error_c=1,
    #                         fix_b=False,
    #                         limit_a=(None, None), limit_b=(0, 10), limit_c=(None, None),
    #                         errordef=1)
    minuit.migrad()
    print(minuit.get_param_states())
    print(minuit.get_fmin())

    #chi2.draw(minuit, axs)
    param_string = ''
    #for p in minuit.parameters:
    #    param_string += '{}: {} +/- {}\n'.format(p, np.format_float_scientific(minuit.values[p], precision=3),
    #                                             np.format_float_scientific(minuit.errors[p], precision=1))
    #axs.text(15, 5000, param_string, fontsize='xx-small')
    #axs.plot(r_sph, func(r_sph, z_array, *minuit.values.values()), color='red', label='fit')

    axs.legend()

    plt.show()


#def func(x, a, b, c, r0):
#    return a*(1/(x-r0)**(2*b)) + c

#def func(x, a, b, c, d, e, f):
#    return a*np.exp(-(x-b)**2/(2*c**2)) + d*(1/x**e) + f

def func(x, a, b, c, d, e, f, g):
    if x <= 16:
        return a*(1/x**b) + c
    else:
        return d*x**3 + e*x**2 + f*x + g
        #return d*(1/x**e) + f
#return a*(1/(x)**(b)) + c
    #return a*x**3 + b*x**2 + c*x +d


def roc_map(roc_num):
    base_index = {'0': 7,
                  '1': 6,
                  '2': 5,
                  '3': 4,
                  '4': 3,
                  '5': 2,
                  '6': 1,
                  '7': 0}
    roc_floor = roc_num // 8
    roc_mod = roc_num % 8
    return base_index[str(roc_mod)] + (roc_floor*8)


if __name__ == "__main__":
    #in_array = read_file("design_0_ge_2pix_nosmear.npy")
    #in_array = read_file("design_0_no_outer_ge_2pix_smear.npy")
    #in_array = read_file("design_z10_no_outer_ge_2pix.npy")
    #in_array = read_file("design_zneg10_no_outer_ge_2pix.npy")
    #in_array = read_file("design_0_nosmear_no_outer_ge_2pix.npy")
    #in_array = read_file("design_0_no_outer_ge_2pix_nosmear_phifix.npy")
    #in_array = read_file("design_0_no_outer_ge_2pix_smear_phifix.npy")
    #in_array = read_file("design_0_no_outer_ge_2pix_smear_phifix.npy")
    #in_array = read_file("design_0_no_outer_all_pix_nosmear_charge.npy")
    #in_array =read_file("design_0_no_outer_all_pix_nosmear_charge10.npy")
    #in_array = read_file("design_0_no_outer_all_pix_nosmear_edge.npy")
    in_array = read_file("design_0_no_outer_all_pix_nosmear_phifix.npy")
    #in_array = read_file("design_0p1_no_outer_all_pix_nosmear_charge.npy")
    #in_array = read_file("design_0_no_outer_all_pix_nosmear_chargel250.npy")
    #in_array = read_file("design_0_no_outer_all_pix_nosmear_chargel200.npy")
    #in_array = read_file("design_0_no_outer_all_pix_nosmear_chargege200.npy")
    remake_arrays(in_array)
