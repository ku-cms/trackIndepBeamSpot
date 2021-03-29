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
    fig, axs = plt.subplots(1, 3, sharex=True, sharey=True, figsize=(20, 20), tight_layout=True) # fraction of rocs and modules
    x, y = 0, 0

    occ_plus = []
    occ_minus = []
    r_sph_plus = []
    r_sph_minus = []

    for roc in roc_index:
        occ_tmp = np.concatenate(array_by_rocs[roc, :, :, 0])
        r = np.concatenate(array_by_rocs[roc, :, :, 1])
        z = np.concatenate(array_by_rocs[roc, :, :, 3])
        z_avg = np.nanmean(z)

        r_sph_tmp = np.sqrt(r**2 + z**2)
        if roc < 32:
            r_sph_minus.append(np.nanmean(r_sph_tmp))
            occ_minus.append(np.nansum(occ_tmp))
        else:
            r_sph_plus.append(np.nanmean(r_sph_tmp))
            occ_plus.append(np.nansum(occ_tmp))

    r_sph_minus = np.array(r_sph_minus)
    occ_minus = np.array(occ_minus)
    sort_minus = np.argsort(r_sph_minus)
    r_sph_minus = r_sph_minus[sort_minus]
    occ_minus = occ_minus[sort_minus]

    r_sph_plus = np.array(r_sph_plus)
    occ_plus = np.array(occ_plus)
    sort_plus = np.argsort(r_sph_plus)
    r_sph_plus = r_sph_plus[sort_plus]
    occ_plus = occ_plus[sort_plus]

    # removing rocs
    r_sph_plus = r_sph_plus[3:-2]
    occ_plus = occ_plus[3:-2]
    r_sph_minus = r_sph_minus[3:-2]
    occ_minus = occ_minus[3:-2]
    ## removing rocs at the edge of the modules
    #mask_plus = np.array([1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1], dtype=bool)
    #r_sph_plus = r_sph_plus[mask_plus]
    #occ_plus = occ_plus[mask_plus]
    #mask_minus = np.array([1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1], dtype=bool)
    #r_sph_minus = r_sph_minus[mask_minus]
    #occ_minus = occ_minus[mask_minus]

    axs[0].plot(r_sph_minus, occ_minus, 'b*', label='minus z')

    axs[1].plot(r_sph_plus, occ_plus, 'b*', label='plus z')

    axs[2].plot(r_sph_minus, occ_minus, 'r*', label='occ minus')
    axs[2].plot(r_sph_plus, occ_plus, 'b*', label='occ plus')

    #bin_centers = np.array([(bins[i] + bins[i+1]) / 2. for i, b in enumerate(bins) if b != bins[-1]])

    # bin_chi2 = pf.costfunc.BinnedChi2(func, np.concatenate(r), weights=np.concatenate(occ),
    #                                  bins=n_r_bins, bound=(r_min_roc, r_max_roc))
    chi2_minus = pf.Chi2Regression(func, r_sph_minus, occ_minus)
    chi2_plus = pf.Chi2Regression(func, r_sph_plus, occ_plus)

    minuit_minus = im.Minuit(chi2_minus, a=2300000, b=2, c=20000,
                             error_a=1, error_b=0.01, error_c=1,
                             fix_b=False,
                             limit_a=(None, None), limit_b=(0, 10), limit_c=(None, None),
                             errordef=1)
    minuit_plus = im.Minuit(chi2_plus, a=2300000, b=2, c=20000,
                            error_a=1, error_b=0.01, error_c=1,
                            fix_b=False,
                            limit_a=(None, None), limit_b=(0, 10), limit_c=(None, None),
                            errordef=1)
    #minuit_minus = im.Minuit(chi2_minus, a=2300000, b=0.5939, c=20000, r0=0,
    #                         error_a=1, error_b=0.01, error_c=1, error_r0=0.01,
    #                         fix_b=True,
    #                         limit_a=(None, None), limit_b=(0, 10), limit_c=(None, None), limit_r0=(None, None),
    #                         errordef=1)
    #minuit_plus = im.Minuit(chi2_plus, a=2300000, b=0.6025, c=20000, r0=0,
    #                        error_a=1, error_b=0.01, error_c=1, error_r0=0.01,
    #                        fix_b=True,
    #                        limit_a=(None, None), limit_b=(0, 10), limit_c=(None, None), limit_r0=(None, None),
    #                        errordef=1)
    #print(minuit.get_param_states())
    minuit_minus.migrad()
    minuit_plus.migrad()
    print(minuit_plus.get_param_states())
    print(minuit_plus.get_fmin())
    print(minuit_minus.get_param_states())
    print(minuit_minus.get_fmin())

    chi2_minus.draw(minuit_minus, axs[0])
    chi2_plus.draw(minuit_plus, axs[1])

    axs[2].plot(r_sph_minus, func(r_sph_minus, *minuit_minus.values.values()), color='red', label='fit minus')
    axs[2].plot(r_sph_plus, func(r_sph_plus, *minuit_plus.values.values()), color='blue', linestyle='dashed', label='fit plus')

    axs[2].legend()
    axs[1].legend()
    axs[0].legend()
    axs[2].grid()

    plt.show()


def func(x, a, b, c):
    return a*(1/(x)**(b)) + c


#def func(x, a, b, c):
#    return a*(1/(x)**(2*b)) + c


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
    #in_array = read_file("design_0_nosmear_ge_2pix.npy")
    #in_array = read_file("design_0_no_outer_ge_2pix_smear.npy")
    #in_array = read_file("design_z10_no_outer_ge_2pix.npy")
    #in_array = read_file("design_zneg10_no_outer_ge_2pix.npy")
    #in_array = read_file("design_0_nosmear_no_outer_ge_2pix.npy")
    #in_array = read_file("design_0_ge_2pix_nosmear.npy")
    remake_arrays(in_array)