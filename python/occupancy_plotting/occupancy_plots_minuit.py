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

    array_by_rocs.argsort()
    roc_index = range(0, n_z_bins)

    # fig, axs = plt.subplots(8, 8, sharex=False, sharey=False, figsize=(160, 160), tight_layout=True) #all rocs and modules
    fig, axs = plt.subplots(1, 1, sharex=True, sharey=True, figsize=(20, 20), tight_layout=True) # fraction of rocs and modules
    x, y = 0, 0

    for roc in roc_index:
        roc_i = roc
        roc = roc_map(roc)
        # if (roc != 16) and (roc != 63) and (roc != 2): continue # looking at a single ROC
        # if roc != 31: continue # looking at a single ROC
        print(roc)
        occ = array_by_rocs[roc, :, :, 0]
        r = array_by_rocs[roc, :, :, 1]
        z = array_by_rocs[roc, :, :, 3]
        z_avg = np.nanmean(np.round(z, 2))

        r_sph = np.sqrt(r**2 + z**2)
        # if z_avg < 0:
        #    r_sph = np.negative(r_sph)
        r_sph_max = np.nanmax(r_sph)
        r_sph_min = np.nanmin(r_sph)
        n_r_sph_bins = int((r_sph_max - r_sph_min) / w_r_bins)

        r_min_roc = np.nanmin(r)
        r_max_roc = np.nanmax(r)
        n_r_bins_roc = int((r_max_roc - r_min_roc) / w_r_bins)

        r_range_min = (r_min_roc - r_min_roc)*100
        r_range_max = (r_max_roc - r_min_roc)*100
        r = (r - r_min_roc)*100

        # bin_values, bins, patches = axs[x][y].hist(np.concatenate(r), bins=n_r_bins_roc,
        #                                           range=(r_range_min, r_range_max),
        #                                           # weights=np.concatenate(occ), density=True)
        #                                           weights = np.concatenate(occ))
        bin_values, bins, patches = axs.hist(np.concatenate(r_sph), bins=n_r_sph_bins,
                                                   range=(r_sph_min, r_sph_max),
                                                   # weights=np.concatenate(occ), density=True)
                                                   weights = np.concatenate(occ))

        bin_centers = np.array([(bins[i] + bins[i+1]) / 2. for i, b in enumerate(bins) if b != bins[-1]])
        # bin_chi2 = pf.costfunc.BinnedChi2(func, np.concatenate(r), weights=np.concatenate(occ),
        #                                  bins=n_r_bins, bound=(r_min_roc, r_max_roc))
        # chi2 = pf.Chi2Regression(func, bin_centers[1:-5], bin_values[1:-5]) # for r_cyl
        # chi2 = pf.Chi2Regression(func, bin_centers[14:-12], bin_values[14:-12]) # deltaR = 0.01
        chi2 = pf.Chi2Regression(func, bin_centers[17:-16], bin_values[17:-16]) # deltaR = 0.02

        minuit = im.Minuit(chi2, a=1, b=0.5, c=1,
                           error_a=1, error_b=0.01, error_c=1,
                           fix_b=False,
                           limit_a=(None, None), limit_b=(0, 4), limit_c=(None, None),
                           errordef=1)
        print(minuit.get_param_states())
        minuit.migrad()
        print(minuit.get_param_states())
        print(minuit.get_fmin())

        # chi2.draw(minuit, axs[x][y])
        chi2.draw(minuit, axs)
        #axs[x][y].text(2.8, 15000, round(z_avg, 2))
        #axs[x][y].legend(fontsize='xx-small')

        x += 1
        if x == 8:
            x = 0
            y += 1
            if y == 8:
                x = 0
                y = 0

    plt.axis((None, None, 0, 20000))
    plt.show()


def func(x, a, b, c):
    return a*(1/x**b) + c




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
    in_array = read_file("design_0_nosmear_no_outer_ge_2pix.npy")
    remake_arrays(in_array)