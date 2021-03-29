# from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as so
import scipy.stats as ss

import iminuit

import string
import pandas as pd

# design_0_nosmear_ge_2pix.npy
alpha_low = string.ascii_lowercase


def read_file(input_file_):
    return np.load(input_file_, allow_pickle=True)


def remake_arrays(input_arr_):

    w_r_bins = 0.005

    # need z-binning corresponding to 1 roc
    w_z_bins = 52  # # of pixels in a roc

    n_z_bins = 3328 / w_z_bins  # 3328 is number of pixels in a ladder row

    inner_array = np.array([row for row in input_arr_ if not np.all(row==None)])
    cleaned_array = np.array([[x if x is not None else [0, np.nan, np.nan, np.nan] for x in row]
                              for row in inner_array])
    r_min = np.nanmin(cleaned_array[:, :, 1])
    r_max = np.nanmax(cleaned_array[:, :, 1])

    n_r_bins = int((r_max - r_min) / w_r_bins)

    array_by_rocs = np.array([cleaned_array[:, i*w_z_bins:(i+1)*w_z_bins] for i in range(n_z_bins)])

    roc_index = range(0, n_z_bins)
    # fig, axs = plt.subplots(8, 8, sharex=True, sharey=True, figsize=(20, 20), tight_layout=True)
    fig, axs = plt.subplots(1, 1, sharex=True, sharey=True, figsize=(20, 20), tight_layout=True)
    x, y = 0, 0

    for roc in roc_index:
        if roc is not max(roc_index)/2: continue

        occ = array_by_rocs[roc, :, :, 0]
        r = array_by_rocs[roc, :, :, 1]
        z = array_by_rocs[roc, :, :, 3]
        z_avg = np.nanmean(np.round(z, 2))
        bin_values, bins, patches = axs[x][y].hist(np.concatenate(r), bins=n_r_bins, range=(r_min, r_max),
                                                   weights=np.concatenate(occ), density=True)
                                                   # weights = np.concatenate(occ))

        bin_centers = np.array([(bins[i] + bins[i+1]) / 2. for i, b in enumerate(bins) if b != bins[-1]])

        # param_init = [1500., 2., -200.]
        # bounds_init = [[-np.inf, 0.0, -np.inf], [np.inf, 3, np.inf]]
        # popt, pcov = so.curve_fit(func, bin_centers[1:-5], bin_values[1:-5], p0=param_init, bounds=bounds_init,
        #                          verbose=1, max_nfev=10000)
        popt, pcov = so.curve_fit(func, bin_centers[1:-5], bin_values[1:-5])

        perr = np.sqrt(np.diag(np.abs(pcov)))

        chi_sq = ss.chisquare(bin_values[1:-5], func(bin_centers[1:-5], *popt))
        ndof = len(bin_centers[1:-5]) - len(popt)

        param_string = ''
        for param, letter, err in zip(popt, alpha_low, perr):
            param_string += '{}={:.2f} +/- {:.2f}\n'.format(letter, param, err)
        param_string += 'chi^2={:.4f}\n'.format(chi_sq[0])
        param_string += 'ndof={}'.format(ndof)

        axs[x][y].plot(bin_centers[1:-5], func(bin_centers[1:-5], *popt), 'r-',
                       label=param_string)

        axs[x][y].text(2.8, 15, round(z_avg, 2))
        axs[x][y].legend(fontsize='xx-small')

        # x += 1
        # if x == 8:
        #     x = 0
        #     y += 1
        #     if y == 8:
        #         x = 0
        #         y = 0

    plt.show()


def func(x, a, b, c, d, e):
    # return a/x**b + c
    return a/x**2 + b/x**2 + c/x**3 + d/x**4 + e


if __name__ == "__main__":
    in_array = read_file("design_0_nosmear_ge_2pix.npy")
    remake_arrays(in_array)