import numpy as np
import probfit as pf
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
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
    phi_array = []
    r_array = []

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

        occ.append(np.sum(occ_tmp))
        x_array.append(np.average(x, weights=occ_tmp))
        y_array.append(np.average(y, weights=occ_tmp))
        z_array.append(np.average(z, weights=occ_tmp))
        phi_array.append(np.average(phi, weights=occ_tmp))
        r_array.append(np.average(r, weights=occ_tmp))

    occ = np.array(occ)
    x_array = np.array(x_array)
    y_array = np.array(y_array)
    z_array = np.array(z_array)
    phi_array = np.array(phi_array)
    r_array = np.array(r_array)

    phi_sort = np.argsort(phi_array)
    z_sort = np.argsort(z_array)

    occ = occ[z_sort]
    x_array = x_array[z_sort]
    y_array = y_array[z_sort]
    z_array = z_array[z_sort]
    phi_array = phi_array[z_sort]
    r_array = r_array[z_sort]

    #for iphi, phi in enumerate(phi_array):
    #    print(r_array[iphi], x_array[iphi], y_array[iphi], phi_array[iphi])
    # removing rocs
    #remove_z = ((z_array > 2) + (z_array < -2)) * (z_array > -25) * (z_array < 25)
    remove_z = (z_array > -25) * (z_array < 25)
    #remove_z = (z_array > -16) * (z_array < 16)
    #remove_z = (z_array < -16) + (z_array > 16)
    #remove_z = ((z_array > 11.5) + (z_array < 8.5)) * (z_array > -25) * (z_array < 25)
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

    remove_hl_low = (phi_array > -2.3) * (phi_array < -2)
    remove_hl_low += (phi_array > -1.3) * (phi_array < -1)
    remove_hl_low += (phi_array > -0.25) * (phi_array < 0)
    remove_hl_low += (phi_array > 0.8) * (phi_array < 1.1)
    remove_hl_low += (phi_array > 1.8) * (phi_array < 2.2)
    remove_hl_low += (phi_array > 2.4) * (phi_array < 2.7)

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
    occ_g0 = occ > 0
    occ = occ[remove_z*remove_blips*occ_g0]
    x_array = x_array[remove_z*remove_blips*occ_g0]
    y_array = y_array[remove_z*remove_blips*occ_g0]
    z_array = z_array[remove_z*remove_blips*occ_g0]
    phi_array = phi_array[remove_z*remove_blips*occ_g0]
    r_array = r_array[remove_z*remove_blips*occ_g0]

    ######## begin condensing for projections ###############
    #z_condense = []
    #occ_condense = []
    #n_half_ladders = 12
    #for iz, z in enumerate(z_array):
    #    if iz % n_half_ladders == 0:
    #        z_section = iz//n_half_ladders
    #        z_condense.append(np.average(z_array[iz:n_half_ladders*(z_section+1)], weights=occ[iz:n_half_ladders*(z_section+1)]))
    #        occ_condense.append(np.sum(occ[iz:n_half_ladders*(z_section+1)]))
    #    else:
    #        continue
    # phi_condense = []
    # occ_condense = []
    # n_z_sections = 64
    # print(phi_array)
    # for iphi, phi in enumerate(phi_array):
    #     if iphi % n_z_sections == 0:
    #         phi_section = iphi//n_z_sections
    #         phi_condense.append(np.mean(phi_array[iphi:n_z_sections*(phi_section+1)]))
    #         occ_condense.append(np.sum(occ[iphi:n_z_sections*(phi_section+1)]))
    #     else:
    #         continue
    # phi_array = phi_condense
    # occ = occ_condense
    ############ end condensing for projections ################

    # print(z_array)


    # removing rocs at the edge of the modules
    # mask_plus = np.array([1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1],
    #                     dtype=bool)
    # r_sph_plus = r_sph_plus[mask_plus]
    # occ_plus = occ_plus[mask_plus]

    #x_array = np.concatenate(x_clean)
    #y_array = np.concatenate(y_clean)
    #z_array = np.concatenate(z_clean)
    #occ = np.concatenate(occ_clean)

    #def nll(x0, y0, z0, n, b1, b2, b3, a1a, a1b, a1c, a3a, a3b, a3c, c1a, c1c, c3a, c3c):
    #def nll(x0, y0, z0, n, b1, b2, b3, a1a, a1b, a1c, a3a, a3b, a3c, c1, c3):
    def nll(x0, y0, z0, n, b1, b2, b3, a1, a3, c1, c3):

        # func_tmp = 0
        ri = np.float64(np.sqrt((x_array - x0) ** 2 + (y_array - y0) ** 2 + (z_array - z0) ** 2))
        phi_cor = np.arctan2(y_array-y0, x_array-x0)
        #a = np.float64(a_par(phi_cor, np.abs(z_array-z0), a1a, a1b, a1c, b2, a3a, a3b, a3c))
        a = np.float64(b_par(phi_cor, a1, b2, a3))
        b = np.float64(b_par(phi_cor, b1, b2, b3))
        #c = np.float64(a_par(phi_cor, np.abs(z_array-z0), c1a, a1b, c1c, b2, c3a, a3b, c3c))
        c = np.float64(b_par(phi_cor, c1, b2, c3))

        func_array = -np.log(1 / (np.sqrt(2 * np.pi * occ))) + (n * func(ri, a, b, c) - occ) ** 2 / (2 * np.pi * occ)

        return np.sum(func_array)

    # axs.plot(phi_array, occ, 'b*')
    # axs.plot(z_condense, occ_condense, 'b*')
    # axs.plot(z_array, phi_array, occ, 'b*')
    # bin_values, bins, patches = axs.hist(z_array, bins=64, range=(np.min(z_array),np.max(z_array)), weights=occ, density=False)
    # patches = axs.hist2d(z_array, phi_array, bins=(62,12), weights=occ, density=False)

    # axs[2].plot(r_sph_minus, occ_minus, 'r*', label='occ minus')
    # axs[2].plot(r_sph_plus, occ_plus, 'b*', label='occ plus')

    # chi2_minus = pf.Chi2Regression(func, r_sph_minus, occ_minus)
    # chi2_plus = pf.Chi2Regression(func, r_sph_plus, occ_plus)

    #minuit = im.Minuit(nll, x0=0, y0=0, z0=0.019, n=1,
    #                    error_x0=0.001, error_y0=0.001, error_z0=0.001, error_n=0.01,
    #                    fix_z0=False,
    #                    limit_x0=(-10, 10), limit_y0=(-10, 10), limit_z0=(-25, 25), limit_n=(None, None),
    #                    errordef=1)
    minuit = im.Minuit(nll, x0=0, y0=0.0, z0=-0.019,
                       n=1,
                       b1=0.01, b2=-2.2, b3=1.182,
                       a1a=1e3, a1b=-0.231, a1c=1e3,
                       a3a=1.7e4, a3b=-0.2118, a3c=3000,
                       c1=190, c3=1645,
                       #c1a=790, c1c=67,
                       #c3a=1e4, c3c=2000,
                       error_x0=0.001, error_y0=0.001, error_z0=0.001,
                       error_n=0.01,
                       error_b1=0.01, error_b2=0.01, error_b3=0.01,
                       #fix_x0=True, fix_y0=True, fix_z0=True,
                       #error_a1=10,
                       #error_a3=0.1,
                       error_c1=1,
                       error_c3=1,
                       fix_n=True,
                       #fix_b1=True, fix_b2=True, fix_b3=True,
                       #fix_a1a=False, fix_a1c=False,
                       #fix_a1b=False,
                       #fix_a3a=True, fix_a3c=True,
                       #fix_a3b=False,
                       #fix_c1a=False, fix_c1c=False,
                       #fix_c3a=True, fix_c3c=True,
                       #fix_c1b=True,
                       #fix_c3b=True,
                       limit_x0=(-1, 1), limit_y0=(-1, 1), limit_z0=(-1, 1),
                       # limit_b1=(0, 2),
                       limit_b2=(-np.pi, np.pi), limit_b3=(1, 3),
                       #limit_a1a=(0.01, None), limit_a1b=(-1, 0), limit_a1c=(0.01, None),
                       #limit_a3a=(0.01, None), limit_a3b=(-1, 0), limit_a3c=(0.01, None),
                       #limit_c1=(50, None), limit_c3=(0.01, None),
                       #limit_c1a=(0.01, None), limit_c1c=(0.01, None),
                       #limit_c3a=(0.01, None), limit_c3c=(0.01, None),
                       # limit_a1=(0., 1e6),
                       # limit_c1=(0, 1e4),
                       errordef=1)

    #minuit = im.Minuit(nll, x0=0, y0=0, z0=0.0, n=1, b1=0.0, b2=0.0, b3=1.167, a1=0.0, a2=0.0, a3=1.49e6, c1=0.0, c2=0.0, c3=3.7e4,
    #                   error_x0=0.001, error_y0=0.001, error_z0=0.001,
    #                   fix_z0=False,
    #                   fix_a1=True, fix_b1=False, fix_c1=False,
    #                   fix_a2=True, fix_b2=False, fix_c2=False,
    #                   limit_x0=(None, None), limit_y0=(None, None), limit_z0=(None, None),
    #                   limit_a1=(None, None), limit_a2=(-np.pi, np.pi), limit_a3=(0, None),
    #                   limit_b1=(None, None), limit_b2=(-np.pi, np.pi), limit_b3=(0, None),
    #                   limit_c1=(None, None), limit_c2=(-np.pi, np.pi), limit_c3=(0, None),
    #                   errordef=1)

    # print(minuit.get_param_states())
    minuit.set_strategy(2)
    minuit.migrad()
    #minuit.hesse()
    #minuit.minos()
    print(minuit.get_param_states())
    print(minuit.get_fmin())

    # axs.legend()

    # plt.show()


def a_par(x, z, a1a, a1b, a1c, a2, a3a, a3b, a3c):
    #return (a1a*np.exp(a1b*z)+a1c)*np.sin(x-a2)+(a3a*np.exp(a3b*z)+a3c)
    return (a1a/z**a1b+a1c)*np.sin(x-a2)+(a3a/z**a3b)+a3c


def b_par(x, b1=0.0, b2=0.0, b3=1.25e6):  #, b4=0.0, b5=0.0, b6=1.167):
    return b1*np.sin(x-b2)+b3


def c_par(x, c1=0.0, c2=0.0, c3=1.25e6):  #, b4=0.0, b5=0.0, b6=1.167):
    return c1*np.sin(x-c2)+c3


def func(x, a, b, c):
    return a * (1 / x ** b) + c


#def func(x, a=-0.8875e5, b=-1.439, c=1.377e7, d=2691, e=-1.618e5, f=2.646e6, g=-2.902e6):
#def func(x, a=-1.015e6, b=-0.684, c=1.626e7, d=0.867e4, e=-0.563e6, f=1.128e7, g=-0.628e8):
#if x <= 16:
#    #    return a*(1/x**b) + c
#    #else:
#    #    return d*x**3 + e*x**2 + f*x + g
#    return np.array([a*(1/x_val**b)+c if x_val <= 16 else d*x_val**3 + e*x_val**2 + f*x_val + g for x_val in x])

#def func(x, a=1.74e5, b=0.2835, c=6.12, d=-1.4e4, e=-0.5622, f=1.561e5):
#def func(x, a=1.74e5, b=0.2835, c=6.12, d=-1.4e4, e=-0.562, f=1.561e5):  # func high occupancy
#    return a*np.exp(-(x-b)**2/(2*c**2)) + d*(1/x**e) + f

#def func(x, a=2.042e6, b=0.9912, c=4.457e4):
#    return a*(1/x**b) + c

#def func(x, a=1.479e6, b=1.162, c=3.378e4):
#    return a*(1/x**b) + c


#def func(x, a=3.839e5, b=-7.195, c=9.847, d=7.152e4):
#    return a*np.exp(-(x-b)**2/(2*c**2)) + d


def func_nll(x, x0, y0, z0, n):
    # func_tmp = 0
    ri = np.sqrt((x[:, 1]-x0)**2 + (x[:, 2]-y0)**2 + (x[:, 3]-z0)**2)

    func_array = -np.log(1/(x[:, 0]*np.sqrt(2*np.pi))) + (n*func(ri) - x[:, 0])**2 / (2*x[:, 0])

    # for ix, x_val in enumerate(x):
    #     ri = np.sqrt((x_val[1]-x0)**2 + (x_val[2]-y0)**2 + (x_val[3]-z0)**2)
    #     func_tmp += -np.log(1/(x_val[0]*np.sqrt(2*np.pi))) + (N*func(ri) - x_val[0])**2 / (2*x_val[0])

    return np.sum(func_array)


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
    #in_array = read_file("design_1p0_no_outer_ge_2pix_smear.npy")
    #in_array = read_file("design_0p1_0p08_no_outer_ge_2pix_smear.npy")
    #in_array = read_file("design_0p1_neg0p08_no_outer_ge_2pix_smear.npy")
    #in_array = read_file("design_z10_no_outer_ge_2pix.npy")
    #in_array = read_file("design_zneg10_no_outer_ge_2pix.npy")
    #in_array = read_file("design_0_nosmear_no_outer_ge_2pix.npy")
    #in_array = read_file("design_0_no_outer_ge_2pix_nosmear_phifix.npy")
    #in_array = read_file("design_z10_0_no_outer_ge_2pix_nosmear_phifix.npy")
    #in_array = read_file("design_0p1_no_outer_ge_2pix_nosmear_phifix.npy")
    #in_array = read_file("design_0_no_outer_ge_2pix_nosmear_phifix.npy")
    #in_array = read_file("design_0p1_no_outer_ge_2pix_smear_phifix.npy")
    #in_array = read_file("design_0_no_outer_all_pix_nosmear_phifix.npy")
    #in_array = read_file("design_0p1_no_outer_all_pix_nosmear_phifix.npy")
    #in_array = read_file("design_0p1_II_no_outer_all_pix_nosmear.npy")
    in_array = read_file("design_0p2_no_outer_all_pix_nosmear.npy")
    #in_array = read_file("design_0p3_no_outer_all_pix_nosmear.npy")
    #in_array = read_file("design_0p1_no_outer_all_pix_nosmear_charge.npy")
    #in_array = read_file("design_0_no_outer_all_pix_nosmear_chargel150_g25.npy")
    #in_array = read_file("design_0_no_outer_all_pix_nosmear_chargel200_g25.npy")
    #in_array = read_file("design_0_no_outer_all_pix_nosmear_chargel200.npy")
    #in_array = read_file("design_0p3_no_outer_all_pix_nosmear_chargel200.npy")
    #in_array = read_file("design_0p1_no_outer_all_pix_nosmear_chargel200.npy")
    #in_array = read_file("design_0p1_no_outer_all_pix_nosmear_charge.npy")
    #in_array = read_file("design_0p1_no_outer_all_pix_nosmear_chargel50.npy")
    #in_array = read_file("design_0p2_no_outer_all_pix_nosmear_chargel50.npy")
    #in_array = read_file("design_0p3_no_outer_all_pix_nosmear_chargel50.npy")
    #in_array = read_file("design_0p3_no_outer_all_pix_nosmear_chargel200.npy")
    remake_arrays(in_array)