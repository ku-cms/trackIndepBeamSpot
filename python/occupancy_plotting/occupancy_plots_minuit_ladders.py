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

    # need phi-binning corresponding to 1 roc (maybe 2?)
    w_phi_bins = 80

    n_z_bins = int(3328 / w_z_bins)  # 3328 is number of pixels in a ladder row

    n_phi_bins = int(960 / w_phi_bins)  # 1440 is number of pixels around phi for all ladders, 960 for inner ladders
    #n_phi_bins = int(1440 / w_phi_bins)  # 1440 is number of pixels around phi for all ladders, 960 for inner ladders

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
    #fig, axs = plt.subplots(12, 2, sharex=True, sharey=True, figsize=(20, 20), tight_layout=False) # fraction of rocs and modules
    fig, axs = plt.subplots(12, sharex=True, sharey=True, figsize=(20, 20), tight_layout=False) # fraction of rocs and modules
    #fig, axs = plt.subplots(3, sharex=True, sharey=True, figsize=(20, 20), tight_layout=True) # fraction of rocs and modules

    # minus - 0-383
    # plus - 384-767
    occ_plus = []
    occ_minus = []
    r_sph_plus = []
    r_sph_minus = []
    phi_plus = []
    phi_minus = []
    z_plus = []
    z_minus = []

    # section off rocs into roc ladders (12 'ladders'), each true ladder is split in half 6 * 2 = 12
    n_ladders = 12

    for x in range(n_ladders):
        occ_plus.append([])
        occ_minus.append([])
        r_sph_plus.append([])
        r_sph_minus.append([])
        phi_plus.append([])
        phi_minus.append([])
        z_plus.append([])
        z_minus.append([])
    i_ladder = 0
    for roc in roc_index:

        occ_tmp = np.concatenate(array_by_rocs[roc, :, :, 0])
        r = np.concatenate(array_by_rocs[roc, :, :, 1])
        phi = np.concatenate(array_by_rocs[roc, :, :, 2])
        z = np.concatenate(array_by_rocs[roc, :, :, 3])
        z_avg = np.nanmean(z)

        r_sph_tmp = np.sqrt(r**2+z**2)

        r_sph_minus[i_ladder].append(np.nanmean(r_sph_tmp))
        occ_minus[i_ladder].append(np.sum(occ_tmp[~np.isnan(occ_tmp)]))
        phi_minus[i_ladder].append(np.nanmean(phi))
        z_minus[i_ladder].append(np.nanmean(z))
        #if roc < 384:
        #    r_sph_minus[i_ladder].append(np.nanmean(r_sph_tmp))
        #    occ_minus[i_ladder].append(np.sum(occ_tmp[~np.isnan(occ_tmp)]))
        #    phi_minus[i_ladder].append(np.nanmean(phi))
        #    z_minus[i_ladder].append(np.nanmean(z))
        #else:
        #    r_sph_plus[i_ladder].append(np.nanmean(r_sph_tmp))
        #    occ_plus[i_ladder].append(np.sum(occ_tmp[~np.isnan(occ_tmp)]))
        #    phi_plus[i_ladder].append(np.nanmean(phi))
        #    z_plus[i_ladder].append(np.nanmean(z))

        i_ladder += 1
        if i_ladder == n_ladders:
            i_ladder = 0

    #avg_phi_plus = np.array([np.mean(iphi) for iphi in phi_plus])
    avg_phi_minus = np.array([np.mean(iphi) for iphi in phi_minus])
    #phi_sort_plus = np.argsort(avg_phi_plus)
    phi_sort_minus = np.argsort(avg_phi_minus)
    # sort ladders by phi
    #r_sph_plus_tmp = np.array(r_sph_plus)
    r_sph_minus_tmp = np.array(r_sph_minus)
    #occ_plus_tmp = np.array(occ_plus)
    occ_minus_tmp = np.array(occ_minus)
    #z_plus_tmp = np.array(z_plus)
    z_minus_tmp = np.array(z_minus)
    #r_sph_plus_tmp = r_sph_plus_tmp[phi_sort_plus]
    r_sph_minus_tmp = r_sph_minus_tmp[phi_sort_minus]
    #occ_plus_tmp = occ_plus_tmp[phi_sort_plus]
    occ_minus_tmp = occ_minus_tmp[phi_sort_minus]
    #z_plus_tmp = z_plus_tmp[phi_sort_plus]
    z_minus_tmp = z_minus_tmp[phi_sort_minus]

    #print(avg_phi_plus[phi_sort_plus])
    print(avg_phi_minus[phi_sort_minus])
    #n_ladders = 6  # for running over modules in a ring instead of rocs
    for x in range(n_ladders):
        #fig, axs = plt.subplots(3, sharex=True, sharey=True, figsize=(20, 20),
        #                        tight_layout=True)  # fraction of rocs and modules
        r_sph_minus[x] = np.array(r_sph_minus_tmp[x])
        occ_minus[x] = np.array(occ_minus_tmp[x])
        z_minus[x] = np.array(z_minus_tmp[x])
        sort_minus = np.argsort(r_sph_minus[x])
        r_sph_minus[x] = r_sph_minus[x][sort_minus]
        occ_minus[x] = occ_minus[x][sort_minus]

        #r_sph_plus[x] = np.array(r_sph_plus_tmp[x])
        #occ_plus[x] = np.array(occ_plus_tmp[x])
        #z_plus[x] = np.array(z_plus_tmp[x])
        #sort_plus = np.argsort(r_sph_plus[x])
        #r_sph_plus[x] = r_sph_plus[x][sort_plus]
        #occ_plus[x] = occ_plus[x][sort_plus]
        #z_plus[x] = z_plus[x][sort_plus]

        # removing rocs
        #remove_z_plus = (z_plus[x] < 25)
        #remove_z_minus = (z_minus[x] > -25)
        #remove_z_plus = (r_sph_plus[x] <= 16)
        #remove_z_minus = (r_sph_minus[x] <= 16)
        remove_z_minus = (r_sph_minus[x] > 12)
        #remove_z_minus = (r_sph_minus[x] <= 12)
        #remove_z_minus = (z_minus[x] > -25) * (z_minus[x] < 25)
        #remove_z_minus = (z_minus[x] >= -16) * (z_minus[x] <= 16)
        #remove_z_minus = (z_minus[x] < -16) + (z_minus[x] > 16)
        # remove_z = ((z_plus[x] > 11.5) + z_plus[x] < 8.5)) * (z_plus[x] > -25) * (z_plus[x] < 25)
        # remove_z = ((z_plus[x] > -8.5) + (z_plus[x] < -11.5)) * (z_plus[x] > -25) * (z_plus[x] < 25)
        remove_blips_minus = (z_minus[x] < 5.75) + (z_minus[x] > 6.5)
        remove_blips_minus *= (z_minus[x] < 12.5) + (z_minus[x] > 13.5)
        remove_blips_minus *= (z_minus[x] < 19) + (z_minus[x] > 20)
        remove_blips_minus *= (z_minus[x] < -21) + (z_minus[x] > -20)
        remove_blips_minus *= (z_minus[x] < -14.5) + (z_minus[x] > -13.5)
        remove_blips_minus *= (z_minus[x] < -7.5) + (z_minus[x] > -6.5)
        remove_blips_minus *= (z_minus[x] < -25) + (z_minus[x] > 25)

        #r_sph_plus[x] = r_sph_plus[x][remove_z_plus*remove_blips_plus]
        r_sph_minus[x] = r_sph_minus[x][remove_z_minus*remove_blips_minus]
        #r_sph_minus[x] = r_sph_minus[x][remove_z_minus]
        #occ_plus[x] = occ_plus[x][remove_z_plus*remove_blips_plus]
        occ_minus[x] = occ_minus[x][remove_z_minus*remove_blips_minus]
        #occ_minus[x] = occ_minus[x][remove_z_minus]

        z_condense = []
        r_sph_condense = []
        occ_condense = []
        for ir, r in enumerate(r_sph_minus[x]):
            if ir % 2 == 0 and ir+1 < len(r_sph_minus[x]):
                r_sph_condense.append((r + r_sph_minus[x][ir + 1]) / 2)
                occ_condense.append(occ_minus[x][ir] + occ_minus[x][ir + 1])
                z_condense.append((z_minus[x][ir] + z_minus[x][ir + 1]) / 2)
            else:
                continue

        r_sph_minus[x] = np.array(r_sph_condense)
        occ_minus[x] = np.array(occ_condense)
        z_minus[x] = np.array(z_condense)

        axs[x].plot(r_sph_minus[x], occ_minus[x], 'b*', label='full z hl '+str(x))
        #axs[x, 0].plot(r_sph_minus[x], occ_minus[x], 'b*', label='minus z - hl '+str(x))
        #axs[x, 1].plot(r_sph_plus[x], occ_plus[x], 'b*', label='plus z - hl '+str(x))
        #axs[0].plot(r_sph_minus[x], occ_minus[x], 'b*', label='minus z - hl '+str(x))
        #axs[1].plot(r_sph_plus[x], occ_plus[x], 'b*', label='plus z - hl '+str(x))

        #axs[2].plot(r_sph_minus[x], occ_minus[x], 'r*', label='occ minus')
        #axs[2].plot(r_sph_plus[x], occ_plus[x], 'b*', label='occ plus')
        chi2_minus = pf.Chi2Regression(func, r_sph_minus[x], occ_minus[x])
        #chi2_plus = pf.Chi2Regression(func, r_sph_plus[x], occ_plus[x])
        #minuit_minus = im.Minuit(chi2_minus, a=0, b=1.16, c=2000,
        #                         error_a=1, error_b=0.01, error_c=1,
        #                         fix_b=False,
        #                         limit_a=(None, None), limit_b=(None, None), limit_c=(None, None),
        #                         errordef=1)
        minuit_minus = im.Minuit(chi2_minus, a=200000, b=-0.4, c=2000,
                                error_a=1, error_b=0.01, error_c=1,
                                fix_b=False,
                                limit_a=(None, None), limit_b=(0, 10), limit_c=(None, None),
                                errordef=1)
        #minuit_minus = im.Minuit(chi2_minus, a=200000, b=0.5989, c=2000, r0=0,
        #                         error_a=1, error_b=0.01, error_c=1, error_r0=0.01,
        #                         fix_b=True,
        #                         limit_a=(None, None), limit_b=(0, 10), limit_c=(None, None), limit_r0=(None, None),
        #                         errordef=1)
        #minuit_plus = im.Minuit(chi2_plus, a=200000, b=0.5989, c=2000, r0=0,
        #                        error_a=1, error_b=0.01, error_c=1, error_r0=0.01,
        #                        fix_b=True,
        #                        limit_a=(None, None), limit_b=(0, 10), limit_c=(None, None), limit_r0=(None, None),
        #                        errordef=1)
        #print(minuit.get_param_states())
        minuit_minus.migrad()
        #minuit_plus.migrad()
        #print(minuit_plus.get_param_states())
        #print(minuit_plus.get_fmin())
        print(minuit_minus.get_param_states())
        print(minuit_minus.get_fmin())

        # chi2_minus.draw(minuit_minus, axs[x, 0])
        # chi2_plus.draw(minuit_plus, axs[x, 1])
        # chi2_minus.draw(minuit_minus, axs[0])
        # chi2_plus.draw(minuit_plus, axs[1])
        param_string = ''
        for p in minuit_minus.parameters:
            param_string += '{}: {} +/- {}\n'.format(p, np.format_float_scientific(minuit_minus.values[p], precision=3),
                                                     np.format_float_scientific(minuit_minus.errors[p], precision=1))
        #axs[x, 0].text(15, 550000, param_string, fontsize='xx-small')
        #axs[x].text(7.5, 20000, param_string, fontsize='xx-small')
        axs[x].text(22, 6500, param_string, fontsize='xx-small')

        #param_string = ''
        #for p in minuit_plus.parameters:
        #    param_string += '{}: {} +/- {}\n'.format(p, np.format_float_scientific(minuit_plus.values[p], precision=3),
        #                                             np.format_float_scientific(minuit_plus.errors[p], precision=1))
        #axs[x, 1].text(15, 550000, param_string, fontsize='xx-small')

        axs[x].plot(r_sph_minus[x], func(r_sph_minus[x], *minuit_minus.values.values()), color='red', label='fit')
        #axs[x, 0].plot(r_sph_minus[x], func(r_sph_minus[x], *minuit_minus.values.values()), color='red', label='fit minus')
        #axs[x, 1].plot(r_sph_plus[x], func(r_sph_plus[x], *minuit_plus.values.values()), color='red', linestyle='solid', label='fit plus')
        #axs[0].plot(r_sph_minus[x], func(r_sph_minus[x], *minuit_minus.values.values()), color='red', label='fit minus')
        #axs[1].plot(r_sph_plus[x], func(r_sph_plus[x], *minuit_plus.values.values()), color='red', linestyle='solid', label='fit plus')
        #axs[2].plot(r_sph_minus[x], func(r_sph_minus[x], *minuit_minus.values.values()), color='red', label='fit minus')
        #axs[2].plot(r_sph_plus[x], func(r_sph_plus[x], *minuit_plus.values.values()), color='blue', linestyle='dashed', label='fit plus')

        axs[x].legend(fontsize='xx-small', loc='upper right')
        #axs[x, 1].legend(fontsize='xx-small', loc='upper right')
        #axs[x, 0].legend(fontsize='xx-small', loc='upper right')
        #axs[2].legend()
        #axs[1].legend(fontsize='xx-small')
        #axs[0].legend(fontsize='xx-small')
        #axs[2].grid()

        #plt.savefig('plots_0_smear/fit_ladders_pm_'+str(x)+'_layer1.png')
        #plt.clf()
    plt.show()


#def func(x, a, b, c, r0):
#    return a*(1/(x-r0)**(2*b)) + c


def func(x, a, b, c):
#def func(x, a, b):

    return a*(1/(x)**(b)) + c
    #return a*x + b


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
    #in_array = read_file("design_0_no_outer_all_pix_nosmear_phifix.npy")
    in_array = read_file("design_0p3_no_outer_all_pix_nosmear.npy")
    #in_array = read_file("design_0p1_II_no_outer_all_pix_nosmear.npy")
    #in_array = read_file("design_0_no_outer_all_pix_nosmear_chargel200.npy")
    #in_array = read_file("design_0p3_no_outer_all_pix_nosmear_chargel200.npy")
    #in_array = read_file("design_0p1_no_outer_all_pix_nosmear_phifix.npy")
    remake_arrays(in_array)