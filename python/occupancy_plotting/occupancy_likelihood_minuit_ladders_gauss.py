import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import iminuit as im
import string

alpha_low = string.ascii_lowercase


def read_file(input_file_):
    return np.load(input_file_, allow_pickle=True, encoding='latin1')


def remake_arrays(input_arr_):
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

    roc_index = range(0, n_z_bins*n_phi_bins)

    fig = plt.figure()
    axs = fig.add_subplot(111, projection='3d')

    occ = []
    x_array = []
    y_array = []
    z_array = []
    phi_array = []
    r_array = []

    # section off rocs into roc ladders (12 'ladders'), each true ladder is split in half 6 * 2 = 12
    for roc in roc_index:
        if roc%12==0: continue
        if roc%12==1: continue
        if roc%12==2: continue
        if roc%12==3: continue
        if roc%12==4: continue
        if roc%12==5: continue
        #if roc%12==6: continue
        #if roc%12==7: continue
        if roc%12==8: continue
        if roc%12==9: continue
        if roc%12==10: continue
        if roc%12==11: continue

        #if 31 < int((roc/12)%64) < 40 and roc%12==0: continue
        #if 31 < int((roc/12)%64) < 40 and roc%12==1: continue

        #if 0 <= int((roc/12)%64) < 8: continue
        #if 55 < int((roc/12)%64) < 64: continue
        ###if 0 <= int((roc/12)%64) < 4: continue
        ###if 55 < int((roc/12)%64) < 60: continue


        #if int((roc/12)%64) == 9 and roc%12==0: continue
        #if int((roc/12)%64) == 16 and roc%12==0: continue
        #if int((roc/12)%64) == 26 and roc%12==0: continue

        #if int((roc/12)%64) == 16 and roc%12==1: continue
        #if int((roc/12)%64) == 40 and roc%12==1: continue
        #if int((roc/12)%64) == 8 and -1 < roc%12 < 2: continue

        #if int((roc/12)%64) == 41 and roc%12==2: continue
        #if int((roc/12)%64) == 16 and roc%12==2: continue

        #if 19 <int((roc/12)%64) < 24 and roc%12==4: continue
        #if int((roc/12)%64) == 8 and roc%12==4: continue
        #if int((roc/12)%64) == 21 and roc%12==4: continue
        #if int((roc/12)%64) == 27 and roc%12==4: continue
        #if int((roc/12)%64) == 28 and roc%12==4: continue
        #if int((roc/12)%64) == 31 and roc%12==4: continue
        #if int((roc/12)%64) == 32 and roc%12==4: continue
        #if int((roc/12)%64) == 35 and roc%12==4: continue

        #if int((roc/12)%64) == 8 and roc%12==5: continue
        #if int((roc/12)%64) == 16 and roc%12==5: continue
        #if int((roc/12)%64) == 20 and roc%12==5: continue
        #if int((roc/12)%64) == 21 and roc%12==5: continue
        #if int((roc/12)%64) == 23 and roc%12==5: continue
        #if int((roc/12)%64) == 28 and roc%12==5: continue
        #if int((roc/12)%64) == 29 and roc%12==5: continue
        #if int((roc/12)%64) == 41 and roc%12==5: continue
        #if int((roc/12)%64) == 40 and roc%12==5: continue

        #if int((roc/12)%64) == 38 and 3 < roc%12 < 6: continue

        #if int((roc/12)%64) == 8 and roc%12==6: continue
        #if int((roc/12)%64) == 23 and roc%12==6: continue
        #if int((roc/12)%64) == 41 and roc%12==6: continue
        #if 35 < int((roc/12)%64) < 40 and roc%12==6: continue
        #if 24 < int((roc/12)%64) < 32 and roc%12==6: continue

        #if int((roc/12)%64) == 8 and roc%12==7: continue
        #if int((roc/12)%64) == 21 and roc%12==7: continue
        #if int((roc/12)%64) == 23 and roc%12==7: continue
        #if 24 < int((roc/12)%64) < 32 and roc%12==7: continue

        #if int((roc/12)%64) == 8 and roc%12==8: continue
        #if int((roc/12)%64) == 14 and roc%12==8: continue
        #if int((roc/12)%64) == 13 and roc%12==8: continue
        #if int((roc/12)%64) == 23 and roc%12==8: continue
        #if int((roc/12)%64) == 50 and roc%12==8: continue

        #if int((roc/12)%64) == 8 and roc%12==9: continue
        #if int((roc/12)%64) == 21 and roc%12==9: continue
        #if int((roc/12)%64) == 22 and roc%12==9: continue
        #if int((roc/12)%64) == 23 and roc%12==9: continue
        #if int((roc/12)%64) == 38 and roc%12==9: continue

        #if int((roc/12)%64) == 11 and roc%12==10: continue
        #if int((roc/12)%64) == 21 and roc%12==10: continue
        #if int((roc/12)%64) == 25 and roc%12==10: continue
        #if int((roc/12)%64) == 28 and roc%12==10: continue

        #if int((roc/12)%64) == 26 and roc%12==11: continue
        #if int((roc/12)%64) == 35 and roc%12==11: continue
        #if int((roc/12)%64) == 32 and roc%12==11: continue
        #if int((roc/12)%64) == 47 and roc%12==11: continue
        #if int((roc/12)%64) == 31 and 9 < roc%12 < 12: continue
        #if int((roc/12)%64) == 8 and 9 < roc%12 < 12: continue

        #if int((roc/12)%64) == 24 and roc%12==0: continue
        #if int((roc/12)%64) == 24 and roc%12==1: continue
        #if int((roc/12)%64) == 24 and roc%12==4: continue
        #if int((roc/12)%64) == 24 and roc%12==6: continue
        #if int((roc/12)%64) == 24 and roc%12==7: continue

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

        #if occ_tmp.size == 0: continue
        occ.append(np.sum(occ_tmp))
        if np.sum(occ_tmp)==0:
            x_array.append(np.average(x))
            y_array.append(np.average(y))
            z_array.append(np.average(z))
            phi_array.append(np.average(phi))
            r_array.append(np.average(r))
        else:
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

    # removing rocs
    remove_z = (z_array > -25) * (z_array < 25)
    #remove_z = (z_array > -20) * (z_array < 20)

    remove_blips = (z_array < -21) + (z_array > -20)
    remove_blips *= (z_array < -14.5) + (z_array > -13.5)
    remove_blips *= (z_array < -7.5) + (z_array > -6.5)
    remove_blips *= (z_array < -1) + (z_array > 0)
    remove_blips *= (z_array < 5.75) + (z_array > 6.5)
    remove_blips *= (z_array < 12.5) + (z_array > 13.5)
    remove_blips *= (z_array < 19) + (z_array > 20)

    occ = occ[remove_z*remove_blips]
    x_array = x_array[remove_z*remove_blips]
    y_array = y_array[remove_z*remove_blips]
    z_array = z_array[remove_z*remove_blips]
    phi_array = phi_array[remove_z*remove_blips]

    #occ = occ[remove_z]
    #x_array = x_array[remove_z]
    #y_array = y_array[remove_z]
    #z_array = z_array[remove_z]
    #phi_array = phi_array[remove_z]
    #def nll(x0, y0, z0, n, b1, b2, b3, a1, a3, c1, c3, ga1, ga2, ga3, gc1, gc3):
    def nll(x0, y0, z0, n, b1, b2, b3, a1, a3, c1, c3, ga1, ga3, gc1, gc3):

        ri = np.float64(np.sqrt((x_array - x0) ** 2 + (y_array - y0) ** 2 + (z_array - z0) ** 2))
        phi_cor = np.arctan2(y_array-y0, x_array-x0)
        a = np.float64(b_par(phi_cor, a1, b2, a3))
        b = np.float64(b_par(phi_cor, b1, b2, b3))
        c = np.float64(b_par(phi_cor, c1, b2, c3))
        #ga = np.float64(b_par(phi_cor, ga1, ga2, ga3))
        #gc = np.float64(b_par(phi_cor, gc1, ga2, gc3))
        ga = np.float64(b_par(phi_cor, ga1, b2-np.pi, ga3))
        gc = np.float64(b_par(phi_cor, gc1, b2-np.pi, gc3))

        # gb = np.sqrt(x0**2 + y0**2 + z0**2)
        gb = z0
        func_array = -np.log(1 / (np.sqrt(2 * np.pi * occ))) + (n * (func(ri, a, b, c)
                                                                     + gauss(z_array, ga, gb, gc))
                                                                - occ) ** 2 / (2 * np.pi * occ)

        return np.sum(func_array)

    axs.plot(z_array, phi_array, occ, 'b*')
    minuit = im.Minuit(nll, x0=0., y0=0., z0=0., n=1,
                       #a1=-6400, a3=89600,
                       #b1=-0.01, b2=0.94, b3=0.4,
                       #c1=-70, c3=1360,
                       #ga1=-1000, ga3=-1200,
                       #ga2=0.1,
                       #gc1=-0.019, gc3=2.153,
                       # (0, 0) fit
                       a1=-0.19e4, a3=0.86e5,
                       b1=0, b2=0, b3=1.212,
                       c1=-110, c3=1260,
                       ga1=-500, ga3=-1.2e4,
                       gc1=-0.031, gc3=2.129,
                       # (0, 0) no smear fit
                       #a1=-0.19e4, a3=0.86e5,
                       #b1=0, b2=0, b3=1.212,
                       #c1=33, c3=1892,
                       #ga1=590, ga3=-2340,
                       # ga2=0,
                       #gc1=0.1, gc3=1.94,

                       # for PU fit
                       #a1=-0.345e5, a3=2.644e5,
                       #b1=0, b2=0, b3=0.322,
                       #c1=1.06e4, c3=-0.576e5,
                       #ga1=-1000, ga3=-1200,
                       #gc1=-0.031, gc3=2.129,
                       error_x0=0.001, error_y0=0.001, error_z0=0.001, error_n=0.01,
                       error_b1=0.01, error_b2=0.01, error_b3=0.01,
                       error_a1=1, error_a3=0.1,
                       error_c1=1, error_c3=0.1,
                       #error_ga=0.1,
                       #error_gb=0.1,
                       #error_gc=0.1,
                       error_ga1=0.1, error_ga3=0.1,
                       #error_ga2=0.01,
                       error_gc1=0.1, error_gc3=0.1,
                       #error_gb1=0.1,
                       #error_gb3=0.1,
                       # error_a1c=0.1, error_c1=0.1,
                       fix_n=True,
                       #fix_gb=True,
                       #fix_x0=True, fix_y0=True,
                       #fix_z0=True,
                       #fix_b1=True, fix_b2=True, fix_b3=True,
                       #fix_a1=True, fix_a3=True,
                       #fix_c1=True, fix_c3=True,
                       #fix_ga1=True, fix_ga3=True,
                       #fix_gc1=True, fix_gc3=True,
                       #fix_z0=True,
                       limit_x0=(-5, 5), limit_y0=(-1, 1), limit_z0=(-10, 10),
                       #limit_b1=(0.2, None),
                       limit_b2=(-np.pi, np.pi),
                       limit_b3=(0, 3),
                       #limit_ga1=(0.01, None), limit_ga3=(0.01, None),
                       #limit_ga2=(-np.pi, np.pi),
                       #limit_ga=(0, None),
                       #limit_gc1=(-10, 10), limit_gc3=(-10, 10),
                       # limit_a1a=(0.01, None), limit_a1b=(0.01, None), limit_a1c=(0.01, None),
                       # limit_a3a=(0.01, None), limit_a3b=(0.01, None), limit_a3c=(0.01, None),
                       # limit_a1=(0., 1e6),
                       # limit_c1=(0, 1e4),
                       errordef=1)

    minuit.migrad()
    minuit.hesse()
    minuit.minos()
    print(minuit.get_param_states())
    print(minuit.get_fmin())

    # axs.legend()

    #plt.show()


def b_par(x, b1=0.0, b2=0.0, b3=1.25e6):
    return b1*np.sin(x-b2)+b3


def func(x, a, b, c):

    return a * (1 / x ** b) + c


def gauss(x, ga, gb, gc):
    return ga * np.exp(-0.5 * (x-gb)**2 / gc**2)


if __name__ == "__main__":
    in_array = read_file("design_0_no_outer_all_pix_smear.npy")
    #in_array = read_file("design_0p1_no_outer_all_pix_smear.npy")
    #in_array = read_file("design_0p1_0p08_no_outer_all_pix_smear.npy")
    #in_array = read_file("design_neg0p1_0p08_no_outer_all_pix_smear.npy")
    #in_array = read_file("design_0p2_no_outer_all_pix_smear.npy")
    #in_array = read_file("design_0p3_no_outer_all_pix_smear.npy")
    #in_array = read_file("design_0p1_z10_no_outer_all_pix_smear.npy")
    #in_array = read_file("design_0p2_no_outer_all_pix_PU.npy")
    #in_array = read_file("design_0p3_no_outer_all_pix_PU.npy")

    #in_array = read_file("singlemu_no_outer_all_pix.npy")
    #in_array = read_file("singlemu_2018AB_no_outer_all_pix.npy")
    #in_array = read_file("singlemu_2018C2_no_outer_all_pix.npy")

    #in_array = read_file("design_0_no_outer_all_pix_nosmear_phifix.npy")
    #in_array = read_file("design_0p01_no_outer_all_pix_nosmear.npy")
    #in_array = read_file("design_0p03_no_outer_all_pix_nosmear.npy")
    #in_array = read_file("design_neg0p08_no_outer_all_pix_nosmear.npy")
    #in_array = read_file("design_0p1_no_outer_all_pix_nosmear_phifix.npy")
    #in_array = read_file("design_0p1_II_no_outer_all_pix_nosmear.npy")
    #in_array = read_file("design_0p2_no_outer_all_pix_nosmear.npy")
    #in_array = read_file("design_0p3_no_outer_all_pix_nosmear.npy")

    #in_array = read_file("design_0p1_no_outer_all_pix_nosmear_50k.npy")
    #in_array = read_file("design_0p1_II_no_outer_all_pix_nosmear_50k.npy")
    #in_array = read_file("design_0p2_no_outer_all_pix_nosmear_50k.npy")
    #in_array = read_file("design_0p3_no_outer_all_pix_nosmear_50k.npy")

    #in_array = read_file("design_0p1_no_outer_all_pix_nosmear_25k.npy")
    #in_array = read_file("design_0p1_II_no_outer_all_pix_nosmear_25k.npy")
    #in_array = read_file("design_0p2_no_outer_all_pix_nosmear_25k.npy")
    #in_array = read_file("design_0p3_no_outer_all_pix_nosmear_25k.npy")

    #in_array = read_file("design_0p1_no_outer_all_pix_smear_50k.npy")
    #in_array = read_file("design_0p2_no_outer_all_pix_smear_50k.npy")
    #in_array = read_file("design_0p3_no_outer_all_pix_smear_50k.npy")

    #in_array = read_file("design_0p1_no_outer_all_pix_smear_25k.npy")
    #in_array = read_file("design_0p2_no_outer_all_pix_smear_25k.npy")
    #in_array = read_file("design_0p3_no_outer_all_pix_smear_25k.npy")
    remake_arrays(in_array)