import numpy as np
import matplotlib.pyplot as plt
import iminuit as im
import string

alpha_low = string.ascii_lowercase


def read_file(input_file_):
    return np.load(input_file_, allow_pickle=True, encoding='latin1')


def remake_arrays(input_arr_, plot_dir, plot_name):
    useWeightedAve = False
    
    if plot_dir[-1] != "/":
        plot_dir += "/"

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

        occ_tmp = np.concatenate(array_by_rocs[roc, :, :, 0])
        r       = np.concatenate(array_by_rocs[roc, :, :, 1])
        phi     = np.concatenate(array_by_rocs[roc, :, :, 2])
        z       = np.concatenate(array_by_rocs[roc, :, :, 3])
        #print("z: {0}".format(z))
        #print("len(z): {0}".format(len(z)))
        #z_avg = np.nanmean(z)

        x = r[~np.isnan(z)] * np.cos(phi[~np.isnan(z)])
        y = r[~np.isnan(z)] * np.sin(phi[~np.isnan(z)])
        r = r[~np.isnan(z)]
        phi = phi[~np.isnan(z)]
        occ_tmp = occ_tmp[~np.isnan(z)]
        z = z[~np.isnan(z)]

        
        if useWeightedAve:
            occ.append(np.sum(occ_tmp))
            x_array.append(np.average(x,        weights=occ_tmp))
            y_array.append(np.average(y,        weights=occ_tmp))
            z_array.append(np.average(z,        weights=occ_tmp))
            phi_array.append(np.average(phi,    weights=occ_tmp))
            r_array.append(np.average(r,        weights=occ_tmp))
        else:
            occ.append(np.sum(occ_tmp))
            x_array.append(np.average(x))
            y_array.append(np.average(y))
            z_array.append(np.average(z))
            phi_array.append(np.average(phi))
            r_array.append(np.average(r))

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

    remove_blips =  (z_array < -21)   + (z_array > -20)
    remove_blips *= (z_array < -14.5) + (z_array > -13.5)
    remove_blips *= (z_array < -7.5)  + (z_array > -6.5)
    remove_blips *= (z_array < 5.75)  + (z_array > 6.5)
    remove_blips *= (z_array < 12.5)  + (z_array > 13.5)
    remove_blips *= (z_array < 19)    + (z_array > 20)

    occ         = occ[remove_z*remove_blips]
    x_array     = x_array[remove_z*remove_blips]
    y_array     = y_array[remove_z*remove_blips]
    z_array     = z_array[remove_z*remove_blips]
    phi_array   = phi_array[remove_z*remove_blips]

    def nll(x0, y0, z0, n, b1, b2, b3, a1, a3, c1, c3):
        ri = np.float64(np.sqrt((x_array - x0) ** 2 + (y_array - y0) ** 2 + (z_array - z0) ** 2))
        phi_cor = np.arctan2(y_array-y0, x_array-x0)
        a = np.float64(b_par(phi_cor, a1, b2, a3))
        b = np.float64(b_par(phi_cor, b1, b2, b3))
        c = np.float64(b_par(phi_cor, c1, b2, c3))

        func_array = -np.log(1 / (np.sqrt(2 * np.pi * occ))) + (n * func(ri, a, b, c) - occ) ** 2 / (2 * np.pi * occ)

        return np.sum(func_array)

    minuit = im.Minuit(nll, x0=0, y0=0, z0=-0.019, n=1,
                       a1=4e4, a3=1.26e5,
                       b1=0.0, b2=3, b3=1.167,
                       c1=1e3, c3=3.12e3,
                       #error_x0=0.001, error_y0=0.001, error_z0=0.001, error_n=0.01,
                       #error_b1=0.01, error_b2=0.01, error_b3=0.01,
                       #error_a1=10, error_a3=0.1,
                       #error_c1=10, error_c3=0.1,
                       #error_a1c=0.1, error_c1=0.1,
                       #fix_n=True,
                       #fix_b3=True,
                       #fix_b2=True,
                       #fix_a3=True, fix_c3=True,
                       #fix_a1=False, fix_c1=False,
                       #fix_z0=True,
                       #limit_x0=(-1, 1), limit_y0=(-1, 1), limit_z0=(-1, 1),
                       #limit_b1=(0, 2),
                       #limit_b2=(-np.pi, np.pi),
                       #limit_b2=(-1, 1),
                       #limit_b3=(1, 3),
                       #limit_a1a=(0.01, None), limit_a1b=(0.01, None), limit_a1c=(0.01, None),
                       #limit_a3a=(0.01, None), limit_a3b=(0.01, None), limit_a3c=(0.01, None),
                       #limit_a1=(0., 1e6),
                       #limit_c1=(0, 1e4),
                       #errordef=1)
                      )

    minuit.migrad()
    #print(minuit.get_param_states())
    #print(minuit.get_fmin())
    
    # plotting
    print("z_array = {0}".format(z_array))
    print("phi_array = {0}".format(phi_array))
    print("occ = {0}".format(occ))
    #print("f = {0}".format(func(np.sqrt(z_array ** 2 + phi_array ** 2), a=1, b=1, c=1)))

    x0 = 0.0
    y0 = 0.0
    z0 = 0.0
    my_r = np.float64(np.sqrt((x_array - x0) ** 2 + (y_array - y0) ** 2 + (z_array - z0) ** 2))
    print("my_r = {0}".format(my_r))
    
    z_2d, phi_2d = np.meshgrid(z_array, phi_array)
    output_2d = 50000 * func(np.sqrt(z_2d ** 2 + phi_2d ** 2), a=1, b=1, c=1)
    
    #axs.contour3D(z_2d, phi_2d, output_2d, 50, cmap='binary')
    axs.plot_wireframe(z_2d, phi_2d, output_2d, color='black')
    #axs.plot_surface(z_2d, phi_2d, output_2d, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
    axs.plot(z_array, phi_array, occ, 'b*')

    # labels
    axs.set_title("Pixel Occupancy", fontsize=20)
    axs.set_xlabel("z",              fontsize=16)
    axs.set_ylabel(r"$\phi$",        fontsize=16)

    # axs.legend()
    #plt.show()
    
    plt.savefig(plot_dir + plot_name + '.png', bbox_inches='tight')
    
    # select view angles for view_init(elev, azim)
    angles = []
    angles.append([0,     0])
    angles.append([270,   0])
    angles.append([0,   270])

    for angle in angles:
        full_name = "{0}{1}_{2}_{3}.png".format(plot_dir, plot_name, angle[0], angle[1])
        # select viewing angle
        axs.view_init(elev=angle[0], azim=angle[1])
        plt.savefig(full_name, bbox_inches='tight')
    

def b_par(x, b1=0.0, b2=0.0, b3=1.25e6):
    return b1*np.sin(x-b2)+b3


def func(x, a, b, c):

    return a * (1 / x ** b) + c


if __name__ == "__main__":
    
    data_dir  = "data/"
    plot_dir  = "plots/"
    
    #in_array = read_file(data_dir + "design_0_no_outer_all_pix_nosmear_phifix.npy")
    #in_array = read_file(data_dir + "design_0p01_no_outer_all_pix_nosmear.npy")
    #in_array = read_file(data_dir + "design_0p03_no_outer_all_pix_nosmear.npy")
    #in_array = read_file(data_dir + "design_neg0p08_no_outer_all_pix_nosmear.npy")
    #in_array = read_file(data_dir + "design_0p1_no_outer_all_pix_nosmear_phifix.npy")
    #in_array = read_file(data_dir + "design_0p1_II_no_outer_all_pix_nosmear.npy")
    #in_array = read_file(data_dir + "design_0p2_no_outer_all_pix_nosmear.npy")
    #in_array = read_file(data_dir + "design_0p3_no_outer_all_pix_nosmear.npy")
    
    #in_array  = read_file(data_dir + "design_0p1_no_outer_all_pix_smear_charge_l1000_size_1_50_v3.npy")
    #plot_name = "fig_v3"
    
    #in_array  = read_file(data_dir + "design_0p1_no_outer_all_pix_smear_charge_l1000_size_1_50_SingleMuon_v1.npy")
    #plot_name = "fig_SingleMuon_v1"
    #remake_arrays(in_array, plot_dir, plot_name)
    #
    #in_array  = read_file(data_dir + "design_0p1_no_outer_all_pix_smear_charge_l1000_size_1_50_ZeroBias_v1.npy")
    #plot_name = "fig_ZeroBias_v1"
    #remake_arrays(in_array, plot_dir, plot_name)
    
    in_array  = read_file(data_dir + "design_0p1_no_outer_all_pix_smear_charge_l1000_size_1_50_SingleMuon_v2.npy")
    plot_name = "fig_SingleMuon_v2"
    remake_arrays(in_array, plot_dir, plot_name)
    
    #in_array  = read_file(data_dir + "design_0p1_no_outer_all_pix_smear_charge_l1000_size_1_50_ZeroBias_v2.npy")
    #plot_name = "fig_ZeroBias_v2"
    #remake_arrays(in_array, plot_dir, plot_name)

