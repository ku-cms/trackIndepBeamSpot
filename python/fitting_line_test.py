import numpy as np
import matplotlib.pyplot as plt
import iminuit as im
import string

alpha_low = string.ascii_lowercase

def read_file(input_file_):
    return np.load(input_file_, allow_pickle=True, encoding='latin1')

def remake_arrays(input_arr_, plot_dir, plot_name):
    useWeightedAve  = False
    drawFit         = False
    maskPhi         = False
    fixPhi          = True
    doMinOccCut     = True
    
    # need z-binning corresponding to 1 roc
    w_z_bins = 52  # # of pixels in a roc

    # need phi-binning corresponding to 1 roc (maybe 2?)
    w_phi_bins = 80

    n_z_bins = int(3328 / w_z_bins)  # 3328 is number of pixels in a ladder row
    #print (n_z_bins)
    n_phi_bins = int(960 / w_phi_bins)  # 1440 is number of pixels around phi for all ladders, 960 for inner ladders
    #print (n_phi_bins)
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
        
        # select in phi (0-11)
        #if not roc % 12 == 3:
        #    continue

        # select in z (0-63 after dividing by 12)
        #if not int((roc/12)%64) == 31:
        #    continue
       
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

        # WARNING: phi = -pi = pi issue for roc % 12 == 3
        # - roc 3 in ladder 3 (roc % 12 == 3) crosses phi = -pi = +pi
        # - do not take a normal average
        # - split into phi < 0 and phi >= 0
        # - for phi < 0, find the different from -pi, then add this to +pi for the average
        # - make sure the final average is within [-pi, +pi]... if avg > pi, then it should be set to -pi < new_avg < 0

        # --- apply fix to phi
        # select only roc == 3 for testing
        #if roc == 3:
        # select all roc 3 (mod 12) to fix phi = -pi = +pi issue
        if fixPhi and (roc % 12 == 3):

            phi_neg = [val for val in phi if val <  0.0]
            phi_pos = [val for val in phi if val >= 0.0]
            # for -pi < phi < 0, find absolute value of different from -pi, and add to +pi
            phi_neg_fixed = [np.pi + abs(-np.pi - val) for val in phi if val < 0.0]
            phi_fixed = phi_pos + phi_neg_fixed
            
            # debugging
            #print("You have chosen roc={0}.".format(roc))
            #print("phi: {0}".format(phi))
            #print("number of phi: {0}".format(len(phi)))
            #print("number of phi <  0: {0}".format(len(phi_neg)))
            #print("number of phi >= 0: {0}".format(len(phi_pos)))
            #print("number of phi_fixed: {0}".format(len(phi_fixed)))
            #print("average for phi <  0: {0}".format(np.average(phi_neg)))
            #print("average for phi >= 0: {0}".format(np.average(phi_pos)))
            #print("average for phi_fixed: {0}".format(np.average(phi_fixed)))
            
        else:
            phi_fixed = phi
        
            
        if useWeightedAve:
            occ.append(np.sum(occ_tmp))
            x_array.append(np.average(x,            weights=occ_tmp))
            y_array.append(np.average(y,            weights=occ_tmp))
            z_array.append(np.average(z,            weights=occ_tmp))
            phi_array.append(np.average(phi_fixed,  weights=occ_tmp))
            r_array.append(np.average(r,            weights=occ_tmp))
        else:
            occ.append(np.sum(occ_tmp))
            x_array.append(np.average(x))
            y_array.append(np.average(y))
            z_array.append(np.average(z))
            phi_array.append(np.average(phi_fixed))
            r_array.append(np.average(r))

        # debugging
        #print("roc {0}: z_avg = {1:.3f}, phi_avg = {2:.3f}, phi_fixed_avg = {3:.3f}".format(roc, np.average(z), np.average(phi), np.average(phi_fixed)))

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

    # logic:
    # "*" are "and"
    # "+" are "or"
    # "True" are kept
    # "False" are removed
    
    # limit -25 < z < 25
    remove_z = (z_array > -25) * (z_array < 25)
    # limit -10 < z < 10
    #remove_z = (z_array > -10) * (z_array < 10)
    # limit -5 < z < 5
    #remove_z = (z_array > -5) * (z_array < 5)
    
    remove_blips =  (z_array < -21) + (z_array > -20)
    remove_blips *= (z_array < -14.5) + (z_array > -13.5)
    remove_blips *= (z_array < -7.5) + (z_array > -6.5)
    remove_blips *= (z_array < 5.75) + (z_array > 6.5)
    remove_blips *= (z_array < 12.5) + (z_array > 13.5)
    remove_blips *= (z_array < 19) + (z_array > 20)
    # remove -10 < z < 10
    #remove_blips *= (z_array > 10) + (z_array < -10)
    # remove -5 < z < 5
    #remove_blips *= (z_array > 5) + (z_array < -5)

    occ = occ[remove_z*remove_blips]
    x_array = x_array[remove_z*remove_blips]
    y_array = y_array[remove_z*remove_blips]
    z_array = z_array[remove_z*remove_blips]
    phi_array = phi_array[remove_z*remove_blips]

    # remove slice in phi
    if maskPhi:
        remove_phi          = (phi_array >= -np.pi) * (phi_array <= np.pi)
        remove_blips_phi    = (phi_array < 1.0) + (phi_array > 1.5)
        
        occ         = occ[          remove_phi * remove_blips_phi]
        x_array     = x_array[      remove_phi * remove_blips_phi]
        y_array     = y_array[      remove_phi * remove_blips_phi]
        z_array     = z_array[      remove_phi * remove_blips_phi]
        phi_array   = phi_array[    remove_phi * remove_blips_phi]

    # remove low occupancy clusters
    if doMinOccCut:
        min_occupancy = 20000 
        occupancy_cut = occ >= min_occupancy
        occ         = occ[       occupancy_cut ]
        x_array     = x_array[   occupancy_cut ]
        y_array     = y_array[   occupancy_cut ]
        z_array     = z_array[   occupancy_cut ]
        phi_array   = phi_array[ occupancy_cut ]

    def nll(x0, y0, z0, n, b1, b2, b3, a1, a3, c1, c3, ga1, ga3, gc1, gc3):
        ri = np.float64(np.sqrt((x_array - x0) ** 2 + (y_array - y0) ** 2 + (z_array - z0) ** 2))
        phi_cor = np.arctan2(y_array-y0, x_array-x0)
        # note: a2 = b2 = c2
        a  = np.float64(b_par(phi_cor, a1, b2, a3))
        b  = np.float64(b_par(phi_cor, b1, b2, b3))
        c  = np.float64(b_par(phi_cor, c1, b2, c3))
        ga = np.float64(b_par(phi_cor, ga1, b2-np.pi, ga3))
        gc = np.float64(b_par(phi_cor, gc1, b2-np.pi, gc3))
        gb = z0

        func_array = -np.log(1 / (np.sqrt(2 * np.pi * occ))) + (n * (func(ri, a, b, c) + line(z_array, ga, gb, gc)) - occ) ** 2 / (2 * np.pi * occ)

        return np.sum(func_array)

    # --- no gauss version from kim --- #
    # minuit = im.Minuit(nll, x0=0, y0=0, z0=0, n=1,
    #                    a1=4e4, a3=1.26e5,
    #                    b1=0.0, b2=0.5, b3=1.167,
    #                    c1=1e3, c3=3.12e3,
    #                    error_x0=0.001, error_y0=0.001, error_z0=0.001, error_n=0.01,
    #                    error_b1=0.01, error_b2=0.01, error_b3=0.01,
    #                    error_a1=10, error_a3=0.1,
    #                    error_c1=10, error_c3=0.1,
    #                    fix_n=True,
    #                    fix_a1=False, fix_c1=False,
    #                    limit_x0=(-1, 1), limit_y0=(-1, 1), limit_z0=(-1, 1),
    #                    limit_b2=(-np.pi, np.pi),
    #                    limit_b3=(0, 3),
    #                    errordef=1
    # )
    
    # --- gauss version from erich --- #
    minuit = im.Minuit(nll,
                       #x0=0., y0=0., z0=0., n=1,
                       x0=0.1, y0=-0.08, z0=-0.3, n=1,
                       # (0, 0) fit
                       a1=-0.19e4,  a3=0.86e5,
                       b1=0,        b2=0,       b3=1.212,
                       c1=-110,     c3=1260,
                       ga1=-500,    ga3=-1.2e4,
                       gc1=-0.031,  gc3=2.129,
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
                       
                       error_x0=0.001,  error_y0=0.001, error_z0=0.001, error_n=0.01,
                       error_b1=0.01,   error_b2=0.01,  error_b3=0.01,
                       error_a1=1,      error_a3=0.1,
                       error_c1=1,      error_c3=0.1,
                       error_ga1=0.1,   error_ga3=0.1,
                       error_gc1=0.1,   error_gc3=0.1,
                       
                       fix_n  = True,
                       fix_x0 = False,
                       fix_y0 = False,
                       fix_z0 = True,
                       
                       limit_x0=(-1, 1), limit_y0=(-1, 1), limit_z0=(-2, 2),
                       limit_b2=(-np.pi, np.pi),
                       limit_b3=(0, 3),
                       errordef=1
    )

    minuit.migrad()
    print(minuit.get_param_states())
    print(minuit.get_fmin())
    
    #print(" - first access") 
    print(minuit.values)
    #print(minuit.errors)
    #print(minuit.covariance)
    #print(minuit.values["x0"])
    
    #print(" - second access") 
    #print(minuit.np_values())
    #print(minuit.np_errors())
    #print(minuit.np_covariance())
    #print(minuit.np_values()[0])

    # get values after fit
    x0  = minuit.values["x0"]
    y0  = minuit.values["y0"]
    z0  = minuit.values["z0"]
    n   = minuit.values["n"]
    b1  = minuit.values["b1"]
    b2  = minuit.values["b2"]
    b3  = minuit.values["b3"]
    a1  = minuit.values["a1"]
    a3  = minuit.values["a3"]
    c1  = minuit.values["c1"]
    c3  = minuit.values["c3"]
    ga1 = minuit.values["ga1"]
    ga3 = minuit.values["ga3"]
    gc1 = minuit.values["gc1"]
    gc3 = minuit.values["gc3"]
    
    #x = np.linspace(-25, 25, 30)
    x = np.linspace(-10, 10, 30)
    y = np.linspace(-np.pi, np.pi, 30)
    
    # X, Y, and Z are 2D matrices
    X, Y = np.meshgrid(x, y)
    Z = func_expanded(X, Y, a1, a3, b1, b2, b3, c1, c3) + line_expanded(X, Y, z0, a1, a3, b1, b2, b3, c1, c3, ga1, ga3, gc1, gc3)
    
    axs.plot(z_array, phi_array, occ, 'b*')
    if drawFit:
        axs.plot_wireframe(X, Y, Z, color='black')
    
    # labels
    axs.set_title("Pixel Occupancy", fontsize=20)
    axs.set_xlabel("z",              fontsize=16)
    axs.set_ylabel(r"$\phi$",        fontsize=16)

    #axs.legend()
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

def gauss(x, ga, gb, gc):
    return ga * np.exp(-0.5 * (x-gb)**2 / gc**2)

def line(x, ga, gb, gc):
    return ga * abs(x - gb) + gc

def func_expanded(x, y, a1, a3, b1, b2, b3, c1, c3):
    # note: a2 = b2 = c2
    a = b_par(y, a1, b2, a3)
    b = b_par(y, b1, b2, b3)
    c = b_par(y, c1, b2, c3)
    # r_2d = sqrt(x ** 2 + y ** 2)
    # r_3d = sqrt(x ** 2 + y ** 2 + z ** 2)
    # r_2d = 1.8 cm
    # r_3d = sqrt(r_2d ** 2 + z ** 2)
    r_2d = 1.8
    r_3d = np.sqrt(r_2d ** 2 + x ** 2)
    return func(r_3d, a, b, c)

def line_expanded(x, y, z0, a1, a3, b1, b2, b3, c1, c3, ga1, ga3, gc1, gc3):
    ga = b_par(y, ga1, b2-np.pi, ga3)
    gc = b_par(y, gc1, b2-np.pi, gc3)
    gb = z0
    return line(x, ga, gb, gc)

def f_example(x, y):
    a1 = 0.003e6
    a3 = 0.076e6
    b1 = 0.47
    b2 = 0.95
    b3 = 1.305
    c1 = 0.8e3
    c3 = 1.0e3
    return (a1*np.sin(y-b2)+a3) * (1 / np.abs(x) ** (b1*np.sin(y-b2)+b3)) + (c1*np.sin(y-b2)+c3)

def runSet(data_dir, plot_dir, data_name, plot_name=""):
    # default: data and plot names are the same
    # if plot name is given, it is used for plots
    custom_name = data_name
    if plot_name:
        custom_name = plot_name
    in_array = read_file(data_dir + data_name + ".npy")
    remake_arrays(in_array, plot_dir, custom_name)

def runZeroBias2017B(data_dir, plot_dir):
    # ZeroBias 2017B
    data_names = [
        "ZeroBias_2017B_AllClusters",
        "ZeroBias_2017B_ClusterSize2_AllClusters",
        "ZeroBias_2017B_ClusterSize2_NumberClusters2000_AllClusters",
    ]
    #plot_names = [
    #    "ZeroBias_2017B_NoOccMin_AllClusters",
    #    "ZeroBias_2017B_ClusterSize2_NoOccMin_AllClusters",
    #    "ZeroBias_2017B_ClusterSize2_NumberClusters2000_NoOccMin_AllClusters",
    #]
    plot_names = [
        "ZeroBias_2017B_MinOcc20000_AllClusters",
        "ZeroBias_2017B_ClusterSize2_MinOcc20000_AllClusters",
        "ZeroBias_2017B_ClusterSize2_NumberClusters2000_MinOcc20000_AllClusters",
    ]
    for i in range(len(data_names)):
        runSet(data_dir, plot_dir, data_names[i], plot_names[i])

def runSingleMuon2017B(data_dir, plot_dir):
    # SingleMuon 2017B
    names = [
        "SingleMuon_2017B_MoreEvents_AllClusters",
        #"SingleMuon_2017B_MoreEvents_ClusterSize2_AllClusters",
    ]
    for name in names:
        runSet(data_dir, plot_dir, name)

def runSingleMuon2018C(data_dir, plot_dir):
    # SingleMuon 2018C
    names = [
        "SingleMuon_MoreEvents_AllClusters",
        #"SingleMuon_MoreEvents_ClusterSize2_AllClusters",
        #"SingleMuon_AllClusters",
        #"SingleMuon_OnTrack",
        #"SingleMuon_OffTrack",
    ]
    for name in names:
        runSet(data_dir, plot_dir, name)

def runZeroBias(data_dir, plot_dir):
    # ZeroBias
    names = [
        "ZeroBias_MoreEvents2_AllClusters",
        "ZeroBias_MoreEvents2_ClusterSize2_AllClusters",
        #"ZeroBias_MoreEvents_AllClusters",
        #"ZeroBias_MoreEvents_ClusterSize2_AllClusters",
        #"ZeroBias_AllClusters",
        #"ZeroBias_OnTrack",
        #"ZeroBias_OffTrack",
    ]
    for name in names:
        runSet(data_dir, plot_dir, name)

def runMinBias2017B(data_dir, plot_dir):
    # MinBias 2017B
    # WARNING: MinBias 2017B does not have good data. Errors can occur.
    names = [
        "MinBias_2017B_AllClusters",
        "MinBias_2017B_OnTrack",
        "MinBias_2017B_OffTrack",
    ]
    for name in names:
        runSet(data_dir, plot_dir, name)

def runMinBias2018C(data_dir, plot_dir):
    # MinBias 2018C
    names = [
        "MinBias_2018C_MoreEvents2_AllClusters",
        "MinBias_2018C_MoreEvents2_ClusterSize2_AllClusters",
        #"MinBias_2018C_MoreEvents_AllClusters",
        #"MinBias_2018C_MoreEvents_ClusterSize2_AllClusters",
        #"MinBias_2018C_AllClusters",
        #"MinBias_2018C_OnTrack",
        #"MinBias_2018C_OffTrack",
    ]
    for name in names:
        runSet(data_dir, plot_dir, name)
    
def runTTBar(data_dir, plot_dir):
    # TTBar
    names = [
        "TTBar_0p2_AllClusters",
        "TTBar_0p2_ClusterSize2_AllClusters",
        "TTBar_pileup_0p2_AllClusters",
        "TTBar_pileup_0p2_ClusterSize2_AllClusters",
        "TTBar_AllClusters",
        "TTBar_OnTrack",
        "TTBar_OffTrack",
        "TTBar_AllClusters_zsmear",
        "TTBar_OnTrack_zsmear",
        "TTBar_OffTrack_zsmear",
    ]
    for name in names:
        runSet(data_dir, plot_dir, name)

if __name__ == "__main__":
    data_dir  = "data/"
    plot_dir  = "plots/"
    
    if data_dir[-1] != "/":
        data_dir += "/"
    
    if plot_dir[-1] != "/":
        plot_dir += "/"
    
    runZeroBias2017B(data_dir, plot_dir)
    #runSingleMuon2017B(data_dir, plot_dir)
    #runZeroBias(data_dir, plot_dir)
    #runSingleMuon2018C(data_dir, plot_dir)
    #runMinBias2017B(data_dir, plot_dir)
    #runMinBias2018C(data_dir, plot_dir)
    #runTTBar(data_dir, plot_dir)
