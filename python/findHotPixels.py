import ROOT as rt
import numpy as np
from pixelMapping_cfi import *
from scipy.ndimage.filters import convolve
import matplotlib.pyplot as plt
import matplotlib.colors as clrs


def make_map(input_files_, directory):
    """
    function to make a digi occupancy map of layer 1 pixels, with 1:1 pixel, index correspondence
    :param input_files: list of PixelTree ntuples to be processed
    :return:
    """
    chain = rt.TChain('pixelTree')
    for f in input_files_:
        chain.Add(directory + f)

    # creating empty array to be used as map
    occ_array = np.zeros((1920, 3328))

    n_events = chain.GetEntries()

    for iev, event in enumerate(chain):
        if iev % 10 == 0:
            print 'Event', iev, ' / ', n_events
        n_digi = event.DgN

        for iDg in xrange(n_digi):

            detid_digi = event.DgDetid[iDg]

            # skipping detids outside of layer 1
            if not(303000000 < detid_digi < 304000000):
                continue  # layer 1


            row_digi = int(event.DgRow[iDg])
            col_digi = int(event.DgCol[iDg])

            # converting detids to a ladder index and module index
            ladder_index = int(ladder_map(detid_digi))
            module_index = int(module_map(detid_digi, ladder_index))

            occ_array[row_digi+(160*ladder_index)][col_digi+(416*module_index)] += int(1)
    return occ_array


def make_roc_map(input_files_, directory_):
    """
    :param input_files_:
    :param directory_:
    :return: occ_array:
    """

    chain = rt.TChain('pixelTree')
    for f in input_files_:
        chain.Add(directory + f)

    # creating empty array to be used as map
    occ_array = np.zeros((24, 64))

    n_events = chain.GetEntries()

    for iev, event in enumerate(chain):
        if iev % 10 == 0:
            print 'Event', iev, ' / ', n_events

        n_digi = event.DgN

        for iDg in xrange(n_digi):

            detid_digi = event.DgDetid[iDg]

            # skipping detids outside of layer 1
            if not(303000000 < detid_digi < 304000000):
                continue  # layer 1

            row_digi = int(event.DgRow[iDg])
            col_digi = int(event.DgCol[iDg])

            # converting detids to a ladder index and module index
            ladder_index = int(ladder_map(detid_digi))
            module_index = int(module_map(detid_digi, ladder_index))

            roc_row_index = (row_digi // 80) + (2*ladder_index)
            roc_col_index = (col_digi // 52) + (8*module_index)

            occ_array[roc_row_index][roc_col_index] += int(1)
    return occ_array


def make_hot_array(input_array, factor):
    """
    make a list of pixels that are outside of the range: factor * std. dev.
    :param input_array: convolved layer 1 map produced from make_map()
    :param factor: number of std. deviations that will be used to define a hot pixel
                   pixels with values beyond this will be considered as 'hot'/'cold'
    :return:
    """
    stdev = input_array.std()

    epsilon = stdev * factor

    # find and store the indices of all pixels which pass as 'hot' or 'cold'
    # i.e. either more active or less active than the average
    hot_indices = np.where(input_array < epsilon)
    cold_indices = np.where(input_array > epsilon)

    hot_pixel_array = []
    cold_pixel_array = []

    # get row, col, module, and ladder indices and convert to style used in PixelTree Ntuples
    for i in xrange(len(hot_indices[0])):
        row = int(hot_indices[0][i])
        col = int(hot_indices[1][i])
        ladder = int(int(row) / 160)
        module = int(int(col) / 416)

        row_true = row - (160*ladder)
        col_true = col - (416*module)

        if 0 <= module <= 3:
            module_true = module - 4
        elif 4 <= module <= 7:
            module_true = module - 3

        if 0 <= ladder <= 5:
            ladder_true = ladder - 6
        elif 6 <= ladder <= 11:
            ladder_true = ladder - 5

        hot_pixel_array.append([row_true, col_true, module_true, ladder_true])
    # get row, col, module, and ladder indices and convert to style used in PixelTree Ntuples
    for i in xrange(len(cold_indices[0])):
        row = int(cold_indices[0][i])
        col = int(cold_indices[1][i])
        ladder = int(int(row) / 160)
        module = int(int(col) / 416)

        row_true = row - (160*ladder)
        col_true = col - (416*module)

        if 0 <= module <= 3:
            module_true = module - 4
        elif 4 <= module <= 7:
            module_true = module - 3

        if 0 <= ladder <= 5:
            ladder_true = ladder - 6
        elif 6 <= ladder <= 11:
            ladder_true = ladder - 5

        cold_pixel_array.append([row_true, col_true, module_true, ladder_true])
    return hot_pixel_array, cold_pixel_array


def make_hot_roc_array(input_array, factor):
    """
    make a list of pixels that are outside of the range: factor * std. dev.
    :param input_array: convolved layer 1 map produced from make_map()
    :param factor: number of std. deviations that will be used to define a hot pixel
                   pixels with values beyond this will be considered as 'hot'/'cold'
    :return:
    """
    stdev = input_array.std()

    epsilon = stdev * factor

    # find and store the indices of all pixels which pass as 'hot' or 'cold'
    # i.e. either more active or less active than 0
    hot_indices = np.where(input_array < epsilon)
    cold_indices = np.where(input_array > epsilon)
    hot_roc_array = []
    cold_roc_array = []

    # get row, col, module, and ladder indices and convert to style used in PixelTree Ntuples
    for i in xrange(len(hot_indices[0])):
        roc_row = int(hot_indices[0][i])
        roc_col = int(hot_indices[1][i])
        ladder = int(int(roc_row) / 2)
        module = int(int(roc_col) / 8)

        row_true = roc_row - (2*ladder)
        col_true = roc_col - (8*module)

        if 0 <= module <= 3:
            module_true = module - 4
        elif 4 <= module <= 7:
            module_true = module - 3

        if 0 <= ladder <= 5:
            ladder_true = ladder - 6
        elif 6 <= ladder <= 11:
            ladder_true = ladder - 5

        hot_roc_array.append([row_true, col_true, module_true, ladder_true])
    # get row, col, module, and ladder indices and convert to style used in PixelTree Ntuples
    for i in xrange(len(cold_indices[0])):
        roc_row = int(cold_indices[0][i])
        roc_col = int(cold_indices[1][i])
        ladder = int(int(roc_row) / 2)
        module = int(int(roc_col) / 8)

        row_true = roc_row - (2*ladder)
        col_true = roc_col - (8*module)

        if 0 <= module <= 3:
            module_true = module - 4
        elif 4 <= module <= 7:
            module_true = module - 3

        if 0 <= ladder <= 5:
            ladder_true = ladder - 6
        elif 6 <= ladder <= 11:
            ladder_true = ladder - 5

        cold_roc_array.append([row_true, col_true, module_true, ladder_true])
    return hot_roc_array, cold_roc_array


def main(inputFiles, directory, factor, sample_type='realistic'):
    occ = make_map(inputFiles, directory)
    kernel = np.array([[1./8., 1./8., 1./8.],
                       [1./8., -1, 1./8.],
                       [1./8., 1./8., 1./8.]])

    plt.imshow(occ, cmap=plt.cm.inferno, norm=clrs.Normalize(vmin=0.9*np.amin(occ), vmax=1.1*np.amax(occ)))
    plt.savefig('hot_map_precv_'+str(factor)+'_'+sample_type+'_18719.png', bbox_inches='tight')
    plt.close()

    occ_convolved = convolve(occ, kernel, mode='wrap')

    plt.imshow(occ_convolved, cmap=plt.cm.inferno, norm=clrs.Normalize(vmin=0.9*np.amin(occ_convolved), vmax=1.1*np.amax(occ_convolved)))
    plt.savefig('hot_map_'+str(factor)+'_'+sample_type+'_18719.png', bbox_inches='tight')
    plt.close()

    hot_pixels, cold_pixels = make_hot_array(occ_convolved, factor)
    return hot_pixels, cold_pixels


def main_roc(inputFiles, directory, factor, sample_type='realistic'):
    occ = make_roc_map(inputFiles, directory)
    kernel = np.array([[1./8., 1./8., 1./8.],
                       [1./8., -1, 1./8.],
                       [1./8., 1./8., 1./8.]])

    plt.imshow(occ, cmap=plt.cm.inferno, norm=clrs.Normalize(vmin=0.9*np.amin(occ), vmax=1.1*np.amax(occ)))
    plt.savefig('hot_roc_map_precv_'+str(factor)+'_'+sample_type+'_18719.png', bbox_inches='tight')
    plt.close()

    occ_convolved = convolve(occ, kernel, mode='wrap')

    plt.imshow(occ_convolved, cmap=plt.cm.inferno, norm=clrs.Normalize(vmin=0.9*np.amin(occ_convolved), vmax=1.1*np.amax(occ_convolved)))
    plt.savefig('hot_roc_map_'+str(factor)+'_'+sample_type+'_18719.png', bbox_inches='tight')
    plt.close()

    hot_rocs, cold_rocs = make_hot_roc_array(occ_convolved, factor)
    return hot_rocs, cold_rocs


if __name__ == "__main__":
    from inputFiles_cfi import *
    # directory = './files/pixeltrees_mc/design/'
    directory = './files/pixeltrees_mc/realistic/'
    # directory = './files/pixeltrees_mc/RelValTTbar_13TeV/crab_RelValTTbar_13TeVrealistic_0/190710_173122/0000/'
    # directory = './files/pixeltrees_mc/RelValTTbar_13TeV/crab_RelValTTbar_13TeVrealistic_0p1/190710_173131/0000/'
    ## directory = './files/324970/'
    hot_pixel_array, cold_pixel_array = main(input_files_realistic, directory, 3, 'realistic_0p01')
    np.save('hot_pixels_realistic_0p01_3_18719.npy', hot_pixel_array)
    np.save('cold_pixels_realistic_0p01_3_18719.npy', cold_pixel_array)
    hot_roc_array, cold_roc_array = main_roc(input_files_realistic, directory, 3, 'realistic_0p01')
    np.save('hot_rocs_realistic_0p01_3_18719.npy', hot_roc_array)
    np.save('cold_rocs_realistic_0p01_3_18719.npy', cold_roc_array)
    # roc_array = main_roc(input_files_realistic, directory, 2, 'realistic')
    # np.save('hot_rocs_realistic_2_18719.npy', roc_array)

    # pixel_array = main(input_files_realistic, directory, 5, 'realistic')
    # np.save('hot_pixels_realistic_5_18719.npy', pixel_array)
    # roc_array = main_roc(input_files_realistic, directory, 5, 'realistic')
    # np.save('hot_rocs_realistic_5_18719.npy', roc_array)
    # directory = './files/pixeltrees_mc/realistic/'
    # pixel_array = main(input_files_realistic, directory, 5)

    # np.save('hot_pixels_realistic_5_12419.npy', pixel_array)


