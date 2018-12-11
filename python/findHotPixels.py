import ROOT as rt
import numpy as np
from pixelMapping_cfi import *
from scipy.ndimage.filters import convolve
import matplotlib.pyplot as plt


def make_map(input_files):
    """
    function to make a digi occupancy map of layer 1 pixels, with 1:1 pixel, index correspondence
    :param input_files: list of PixelTree ntuples to be processed
    :return:
    """
    chain = rt.TChain('pixelTree')
    for f in input_files:
        chain.Add(f)

    # creating empty array to be used as map
    occ_array = np.zeros((1920, 3328))

    for iev, event in enumerate(chain):
        if iev % 1 == 0:
            print 'Event', iev

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

            occ_array[row_digi+(160*ladder_index)][col_digi+(416*module_index)] += 1
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
    indices = np.where(abs(input_array) > epsilon)

    pixel_array = []

    # get row, col, module, and ladder indices and convert to style used in PixelTree Ntuples
    for i in xrange(len(indices[0])):
        row = int(indices[0][i])
        col = int(indices[1][i])
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

        pixel_array.append([row_true, col_true, module_true, ladder_true])
    return pixel_array


def main(inputFiles, factor):
    occ = make_map(inputFiles)
    kernel = np.array([[1./8., 1./8., 1./8.],
                       [1./8., -1, 1./8.],
                       [1./8., 1./8., 1./8.]])

    occ_convolved = convolve(occ, kernel, mode='constant', cval=0.0)

    plt.imshow(occ_convolved, cmap=plt.cm.gray)
    plt.savefig('hot_map_'+str(factor)+'_61218.jpg', bbox_inches='tight')
    plt.close()

    pixels = make_hot_array(occ_convolved, factor)
    return pixels


if __name__ == "__main__":
    from inputFiles_cfi import *

    pixel_array = main(input_files_2, 5)

    np.save('hot_pixels_5_61218.npy', pixel_array)



