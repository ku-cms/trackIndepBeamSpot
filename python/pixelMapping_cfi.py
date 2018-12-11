#!/usr/bin/env python
import numpy as np

ladder_module_map = np.array([[32, 36, 40, 44, 48, 52, 56, 60, 64, 20, 24, 28],
                              [36, 40, 44, 48, 52, 56, 60, 64, 68, 24, 28, 32],
                              [40, 44, 48, 52, 56, 60, 64, 68, 72, 28, 32, 36],
                              [44, 48, 52, 56, 60, 64, 68, 72, 76, 32, 36, 40],
                              [48, 52, 56, 60, 64, 68, 72, 76, 80, 36, 40, 44],
                              [52, 56, 60, 64, 68, 72, 76, 80, 84, 40, 44, 48],
                              [56, 60, 64, 68, 72, 76, 80, 84, 88, 44, 48, 52],
                              [60, 64, 68, 72, 76, 80, 84, 88, 92, 48, 52, 56]])

bin_edge_map = np.array([[-3, -2, 1.08],
                         [-2, -1, 1.35],
                         [-1, 0.2, 1.58],
                         [0.2, 1.1, 1.54],
                         [1.1, 2.2, 1.33],
                         [2.2, np.pi, 1.35]])


def module_map(detid, ladder_index):
    last_digits = detid % 100
    module = ladder_module_map[:, ladder_index].tolist().index(last_digits)
    return module


def r_color_map(r):
    return -12. * r + 279.


def ladder_map(detid):
    if 303075300 < detid < 303075400:
        ladder = 0
    elif 303071200 < detid < 303071300:
        ladder = 1
    elif 303067100 < detid < 303067200:
        ladder = 2
    elif 303063000 < detid < 303063100:
        ladder = 3
    elif 303058900 < detid < 303059000:
        ladder = 4
    elif 303054800 < detid < 303054900:
        ladder = 5
    elif 303050700 < detid < 303050800:
        ladder = 6
    elif 303046600 < detid < 303046700:
        ladder = 7
    elif 303042500 < detid < 303042600:
        ladder = 8
    elif 303087600 < detid < 303087700:
        ladder = 9
    elif 303083500 < detid < 303083600:
        ladder = 10
    elif 303079400 < detid < 303079500:
        ladder = 11

    return ladder


def rocPhi(roc, ladder):
    link = int(roc)/8
    return float(ladder) - 0.5 + (0.5/2.) + (float(link) * 0.5)


def rocID(col, row):
    rocRow = int(row)/80
    rocCol = int(col)/52
    rocID = rocCol + rocRow*8
    return int(rocID)


def ladderType(ladder):
    if ladder % 2 == 0:
        ladType = 'inner'
    else:
        ladType = 'outer'


