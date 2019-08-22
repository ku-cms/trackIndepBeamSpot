#!/usr/bin/env python

import ROOT as rt


def read_event_clusters(input_files):
    """

    :param input_files: a list of files names
    """
    tChain = rt.TChain('pixelTree')
    for f in input_files:
        tChain.Add(f)

    of = rt.TFile(input_files[0].split('_')[0] + '_out_fitting.root', 'recreate')

    nEvents = tChain.GetEntries()

    # clustRZ_all = []
    # clustRPhi_all = []
    # clustPhiZ_all = []
    clustPhiZ_array = []
    for iev in xrange(nEvents):
        #     clustRZ_all.append(rt.TH2F('hClustRZ_all_'+str(iev), 'Clusters R vs Z', 500, -28., 28, 500, 0, 20.))
        #     clustRPhi_all.append(rt.TH2F('hClustRPhi_all_'+str(iev), 'Clusters R vs phi', 500, -rt.TMath.Pi()-1, rt.TMath.Pi()+1, 500, 0, 20.))
        #     clustPhiZ_all.append(rt.TH2F('hClustPhiZ_all_'+str(iev), 'Clusters phi vs Z', 500, -rt.TMath.Pi()-1, rt.TMath.Pi()+1,  500, -28., 28.,))
        clustPhiZ_array.append([])

    nClusTot = 0.
    for iev in xrange(5):
        if iev % 1 == 0:
            print('Event: ', iev)
        tChain.GetEntry(iev)
        nClus = tChain.ClN
        nClusTot += nClus
        if iev == 0:
            print(nClus)
        for iCl in xrange(nClus):
            layer = tChain.ClLayer[iCl]
            clusGX = tChain.ClGx[iCl]
            clusGY = tChain.ClGy[iCl]
            clusGZ = tChain.ClGz[iCl]
            clPhi = np.arctan(clusGX / clusGY)

            clR = np.sqrt(clusGX ** 2 + clusGY ** 2)
            if clusGX < 0 < clusGY:
                clPhi += np.pi
            elif clusGX < 0 and clusGY < 0:
                clPhi -= np.pi

            if layer == 2 or layer == 3 or layer == 4:
                #if 0 < clusGZ < 2:
                #clustRPhi_all[iev].Fill(clPhi, clR)
                #if 0 < clPhi < 0.5:
                #clustRZ_all[iev].Fill(clusGZ, clR)
                #clustPhiZ_all[iev].Fill(clPhi, clusGZ)
                clustPhiZ_array[iev].append([clPhi, clusGZ, clR])
                continue
            if layer != 1: continue
            #if 0 < clusGZ < 2:
            #clustRPhi_all[iev].Fill(clPhi, clR)
            #if 0 < clPhi < 0.5:
            #clustRZ_all[iev].Fill(clusGZ, clR)
            #clustPhiZ_all[iev].Fill(clPhi, clusGZ)
            clustPhiZ_array[iev].append([clPhi, clusGZ, clR])

    of.Write()
    of.Close()
    np.save('phi_z_r_array.npy', np.array(clustPhiZ_array))

if __name__ == "__main__":
    input_files = ['PixelTree_317392_0.root']
    read_event_clusters(input_files)


