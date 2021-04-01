def get_list_of_files(dir_):
    import os
    file_list_tmp = [os.path.join(dir_, f) for f in os.listdir(dir_) if (os.path.isfile(os.path.join(dir_, f)) and ('.root' in f))]
    #ind = 0
    #for f in xrange(len(file_list_tmp)):
    #    if file_list_tmp[f].split('/')[-1] == 'PixelTree_92.root':
    #        ind = f
    #file_list_tmp = file_list_tmp[ind:]
    return file_list_tmp
