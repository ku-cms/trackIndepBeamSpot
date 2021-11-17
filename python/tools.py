# tools.py

import os
import glob

# creates directory if it does not exist
def makeDir(dir_name):
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)

# get list of local files
def get_file_list(dir_):
    file_list_tmp = [os.path.join(dir_, f) for f in os.listdir(dir_) if (os.path.isfile(os.path.join(dir_, f)) and ('.root' in f))]
    return file_list_tmp

# get list of local files using glob
def get_file_list_glob(directory, pattern="*.root"):
    d = directory
    if d[-1] != "/":
        d += "/"
    return glob.glob(d + pattern)

# get list of EOS files
def get_eos_file_list(path, eosurl = "root://cmseos.fnal.gov"):
    output = [] 
    with eosls(path, "", eosurl) as files:
        for f in files:
            name = f.strip()
            # add ROOT files to list
            if name.endswith(".root"):
                full_name = "{0}/{1}".format(eosurl, name)
                output.append(full_name)
    return output

# eosls command
def eosls(path, option = "", eosurl = "root://cmseos.fnal.gov"):
    return os.popen("xrdfs %s ls %s %s"%(eosurl, option, path))


