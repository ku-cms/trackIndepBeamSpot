import os

def get_file_list(dir_):
    file_list_tmp = [os.path.join(dir_, f) for f in os.listdir(dir_) if (os.path.isfile(os.path.join(dir_, f)) and ('.root' in f))]
    return file_list_tmp

def eosls(path, option = "", eosurl = "root://cmseos.fnal.gov"):
    return os.popen("xrdfs %s ls %s %s"%(eosurl, option, path))

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

