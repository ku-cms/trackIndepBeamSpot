import os

def get_list_of_files(dir_):
    file_list_tmp = [os.path.join(dir_, f) for f in os.listdir(dir_) if (os.path.isfile(os.path.join(dir_, f)) and ('.root' in f))]
    return file_list_tmp

