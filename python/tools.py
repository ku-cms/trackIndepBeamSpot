# tools.py

import os
import glob

# get list of local files
def get_list_of_files(directory, pattern="*.root"):
    d = directory
    if d[-1] != "/":
        d += "/"
    return glob.glob(d + pattern)

# creates directory if it does not exist
def makeDir(dir_name):
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)

