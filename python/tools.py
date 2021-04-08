# tools.py
import glob

def get_list_of_files(directory, pattern="*.root"):
    d = directory
    if d[-1] != "/":
        d += "/"

    return glob.glob(d + pattern)

