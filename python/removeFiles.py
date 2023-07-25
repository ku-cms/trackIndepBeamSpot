# removeFiles.py

import time
import tools
import argparse

def run():
    # options
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--directory",  "-d", default="", help="directory containing root files")
    parser.add_argument("--pattern",    "-p", default="", help="pattern for root file names")

    options     = parser.parse_args()
    directory   = options.directory
    pattern     = options.pattern

    if not directory:
        print("ERROR: 'directory' is not set. Please provide a directory using the -d option.")
        return
    
    if not pattern:
        print("ERROR: 'pattern' is not set. Please provide a pattern using the -p option.")
        return
    
    files = tools.get_eos_file_list(directory)
    n_files = len(files)
    print("Number of files: {0}".format(n_files))

def main():
    t_start = time.time()
    run()
    t_stop  = time.time()
    t_run   = t_stop - t_start
    print("run time (sec): {0:.3f}".format(t_run))

if __name__ == "__main__":
    main()

