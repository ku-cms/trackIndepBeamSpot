# removeFiles.py

import os
import time
import tools
import argparse

# TODO
# - Add min/max number range as arguments for files to delete.
# - Print files in min/max range and count files.
# - Add flag to run deletion (e.g. -f). 
# - Add eosrm command.
# DONE
# - List all files in a directory.
# - Extract number from file name.

# Get file number from the file name.
def getFileNumber(file_name):
    # Remove file extension from the file name.
    file_name_no_ext = os.path.splitext(file_name)[0]
    
    # Get numbers from file name without the file extension.
    delimiter = '_'
    numbers = tools.getNumbers(file_name_no_ext, delimiter)
    n_numbers = len(numbers)
    
    # Check for exactly one number.
    if n_numbers == 0:
        print("ERROR: No numbers found in the file name '{0}' using the delimiter '{1}'.".format(file_name, delimiter))
        return -999
    elif n_numbers > 1:
        print("ERROR: Multiple numbers ({0} numbers) found in the file name '{1}' using the delimiter '{2}'.".format(n_numbers, file_name, delimiter))
        return -999
    else:
        return numbers[0]

def removeEOSFiles():
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
    
    # match pattern to base file name, but save full path to file
    files = tools.get_eos_file_list(directory)
    files_matching = [f for f in files if pattern in os.path.basename(f)]
    
    n_files_total       = len(files)
    n_files_matching    = len(files_matching)

    print("directory: {0}".format(directory))
    print("pattern: {0}".format(pattern))
    print("----------------------------")
    
    print("matching files:")
    for f in files_matching:
        file_name = os.path.basename(f)
        file_number = getFileNumber(file_name)
        print(" - {0}; file number: {1}".format(file_name, file_number))
    
    print("----------------------------")
    print("Number of files (total): {0}".format(n_files_total))
    print("Number of files containing the pattern '{0}': {1}".format(pattern, n_files_matching))

def main():
    t_start = time.time()
    removeEOSFiles()
    t_stop  = time.time()
    t_run   = t_stop - t_start
    print("run time (sec): {0:.3f}".format(t_run))

if __name__ == "__main__":
    main()

