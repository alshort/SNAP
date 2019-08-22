#!/usr/bin/python
"""
Script to execute and compare runs
Andrew Lamzed-Short
u1897268
"""

import argparse
import datetime
import os
import re
import subprocess
import sys


gmon_regex = "gmon\.out-.*"


def delete_files(dir, pattern):
    for f in os.listdir(dir):
        if re.match(pattern, f):
            os.remove(os.path.join(dir, f))


# Setup argument parsing
parser = argparse.ArgumentParser(description="Executes snap with given parameters.")

parser.add_argument("--setup", action="store_true")
parser.add_argument("--clean", action="store_true")

parser.add_argument("--exe", help="Location of executable to run")
parser.add_argument("--fi", help="Path to input file to read")
parser.add_argument("--fo", help="Path to output file to write")

args = parser.parse_args()


if __name__ == "__main__":
    print("Executing main program")

    if args.setup:
        print("Initiate setup")
        subprocess.call([".", "~/intel/bin/compilervars.sh", "-arch", "intel64", "-platform", "linux"], shell=True)
        subprocess.call([".", "~/intel/mkl/bin/mklvars.sh", "intel64"], shell=True)
        subprocess.call([".", "~/intel/bin/iccvars.sh", "-arch", "intel64", "-platform", "linux"], shell=True)
        subprocess.call([".", "~/intel/compilers_and_libraries/linux/mpi/intel64/bin/mpivars.sh"], shell=True)

    if args.clean:
        delete_files(os.getcwd(), gmon_regex)

    # Validation
    if args.exe is None:
        print("exe not set")
    if args.fi is None:
        print("fi not set")
    if args.fo is None:
        print("fo not set")

    sys.exit(0)