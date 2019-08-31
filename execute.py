#!/usr/bin/python
"""
Script to execute and compare runs
Andrew Lamzed-Short
u1897268
"""

import argparse
import csv
import datetime
import os
import re
import subprocess
import sys
import time


gmon_regex = "gmon\.out.*"


def delete_files(dir, pattern):
    for f in os.listdir(dir):
        if re.match(pattern, f):
            os.remove(os.path.join(dir, f))

def avg_times():
    file_regex = ".*time_\d\.txt$"
    time_regex = "Elapsed time: (?P<time>\d+\.\d+) seconds"

    cd = os.getcwd()
    out_path = os.path.join(cd, "out")

    paths_to_process = []
    partitioned_sets = []

    # Filter out the files we need
    for f in os.listdir(out_path):
        if re.match(file_regex, f):
            paths_to_process.append(os.path.join(out_path, f))

    # Partition into related groups
    paths_to_process = sorted(paths_to_process)
    size = 3
    partitioned_sets = [paths_to_process[i:i + size] for i  in range(0, len(paths_to_process), size)]
    # print(partitioned_sets)

    # Read and average the times
    # Output the result
    with open("execution.csv", "wb") as csv_file:
        writer = csv.writer(csv_file)
        i = 1
        for g in partitioned_sets:
            t = 0
            for f in g:
                with open(f, "r") as fi:
                    match = re.search(time_regex, fi.read())
                    if match:
                        t += float(match.group('time'))
            writer.writerow([i, t / 3])
            i += 1

    


# Setup argument parsing
parser = argparse.ArgumentParser(description="Executes snap with given parameters.")

parser.add_argument("--setup", action="store_true")
parser.add_argument("--clean", action="store_true")

parser.add_argument("--avg_times", action="store_true")

parser.add_argument("--exe", help="Location of executable to run")
parser.add_argument("--np", help="Number of processes to run")
parser.add_argument("--fi", help="Path to input file to read")
parser.add_argument("--fo", help="Path to output file to write")

args = parser.parse_args()


if __name__ == "__main__":
    print("Executing main program")


    if args.avg_times:
        avg_times()
        sys.exit(0)

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
        sys.exit(1)

    if args.np is None:
        print("np not set")
        sys.exit(1)

    if args.fi is None:
        print("fi not set")
        sys.exit(1)

    if args.fo is None:
        print("fo not set")
        sys.exit(1)

    # Execute and time
    start_time = time.time()
    subprocess.call(["mpirun", "-np", args.np, args.exe, "--fi", args.fi, "--fo", args.fo])
    end_time   = time.time()
    print("Elapsed time: {} seconds".format(end_time - start_time))

    sys.exit(0)