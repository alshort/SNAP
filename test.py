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

VERBOSE = False

gmon_regex = "gmon\.out.*"

cpu_path = "src/snap_cpu"
gpu_path = "src/snap"


def printv(*str):
    if VERBOSE:
        print(*str)

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
    printv(partitioned_sets)

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
parent_parser = argparse.ArgumentParser(description="Executes snap with given parameters.")
parent_parser.add_argument("-v", "--verbose", action="store_true", help="Use verbose output")
subparsers = parent_parser.add_subparsers(title="commands", dest="command")

parser_setup = subparsers.add_parser("setup", add_help=False,
    description="setup", help="sets up the executing environment")

parser_clean = subparsers.add_parser("clean", add_help=False,
    description="clean", help="cleans the working directory")

parser_run = subparsers.add_parser("run",
    description="run", help="runs one of the SNAP program")
parser_run.add_argument("--cpu", action="store_true", help="Use CPU-based build")
parser_run.add_argument("--gpu", action="store_true", help="Use GPU-based build")
parser_run.add_argument("-n", "--np", help="Number of processes to run")
parser_run.add_argument("-i", "--fi", help="Path to input file to read")
parser_run.add_argument("-o", "--fo", help="Path to output file to write")

parser_analyse = subparsers.add_parser("analyse",
    description="analyse", help="performs analysis on the pre-generated output of SNAP executions")

args = parent_parser.parse_args()


def run():
    # Validation
    if args.np is None:
        print("ERROR: --np is not set.")
        parser_run.print_help()
        sys.exit(1)

    if args.fi is None:
        print("ERROR: --fi is not set.")
        parser_run.print_help()
        sys.exit(1)

    if args.fo is None:
        print("ERROR: --fo is not set.")
        parser_run.print_help()
        sys.exit(1)

    # Execute and time
    app = cpu_path
    if (args.cpu is None and args.gpu is None) or (args.cpu is not None and args.gpu is not None):
        print("ERROR: Must set one of either --cpu and --gpu, not both nor neither")
        parser_run.print_help()
        sys.exit(1)

    if args.cpu:
        app = cpu_path
    elif args.gpu:
        app = gpu_path

    start_time = time.time()
    subprocess.call(["mpirun", "-np", args.np, app, "--fi", args.fi, "--fo", args.fo])
    end_time   = time.time()
    print("Elapsed time: {} seconds".format(end_time - start_time))

def analyse():
    avg_times()


if __name__ == "__main__":
    if args.verbose: VERBOSE = True

    printv("Command:", args.command)

    if args.command == None:
        parent_parser.print_help()
        sys.exit(1)

    if args.command == "setup":
        printv("executing setup")
        subprocess.call([".", "~/intel/bin/compilervars.sh", "-arch", "intel64", "-platform", "linux"], shell=True)
        subprocess.call([".", "~/intel/mkl/bin/mklvars.sh", "intel64"], shell=True)
        subprocess.call([".", "~/intel/bin/iccvars.sh", "-arch", "intel64", "-platform", "linux"], shell=True)
        subprocess.call([".", "~/intel/compilers_and_libraries/linux/mpi/intel64/bin/mpivars.sh"], shell=True)
        sys.exit(0)

    elif args.command == "clean":
        printv("executing clean")
        delete_files(os.getcwd(), gmon_regex)
        sys.exit(0)

    elif args.command == "run":
        printv("executing run")
        run()
        sys.exit(0)

    elif args.command == "analyse":
        printv("executing analyse")
        analyse()
        sys.exit(0)
