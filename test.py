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

def exist_in_folder(folder, regex):
    for f in os.listdir(folder):
        if re.match(regex, f): return True
    return False

def count_matching_regex_in_folder(folder, regex):
    count = 0
    for f in os.listdir(folder):
        if re.match(regex, f): count += 1
    return count


def parse_arguments():
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

    parser_analyse = subparsers.add_parser("analyse",
        description="analyse", help="performs analysis on the pre-generated output of SNAP executions")

    return parent_parser #.parse_args()

def setup(args):
    subprocess.call([".", "~/intel/bin/compilervars.sh", "-arch", "intel64", "-platform", "linux"], shell=True)
    subprocess.call([".", "~/intel/mkl/bin/mklvars.sh", "intel64"], shell=True)
    subprocess.call([".", "~/intel/bin/iccvars.sh", "-arch", "intel64", "-platform", "linux"], shell=True)
    subprocess.call([".", "~/intel/compilers_and_libraries/linux/mpi/intel64/bin/mpivars.sh"], shell=True)

def clean(args):
    delete_files(os.getcwd(), gmon_regex)

def delete_files(dir, pattern):
    for f in os.listdir(dir):
        if re.match(pattern, f):
            os.remove(os.path.join(dir, f))

def run(args, parser):
    OUT_DIR = "./out/"

    # Validation
    if args.np is None:
        print("ERROR: --np is not set.")
        parser.print_help()
        sys.exit(1)

    if args.fi is None:
        print("ERROR: --fi is not set.")
        parser.print_help()
        sys.exit(1)

    # Execute and time
    app = cpu_path
    if ((args.cpu is False) and (args.gpu is False)) or ((args.cpu is True) and (args.gpu is True)):
        print("ERROR: Must set one of either --cpu and --gpu, not both nor neither")
        parser.print_help()
        sys.exit(1)

    if args.cpu:
        app = cpu_path
    elif args.gpu:
        app = gpu_path

    # Calculate the path of the file the output to
    in_path = os.path.basename(args.fi)
    input_regex = "in(\d\d)"
    out_num = re.match(input_regex, in_path).group(1)

    out_name = "out{}_{}_n{}".format(out_num, "cpu" if args.cpu else "gpu", args.np)
    matching = count_matching_regex_in_folder(OUT_DIR, out_name + "_\d")
    run_id = 1 if matching == 0 else matching + 1

    out_time_name = out_name
    if matching == 0:
        out_name = out_name + "_1")
        out_time_name = out_time_name + "_time_1"
    else:
        printv(out_name + "_{}".format(matching + 1))    

    start_time = time.time()
    subprocess.call(["mpirun", "-np", args.np, app, "--fi", args.fi, "--fo", os.path.join(OUT_DIR, out_name)])
    end_time   = time.time()
    print("Elapsed time: {} seconds".format(end_time - start_time))

    # Times
    if not os.path.isfile("times.csv"):
        with open("times.csv", "wb") as csv_file:
            writer = csv.writer(csv_file)
            writer.writerow(["input_file,run,np,cpu,gpu,block_num,tpb,time"])

    with open("times.csv", "wb") as csv_file:
        writer = csv.writer(csv_file)
        c = 0 if args.gpu else 1
        g = 0 if args.cpu else 1
        writer.writerow([in_path, run_id, args.np, c, g, 0, 0, end_time-start_time])

def analyse(args):
    avg_times()

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


if __name__ == "__main__":
    parser = parse_arguments()
    args = parser.parse_args()
    if args.verbose: VERBOSE = True

    printv("Command:", args.command)

    if args.command == None:
        parent_parser.print_help()
        sys.exit(1)

    if args.command == "setup":
        printv("executing setup")
        setup(args)
        sys.exit(0)

    elif args.command == "clean":
        printv("executing clean")
        clean(args)
        sys.exit(0)

    elif args.command == "run":
        printv("executing run")
        run(args, parser)
        sys.exit(0)

    elif args.command == "analyse":
        printv("executing analyse")
        analyse(args)
        sys.exit(0)
