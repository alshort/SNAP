#!/usr/bin/python3
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

from functools import reduce


VERBOSE = False

gmon_regex = "gmon\.out.*"

OUT_DIR = "out/"

cpu_path = "src/snap_cpu"
gpu_path = "src/snap"

test_data_path = "testing.csv"
schema = ["input_file", "run", "np", "cpu", "gpu", "block_num", "tpb", "time", "out_file"]


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
    parser_run.add_argument("-n", "--np", type=int, help="Number of processes to run")
    parser_run.add_argument("-i", "--fi", help="Path to input file to read")
    parser_run.add_argument("--clean", action="store_true", help="Cleans up files from a previous run")

    parser_average = subparsers.add_parser("average",
        description="average", help="averages executions times for a particular configuration of SNAP run")    
    parser_average.add_argument("-i", help="Name of input file")
    parser_average.add_argument("-n", "--np", type=int, help="Number of processes")
    parser_average.add_argument("--cpu", action="store_true", help="Avg cpu builds")
    parser_average.add_argument("--gpu", action="store_true", help="Avg gpu builds")
    parser_average.add_argument("-b", "--blocks", type=int, help="Number of blocks")
    parser_average.add_argument("-t", "--threads", type=int, help="Number of threads")

    parser_compare = subparsers.add_parser("compare",
        description="compare", help="compares two given files for likeness")
    parser_compare.add_argument("-f1", "--file1", help="name of output file 1")
    parser_compare.add_argument("-f2", "--file2", help="name of output file 2")

    return parent_parser

def setup(args):
    subprocess.call([".", "~/intel/bin/compilervars.sh", "-arch", "intel64", "-platform", "linux"], shell=True)
    subprocess.call([".", "~/intel/mkl/bin/mklvars.sh", "intel64"], shell=True)
    subprocess.call([".", "~/intel/bin/iccvars.sh", "-arch", "intel64", "-platform", "linux"], shell=True)
    subprocess.call([".", "~/intel/compilers_and_libraries/linux/mpi/intel64/bin/mpivars.sh"], shell=True)

def clean(args):
    delete_files(os.getcwd(), gmon_regex)

def delete_files(dir, pattern):
    for f in find_files_pattern(dir, pattern):
        os.remove(os.path.join(dir, f))

def find_files_pattern(dir, pattern):
    return [f for f in os.listdir(dir) if re.match(pattern, f)]

def run(args, parser):
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
    if args.clean:
        delete_files(os.getcwd(), gmon_regex)

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
    matching = count_matching_regex_in_folder(OUT_DIR, out_name + "_\d$")
    run_id = 1 if matching == 0 else matching + 1
    out_name = "{}_{}".format(out_name, run_id)

    printv("Output file: {}".format(out_name))
    out_path = os.path.join(OUT_DIR, out_name)

    if not os.path.isdir(OUT_DIR):
        os.mkdir(OUT_DIR)

    # print(args.np, type(args.np))
    # print(app, type(app))
    # print(args.fi, type(args.fi))
    # print(out_path, type(out_path))

    start_time = time.time()
    returncode = subprocess.call(["mpirun", "-np", str(args.np), app, "--fi", args.fi, "--fo", out_path])
    end_time = time.time()
    if not returncode == 0:
        print("ERROR: subprocess.call returned {}. See standard output for further information.".format(returncode))
        if os.path.isfile(out_path):
            os.remove(out_path)
        sys.exit(returncode)

    printv("Elapsed time: {} seconds".format(end_time - start_time))

    # Times
    write_headers = False
    if not os.path.isfile(test_data_path):
        write_headers = True

    with open(test_data_path, "a") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=schema)

        if write_headers:
            writer.writeheader()

        c = 0 if args.gpu else 1
        g = 0 if args.cpu else 1
        data = [in_path, run_id, args.np, c, g, 0, 0, end_time - start_time, out_path]
        writer.writerow(dict(zip(schema, data)))

    # Post-processing
    remove_line_from_file(out_path, 4)
    split_file_n_rows_from_end(out_path, ".timing", 22)

def compare(out_file_1, out_file_2):
    f1s = []
    f2s = []
    with open(os.path.join(OUT_DIR, out_file_1), "r") as f1:
        for r1 in f1: f1s.append(r1)
    with open(os.path.join(OUT_DIR, out_file_2), "r") as f2:
        for r2 in f2: f2s.append(r2)

    diffs = 0
    for i, r in enumerate(f1s):
        if not r == f2s[i]:
            diffs += 1
    
    difference = (diffs / len(f1s)) * 100
    print("Percentage difference: {}".format(difference))
    if difference < 0.005:
        print("Outputs are identical")
    else:
        print("Outputs differ")

def average(args):
    average_time_for(args.i, args.np, args.cpu, args.gpu, args.blocks, args.threads)

def average_time_for(input_file, proc, cpu, gpu, blocks, threads):
    if ((cpu is False) and (gpu is False)) or ((cpu is True) and (gpu is True)):
        print("ERROR: Cannot avg. Must set one of either cpu and gpu, not both nor neither")
        sys.exit(1)

    data = []
    with open("testing.csv", "r") as f:
        reader = csv.DictReader(f, fieldnames=schema)
        headers = next(reader)
        data = list(reader)
    
    data = list(filter(lambda x: filter_row(x, input_file, proc, cpu, gpu, blocks, threads), data))
    if len(data) == 0:
        print("ERROR: No files found matching that configuration.")
        sys.exit(1)

    for dct in data:
        printv(dct)
    
    printv("Averaging...")
    print(reduce(lambda x, y: x + y, [float(d["time"]) for d in data]) / len(data))
    
def filter_row(row, input_file, proc, cpu, gpu, blocks, threads):
    c = int(row["cpu"])
    g = int(row["gpu"])
    
    return (row["input_file"] == input_file and 
        int(row["np"]) == proc and 
        bool(c) == cpu and 
        bool(g) == gpu and
        int(row["block_num"]) == blocks and
        int(row["tpb"]) == threads)

def avg_times():
    test_data = []
    sets = dict()

    with open(test_data_path, "r") as data:
        reader = csv.DictReader(data, fieldnames=schema)
        headers = next(reader)
        for row in reader:
            test_data.append(row)    

    for d in test_data:
        in_file = d["input_file"]
        if in_file in sets:
            sets[in_file].append(d)
        else:
            sets[in_file] = [d]

    for f in sets:
        runs = sets[f]
        cpu_runs
        avg_time = reduce(lambda x, y: x + y, map(lambda x: float(x["time"]), runs)) / len(runs)
        print(f, avg_time)

# Inspired by https://stackoverflow.com/questions/4710067/using-python-for-deleting-a-specific-line-in-a-file
# lines is 1-indexed
def remove_line_from_file(file, *lines):
    with open(file, "r+") as f:
        new_f = f.readlines()
        f.seek(0)
        for i, row in enumerate(new_f):
            if (i + 1) not in lines:
                f.write(row)
        f.truncate()

def split_file_n_rows_from_end(file, new_file_suffix, n):
    with open(file, "r+") as f:
        with open(file + new_file_suffix, "w") as nf:
            new_f = f.readlines()
            f.seek(0)
            for i, row in enumerate(new_f):
                if (i + 1) < len(new_f) - n:
                    f.write(row)
                else:
                    nf.write(row)
            f.truncate()

if __name__ == "__main__":
    parser = parse_arguments()
    args = parser.parse_args()
    if args.verbose: VERBOSE = True

    printv("Command:", args.command)

    if args.command == None:
        parser.print_help()
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

    elif args.command == "average":
        printv("executing average")
        average(args)
        sys.exit(0)
    
    elif args.command == "compare":
        printv("executing compare")
        compare(args.file1, args.file2)
        sys.exit(0)
