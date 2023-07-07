import os
import subprocess
import sys
import re
import glob
import argparse
import multiprocessing as mp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import hatchet as ht
import thicket as tt

def main():
    parser = argparse.ArgumentParser(description='Compare Sundials performance results against previous results')

    parser.add_argument('--calidir', dest='caliDir', type=str,help='name of folder containing caliper files', default="/usr/workspace/pan13/shrimp")

    parser.add_argument('--outpath', dest='outpath', type=str, help='path to file to write results to', default="/g/g20/pan13/shrimp/output.out")

    args = parser.parse_args()

    caliDir = args.caliDir
    outPath = args.outpath

    # add a thing to check if the output file exists and make sure it does

    # get the latest test run
    runDirs = glob.glob("%s/Testing/*" % caliDir, recursive = True)
    runDirs.sort(key=os.path.getmtime, reverse=True)
    runDir = runDirs[0]

    # Gather files to process
    testFiles = glob.glob("%s/*.cali" % (runDir))

    # Compare test results against past runs. If a test performs below a threshold, output test name to file.
    outFile = open(outPath, 'w')
    with mp.Pool() as pool:
        for res in pool.starmap(process_test, [(caliDir, i) for i in testFiles]):
            if res:
                outFile.write(res + "\n")
    outFile.close()

    outFile = open(outPath, 'r')
    try:
        outLines = outFile.readlines()
    finally:
        outFile.close()

    if (len(outLines) == 0):
        return -1 
    return 0

def process_test(caliDir, file):

    th = tt.Thicket.from_caliperreader(file)

    testName = th.metadata['env.TEST_NAME'].values[0]
    files = glob.glob("%s/Testing/*/%s.*.cali" % (caliDir, testName), recursive=True)
    th_files = tt.Thicket.from_caliperreader(files)

    metrics = ['Max time/rank']
    tt.mean(th_files, columns=metrics)
    tt.mean(th, columns=metrics)

    ratio = th.statsframe.dataframe['Max time/rank_mean'] / th_files.statsframe.dataframe['Max time/rank_mean']

    if 1 - ratio[0] < 0.1:
        return testName

if __name__ == "__main__":
    main()
