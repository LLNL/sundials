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

    parser.add_argument('--calidir', dest='caliDir', type=str,help='name of folder containing caliper files', default="/usr/workspace/pan13/shrimp2")

    parser.add_argument('--outdir', dest='outDir', type=str, help='path to directory to write results to', default="/g/g20/pan13/shrimp")

    args = parser.parse_args()

    caliDir = args.caliDir
    outPath = args.outpath

    # get the latest test run
    runDirs = glob.glob("%s/Testing/*" % caliDir, recursive = True)
    runDirs.sort(key=os.path.getmtime, reverse=True)
    runDir = runDirs[0]

    # Gather files to process
    testFiles = glob.glob("%s/*.cali" % (runDir))

    if not os.path.exists(outPath):
        os.mkdirs(outPath)
    outFile = open("%s/output.out" % outPath, 'w')
    

    # Compare test results against past runs. If a test performs below a threshold, output test name to outFile.
    with mp.Pool() as pool:
        for res in pool.starmap(process_test, [(caliDir, i) for i in testFiles]):
            if res:
                outFile.write(res + "\n")
    outFile.close()

    # for file in testFiles:
    #     process_test(caliDir, file)

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

    # Gather the last numRuns of tests and remove the current run from set of runs to aggregate
    numRuns = 5
    files = glob.glob("%s/Testing/*/%s.*.cali" % (caliDir, testName), recursive=True)
    files.sort(key=os.path.getmtime, reverse=True)
    files.pop(0) # remove current run
    files = files[:numRuns]

    files.sort(key=os.path.ge)

    th_files = tt.Thicket.from_caliperreader(files)
    metrics = ['Max time/rank']
    tt.mean(th_files, columns=metrics)
    tt.mean(th, columns=metrics)

    ratio = th.statsframe.dataframe['Max time/rank_mean'] / th_files.statsframe.dataframe['Max time/rank_mean']

    if 1 - ratio[0] < 0.1:
        return testName

if __name__ == "__main__":
    main()
