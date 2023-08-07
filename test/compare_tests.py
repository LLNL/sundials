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

    parser.add_argument('--release', dest='release', action='store_true', help='indicate if the current run to process is a release')

    parser.add_argument('--calidir', dest='caliDir', type=str, help='path to directory containing caliper files', default="/usr/workspace/sundials/califiles")

    parser.add_argument('--releasedir', dest='releaseDir', type=str, help='path to directory containing release caliper files', default="/usr/workspace/sundials/califiles/Release")

    parser.add_argument('--outpath', dest='outPath', type=str, help='path to directory to write results to', default="/g/g20/pan13/shrimp")

    args = parser.parse_args()

    release = args.release
    releaseDir = args.releaseDir
    caliDir = args.caliDir
    outPath = args.outPath

    # Get the latest test run
    runDirs = glob.glob("%s/Testing/*" % caliDir, recursive = True)
    runDirs.sort(key=os.path.getmtime, reverse=True)
    runDir = runDirs[0]

    runFile = glob.glob(runDir)[0]
    # Load in file to access release value
    th_temp = tt.Thicket.from_caliperreader(runFile) 
    if release:
        # Compare against the last release
        versionDirs = glob.glob("%s/*" % releaseDir)
        versionDirs.sort(key=os.path.getmtime, reverse=True)
        versionDir = versionDirs[1]
    else:
        # Compare against latest release
        version = th_temp.metadata['sundials_version'].values[0]
        versionDir = "%s/%s" % (releaseDir, version)

    # Gather files to process
    runFiles = glob.glob("%s/*.cali" % (runDir))

    if not os.path.exists(outPath):
        os.makedirs(outPath)
    outFile = open("%s/output.out" % outPath, 'w')

    # Compare test results against past runs. If a test performs below a threshold, output test name to outFile.
    with mp.Pool() as pool:
        for res in pool.starmap(compare_against_release, [(versionDir, i) for i in runFiles]):
            if res:
                outFile.write(res + "\n")
    outFile.close()

    outFile = open("%s/output.out" % outPath, 'r')
    try:
        outLines = outFile.readlines()
    finally:
        outFile.close()

    if (len(outLines) == 0):
        return -1 
    return 0

def compare_against_release(releaseDir, file):
    th = tt.Thicket.from_caliperreader(file)

    testName = th.metadata['env.TEST_NAME'].values[0]

    # Gather release run
    release_file = glob.glob("%s/Testing/*/%s.*.cali" % (releaseDir, testName), recursive=True)
    release_th = tt.Thicket.from_caliperreader(release_file)

    metrics = ['Max time/rank']
    tt.mean(release_th, columns=metrics)
    tt.mean(th, columns=metrics)

    ratio = th.statsframe.dataframe['Max time/rank_mean'] / release_th.statsframe.dataframe['Max time/rank_mean']
    print(ratio[0])
    tolerance = 0.0002
    if 1 - ratio[0] < tolerance:
        return testName

if __name__ == "__main__":
    main()
