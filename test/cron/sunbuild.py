#! /usr/bin/env python3
# -------------------------------------------------------------------------------
# Programmer(s): Eddy Banks and David J. Gardner @ LLNL
# -------------------------------------------------------------------------------
# This is a python script for running nightly tests on the most recent commit to
# the SUNDIALS git repository. Results are stored in /usr/casc/sundials/devtest
# and the main log file is sent to the SUNDIALS developers mailing list
# "sundials-devs@llnl.gov" after all tests are completed.
# -------------------------------------------------------------------------------

import sys, os, shutil
import subprocess
import datetime, time
import smtplib
from email.MIMEText import MIMEText
from optparse import OptionParser

def main():

    maxCommitDirs = 5               # max number of commit directories to keep
    maxBuildsPerPreviousCommit = 3  # max number of builds to keep per previous commit
    maxBuildsPerCurrentCommit  = 10 # max number of builds to keep for current commit

    sunGitRepo = "ssh://git@mystash.llnl.gov:7999/sundials/sunrepo.git"

    sunBaseDir = "/usr/casc/sundials/devtest"
    sunNightlyDir = os.path.join(sunBaseDir,"nightly")

    # main log file for this script
    buildLogFileName = "build.log"
    buildLogFile = os.path.join(sunNightlyDir, buildLogFileName)
    sys.stdout = Logger(buildLogFile)

    # other log files created in this script
    cmakeLogFileName = "cmake.log"
    makeLogFileName  = "make.log"

    # log file created by the testRunner script
    testLogFileName  = "build/Testing/Temporary/LastTest.log"

    # place holders in case an error is encountered before these variables are created
    makeLogFile  = "Oops - didn't get this far!"
    cmakeLogFile = "Oops - didn't get this far!"
    testLogFile  = "Oops - didn't get this far!"
    sunTestDir   = "Oops - didn't get this far!"
    msg          = "Oops - didn't get this far!"

    todayDate = time.strftime("%Y-%m-%d")

    print "*** SUNDIALS Automated Build/Test ***"
    startTime = datetime.datetime.now()
    print "\nStart time:", startTime.ctime()

    try:
        # get most recent commit hash
        cmd = "git ls-remote " + sunGitRepo + " | grep HEAD"
        cmdout = runCommand(cmd)

        commitHashLong  = cmdout.split()[0].strip()
        commitHashShort = commitHashLong[0:7]
        print "Git commit:",commitHashLong

        # create commit test directory
        sunTestDir = os.path.join(sunNightlyDir,
                                  "commit_"+commitHashShort,
                                  todayDate)
        print "\n*** Creating commit test directory: " + sunTestDir
        if not os.path.exists(sunTestDir):
            os.makedirs(sunTestDir)

        # change working directory to test dir
        os.chdir(sunTestDir)

        # if commit test directory is not empty - exit
        if os.listdir(sunTestDir) != []:
            msg = "'" + sunTestDir + "' is not empty."
            raise Exception(msg)

        # checkout only the most recent commit to the repo
        cmd = "git clone --quiet --depth 1 " + sunGitRepo
        print "\n*** Cloning most recent commit: " + cmd
        cmdout = runCommand(cmd)

        # change working directory to clone of sundials repo
        sunSrcDir = os.path.join(sunTestDir,"sunrepo")
        os.chdir(sunSrcDir)

        # print most recent commit message
        cmd = "git log -1"
        print "\n*** Most recent commit message: " + cmd + " \n"
        cmdout = runCommand(cmd)
        print cmdout

        # change working directory to test dir
        os.chdir(sunTestDir)

        # create build directory
        sunBuildDir = os.path.join(sunTestDir, "build")
        print "\n*** Creating build directory: " + sunBuildDir
        if not os.path.exists(sunBuildDir):
            os.makedirs(sunBuildDir)

        # if build directory is not empty - exit
        if os.listdir(sunBuildDir) != []:
            msg = "'" + sunBuildDir + "' is not empty."
            raise Exception(msg)

        # change working directory to build dir
        os.chdir(sunBuildDir)

        # run CMake to configure
        cmakeLogFile = os.path.join(sunBuildDir, cmakeLogFileName)
        cmd = "cmake \ \n"
        # set compiler flags to check for non-standard code
        # -ansi  OR  -std-c89  OR  -std=c99  OR  -std=c11
        # NOTE: PETSC requires -std=c99 or newer
        cmd = cmd + "-DCMAKE_C_FLAGS='-Wall -std=c99 -pedantic' \ \n"
        # enable mpi
        cmd = cmd + "-DENABLE_MPI=ON \ \n"
        # enable C++
        cmd = cmd + "-DCXX_ENABLE=TRUE \ \n"
        # enable lapack   (NOTE: will find libraries in LD_LIBRARY_PATH)
        cmd = cmd + "-DENABLE_LAPACK=ON \ \n"
        # enable klu
        cmd = cmd + "-DENABLE_KLU=ON \ \n"
        cmd = cmd + "-DKLU_INCLUDE_DIR=/usr/casc/sundials/apps/rh6/suitesparse/4.5.3/include \ \n"
        cmd = cmd + "-DKLU_LIBRARY_DIR=/usr/casc/sundials/apps/rh6/suitesparse/4.5.3/lib \ \n"
        # enable hypre
        cmd = cmd + "-DENABLE_HYPRE=ON \ \n"
        cmd = cmd + "-DHYPRE_INCLUDE_DIR=/usr/casc/sundials/apps/rh6/hypre/2.11.1_long_int_fpic/include \ \n"
        cmd = cmd + "-DHYPRE_LIBRARY_DIR=/usr/casc/sundials/apps/rh6/hypre/2.11.1_long_int_fpic/lib \ \n"
        # enable PETSc
        cmd = cmd + "-DENABLE_PETSC=ON \ \n"
        cmd = cmd + "-DPETSC_INCLUDE_DIR=/usr/casc/sundials/apps/rh6/petsc/3.7.2_long_int/include \ \n"
        cmd = cmd + "-DPETSC_LIBRARY_DIR=/usr/casc/sundials/apps/rh6/petsc/3.7.2_long_int/lib \ \n"
        # enable openmp
        cmd = cmd + "-DENABLE_OPENMP=ON \ \n"
        # enable pthreads
        cmd = cmd + "-DENABLE_PTHREAD=ON \ \n"
        # enable SUPERLU_MT
        cmd = cmd + "-DENABLE_SUPERLUMT=ON \ \n"
        cmd = cmd + "-DSUPERLUMT_INCLUDE_DIR=/usr/casc/sundials/apps/rh6/superlu_mt/SuperLU_MT_3.1_long_int_fpic/SRC \ \n"
        cmd = cmd + "-DSUPERLUMT_LIBRARY_DIR=/usr/casc/sundials/apps/rh6/superlu_mt/SuperLU_MT_3.1_long_int_fpic/lib \ \n"
        cmd = cmd + "-DSUPERLUMT_THREAD_TYPE=Pthread \ \n"
        # turn on development tests
        cmd = cmd + "-DSUNDIALS_TEST_DEVTESTS=ON \ \n"
        # specify source
        cmd = cmd + sunSrcDir
        # redirect output to config log file
        cmd = cmd + " &> " + cmakeLogFile

        print "\n*** Running CMake:\n" + cmd

        # remove newlines from cmd
        cmd = cmd.replace('\n','')
        cmd = cmd.replace('\\','')
        cmdout = runCommand(cmd)
        print cmdout

        # run make to build libs and executables
        makeLogFile = os.path.join(sunBuildDir, makeLogFileName)
        cmd = "make -j4 &> " + makeLogFile
        print "\n*** Building with: " + cmd
        cmdout = runCommand(cmd)
        print cmdout

        # set location of test runner log file
        testLogFile = os.path.join(sunBuildDir, testLogFileName)

        # run 'make test' to test examples
        cmd = "make test"
        print "\n*** Testing with: " + cmd
        cmdout = runCommand(cmd)
        print cmdout

        # find the percentage of successful tests
        percentLine = findPercentSuccessfulLine(cmdout)

        # tests finished
        msg = "Finished: Commit: " + commitHashShort + " | " + percentLine

        # purge nightly build directory (keep maxRevDirs number of revisions)
        print "\n*** Purging nightly build directory"
        commitDirList = purgeNightlyDir(sunNightlyDir, maxCommitDirs)

        # now purge each revision directory (except most current)
        for commitDir in commitDirList[:-1]:
            purgeNightlyDir(commitDir, maxBuildsPerPreviousCommit)

        # now purge current revision directory (last entry in list)
        purgeNightlyDir(commitDirList[-1], maxBuildsPerCurrentCommit)

    except Exception as e:
        msg = "FAILED: Commit: " + commitHashShort + " | " + str(e)
        return 1

    finally:
        # clean up test after success or fail
        cleanup(msg, startTime, sunTestDir, buildLogFile, buildLogFileName,
                cmakeLogFile, makeLogFile, testLogFile)

    return 0

# ===============================================================================

#
# class to redirect system output to log file and screen
#
class Logger(object):
    def __init__(self, logFile):
        self.terminal = sys.stdout
        self.log = open(logFile, "w")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        self.log.flush()

    def close(self):
        self.log.close()

#
# run external command
#
def runCommand(cmd):
    cmdout = subprocess.check_output(cmd, shell=True)

    return(cmdout)

#
# find % test successful line in 'make test' output
#
def findPercentSuccessfulLine(cmdout):
    # search for line with first word containing %
    percentLine = ""
    for line in cmdout.splitlines():
        if '%' in line.split(' ')[0]:
            percentLine = line
            break

    return percentLine

#
# sort directories by modification time
#
def sortDirByTime(rootDir):
    mtime   = lambda x: os.stat(os.path.join(rootDir, x)).st_mtime
    dirList = sorted(os.listdir(rootDir), key=mtime)
    return dirList

#
# purge the nightly directory of older builds
#
def purgeNightlyDir(rootDir, maxKeepDirs):
    dirList = []
    for tmpName in sortDirByTime(rootDir):
        tmpSpec = os.path.join(rootDir, tmpName)
        if os.path.isdir(tmpSpec):
            dirExtension = os.path.splitext(tmpSpec)[1]
            # only purge directories without extensions
            # (e.g. allows user to add '.save' to specific build)
            if dirExtension == '':
                dirList.append(tmpSpec)
                #print tmpSpec

    # if more than maxRevDirs - remove them
    numDirs = len(dirList)
    removeDirCount = numDirs - maxKeepDirs
    for i in xrange(removeDirCount):
        shutil.rmtree(dirList[i])

    # return list of remaining directories
    return dirList[removeDirCount:numDirs]

#
# cleanup after test
#
def cleanup(msg, startTime, sunTestDir, buildLogFile, buildLogFileName,
            cmakeLogFile, makeLogFile, testLogFile):

    # move log file to test directory
    if (os.path.isdir(sunTestDir)):
        finalBuildLogFile = os.path.join(sunTestDir, buildLogFileName)
        os.rename(buildLogFile, finalBuildLogFile)
    else:
        finalBuildLogFile = buildLogFile

    # print closing info
    print "\n" + msg
    print "For build script details see:    " + finalBuildLogFile
    print "For CMake specific details see:  " + cmakeLogFile
    print "For make specific details see:   " + makeLogFile
    print "For test details see:            " + testLogFile

    endTime = datetime.datetime.now()
    print "\nEnd time:", endTime.ctime()
    elapsedTime = endTime - startTime
    print "Elapsed Time:", elapsedTime

    sys.stdout.flush()
    sys.stdout.close()

    # send email
    sendEmail(finalBuildLogFile, msg)

    return

#
# send email to SUNDIALS dev list
#
def sendEmail(logFile, subject):
    # Open a plain text file for reading.  For this example, assume that
    # the text file contains only ASCII characters.
    fp = open(logFile, 'rb')

    # Create a text/plain message
    msg = MIMEText(fp.read())
    fp.close()

    ### me == the sender's email address
    me = "SUNDIALS.sunbuild@llnl.gov"
    ### you == the recipient's email address
    you = "sundials-devs@llnl.gov"
    msg['Subject'] = subject
    msg['From'] = me
    msg['To'] = you

    # Send the message via our own SMTP server, but don't include the
    # envelope header.
    s = smtplib.SMTP('smtp.llnl.gov')
    s.sendmail(me, [you], msg.as_string())
    s.quit()

# ===============================================================================

if __name__ == "__main__":
    sys.exit(main())
