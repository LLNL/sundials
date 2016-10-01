#!/usr/bin/env python

import sys, os, shutil
import subprocess
import datetime, time
import smtplib
from email.MIMEText import MIMEText
from optparse import OptionParser


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

def main():
    maxRevDirs = 5                  # maximum revision directories to keep
    maxBuildsPerPreviousRev = 3     # maximum number of builds to keep per previous revisions
    maxBuildsPerCurrentRev = 10     # maximum number of builds to keep for current revision
    
    svnRepo = "file:///usr/casc/sundials/svnrepo/trunk"
    #sunBaseDir = "/usr/casc/sundials/etest/devtest"
    sunBaseDir = "/usr/casc/sundials/devtest"
    sunBinDir = os.path.join(sunBaseDir, "bin")
    sunNightlyDir = os.path.join(sunBaseDir, "nightly")
    todayDate = time.strftime("%Y-%m-%d")
    logFileName = "build.log"
    cmakeLogFileName = "cmake.log"
    makeLogFileName = "make.log"
    makeLogFile = "Oops - didn't get this far!"
    cmakeLogFile = "Oops - didn't get this far!"
    tmpLogFile = os.path.join(sunNightlyDir, logFileName)
    sys.stdout = Logger(tmpLogFile)
    
    print "*** SUNDIALS Automated Build/Test ***"
    startTime = datetime.datetime.now()
    print startTime.ctime()

    # parse command line options
    parser = OptionParser()
    parser.add_option("-r", "--rev", dest="rev", help="svn revision to checkout", metavar="number", default=0) # 0 -> latest

    (options, args) = parser.parse_args()
    
    # check for required options (none required at this time...could easily take specific revision# to build)
    try:
        try:    
            # determine revision to checkout
            svnrev = options.rev
            if (svnrev == 0):
                # get revision number from "Revision :  nnnn"
                cmd = "svn info " + svnRepo + " | grep Revision"
                svnrev = runCommandPopen(cmd).replace(' ','').split(':')[1].strip()
                print svnrev
            print "svn revision number: " + svnrev
            
            # create directory for "today's" checkout
            sunRevDir = os.path.join(sunNightlyDir, "rev"+svnrev)
            sunCheckoutDir = os.path.join(sunRevDir, todayDate)
            print "\n*** Creating checkout directory: " + sunCheckoutDir + " ..."
            if not os.path.exists(sunCheckoutDir):
                os.makedirs(sunCheckoutDir)
    
            # if directory is not empty - exit
            if os.listdir(sunCheckoutDir) != []:
                msg = "Revision 'rev"+svnrev+"/"+todayDate + "' already present."
                print msg
                raise Exception(msg)
    
            # checkout desired revision
            sunSrcDir = os.path.join(sunCheckoutDir,"trunk")
            cmd = "svn checkout -r " + svnrev + " " + svnRepo + " " + sunSrcDir
            print "\n*** Checkout revision : " + cmd + " ..."
            cmdout = runCommandPopen(cmd)
            print cmdout
                
            # create build directory
            sunBuildDir = os.path.join(sunCheckoutDir, "build")
            print "\n*** Creating build directory: " + sunBuildDir
            if not os.path.exists(sunBuildDir):
                os.makedirs(sunBuildDir)
    
            # if directory is not empty - exit
            if os.listdir(sunBuildDir) != []:
                msg = "Build directory for: 'rev" + svnrev + "/" + todayDate + "' is not empty."
                print msg
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
            cmd = cmd + "-DMPI_ENABLE=ON \ \n"
            # enable FCMIX
            cmd = cmd + "-DFCMIX_ENABLE=TRUE \ \n"
            # enable F90
            cmd = cmd + "-DF90_ENABLE=TRUE \ \n"
            # enable C++
            cmd = cmd + "-DCXX_ENABLE=TRUE \ \n"
            # enable lapack   (NOTE: will find libraries in LD_LIBRARY_PATH)
            cmd = cmd + "-DLAPACK_ENABLE=ON \ \n"
            # enable klu
            cmd = cmd + "-DKLU_ENABLE=ON \ \n"
            cmd = cmd + "-DKLU_INCLUDE_DIR=/usr/casc/sundials/apps/rh6/suitesparse/4.5.3/include \ \n"
            cmd = cmd + "-DKLU_LIBRARY_DIR=/usr/casc/sundials/apps/rh6/suitesparse/4.5.3/lib \ \n"
            # enable hypre
            cmd = cmd + "-DHYPRE_ENABLE=ON \ \n"
            cmd = cmd + "-DHYPRE_INCLUDE_DIR=/usr/casc/sundials/apps/rh6/hypre/2.11.1/include \ \n"
            cmd = cmd + "-DHYPRE_LIBRARY=/usr/casc/sundials/apps/rh6/hypre/2.11.1/lib/libHYPRE.a \ \n"
            # enable PETSc
            cmd = cmd + "-DPETSC_ENABLE=ON \ \n"
            cmd = cmd + "-DPETSC_INCLUDE_DIR=/usr/casc/sundials/apps/rh6/petsc/3.7.2/include \ \n"
            cmd = cmd + "-DPETSC_LIBRARY_DIR=/usr/casc/sundials/apps/rh6/petsc/3.7.2/lib \ \n"
            # enable openmp
            cmd = cmd + "-DOPENMP_ENABLE=ON \ \n"
            # enable pthreads
            cmd = cmd + "-DPTHREAD_ENABLE=ON \ \n"
            # enable SUPERLU_MT
            cmd = cmd + "-DSUPERLUMT_ENABLE=ON \ \n"
            cmd = cmd + "-DSUPERLUMT_INCLUDE_DIR=/usr/casc/sundials/apps/rh6/superlu_mt/SuperLU_MT_3.1/SRC \ \n"
            cmd = cmd + "-DSUPERLUMT_LIBRARY_DIR=/usr/casc/sundials/apps/rh6/superlu_mt/SuperLU_MT_3.1/lib \ \n"
            cmd = cmd + "-DSUPERLUMT_THREAD_TYPE=Pthread \ \n"
            # specify source
            cmd = cmd + sunSrcDir
            # redirect output to config log file
            cmd = cmd + " &> " + cmakeLogFile
            
            print "\n*** Running CMake:\n" + cmd + " ..."
            
            # remove newlines from cmd
            cmd = cmd.replace('\n','')
            cmd = cmd.replace('\\','')
            #print "\n*** FYI: " + cmd
            cmdout = runCommandPopen(cmd)
            print cmdout
            
            # run make to build libs and executables
            makeLogFile = os.path.join(sunBuildDir, makeLogFileName)
            cmd = "make -j4 &> " + makeLogFile
            print "\n*** Building with:  " + cmd + " ..."
            cmdout = runCommandPopen(cmd)
            print cmdout
    
            # run 'make test' to test examples
            cmd = "make test"
            print "\n*** Testing with:  " + cmd + " ..."
            (msg, cmdout) = runTestCommandPopen(cmd)
            print cmdout
            
            # Parse cmake log file for warning
                
            # purge nightly build directory (keep maxRevDirs number of revisions)
            print "\n*** Purging nightly build directory..."
            revDirList = purgeNightlyDir(sunNightlyDir, maxRevDirs)
            
            # now purge each revision directory (except most current)
            for revDir in revDirList[:-1]:
                purgeNightlyDir(revDir, maxBuildsPerPreviousRev)
            
            # now purge current revision directory
            purgeNightlyDir(revDirList[-1], maxBuildsPerCurrentRev) # most current rev is last on list
            
        except Exception, e:
            msg = "FAILED: " + str(e)
            cleanup(msg, startTime, sunCheckoutDir, tmpLogFile, logFileName, cmakeLogFile, makeLogFile)
            return 2
    
        # build successful - cleanup
        msg = "Finished: Rev: " + svnrev + " | " + msg
        cleanup(msg, startTime, sunCheckoutDir, tmpLogFile, logFileName, cmakeLogFile, makeLogFile)
    
    # catch any other exception
    except:
        # send log file in progress...
        #print "Sending email from file: ", tmpLogFile
        msg = "ABORT: Build failed for rev: " + str(svnrev)
        cleanup(msg, startTime, sunCheckoutDir, tmpLogFile, logFileName, cmakeLogFile, makeLogFile)
        raise Exception("ABORT")
    
    # nothing else to do
    return 0

############################
# functions
############################
#
# purge the nightly directory of older builds
#
def purgeNightlyDir(rootDir, maxKeepDirs):
    dirList = []
    for tmpName in sorted(os.listdir(rootDir)):
        tmpSpec = os.path.join(rootDir, tmpName)
        if os.path.isdir(tmpSpec):
            dirExtension = os.path.splitext(tmpSpec)[1]
            # only purge directories without extensions (e.g. allows user to add '.save' to specific build) 
            if dirExtension == '':
                dirList.append(tmpSpec)
                #print tmpSpec
    
    # if more than maxRevDirs - remove them
    numDirs = len(dirList)
    removeDirCount = numDirs - maxKeepDirs
    for  i in xrange(removeDirCount):
        shutil.rmtree(dirList[i])
    
    # return list of remaining directories
    return dirList[removeDirCount:numDirs]

#
# cleanup
#
def cleanup(msg, startTime, sunCheckoutDir, tmpLogFile, logFileName, cmakeLogFile, makeLogFile):

    # file spec of final log file in revision checkout directory
    revLogFile = os.path.join(sunCheckoutDir, logFileName)
    
    # print closing info
    print "\n" + msg
    print "For build script details see:    " + revLogFile
    print "For CMake specific details see:  " + cmakeLogFile
    print "For make specific details see:   " + makeLogFile
    print "For test details see: " + sunCheckoutDir + "/build/Testing/Temporary/LastTest.log"
    
    endTime = datetime.datetime.now()
    print "\nEnd time: ", endTime.ctime()
    elapsedTime = endTime - startTime
    print "Elapsed Time:", elapsedTime

    sys.stdout.flush()
    sys.stdout.close()

    # move log file, don't overwrite a previous log file (since this is daily - this should never happen. But...)
    if fileExists(revLogFile):
        revLogFile = tmpLogFile
    else:
        os.rename(tmpLogFile, revLogFile)

    # send email
    sendEmail(revLogFile, msg)
    

    return


#
# check if file exists
#
def fileExists(fname):
    if os.path.exists(fname):
        return True
    return False

#
# run external command 
#
def runCommandPopen(cmd):
    cmdout = subprocess.check_output(cmd, shell=True)
    # return last line of output
    return ''.join(cmdout.splitlines()[-1:])
#
# run external command 
#
def runTestCommandPopen(cmd):
    cmdout = ""
    try: 
        cmdout = subprocess.check_output(cmd, shell=True)
    except subprocess.CalledProcessError, e:
        cmdout = e.output
    finally:
        percentLine = findPercentSuccessfulLine(cmdout)
        return (percentLine, cmdout)

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
    #you = "banks12@llnl.gov" # Eddy for testing
    msg['Subject'] = subject
    msg['From'] = me
    msg['To'] = you

    # Send the message via our own SMTP server, but don't include the
    # envelope header.
    s = smtplib.SMTP('smtp.llnl.gov')
    s.sendmail(me, [you], msg.as_string())
    s.quit()

############################

if __name__ == "__main__":
    sys.exit(main())

