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
    sunBaseDir = "/usr/casc/sundials/devtest"
    sunBinDir = os.path.join(sunBaseDir, "bin")
    sunNightlyDir = os.path.join(sunBaseDir, "nightly")
    todayDate = time.strftime("%Y-%m-%d")
    logFileName = "build.log"
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
            cmd = "cmake -DMPI_ENABLE=ON " + sunSrcDir
            print "\n*** Configuring with:  " + cmd + " ..."
            cmdout = runCommandPopen(cmd)
            print cmdout
    			
            # run make to build libs and executables
            cmd = "make"
            print "\n*** Building with:  " + cmd + " ..."
            cmdout = runCommandPopen(cmd)
            print cmdout
    
            # run 'make test' to test examples
            cmd = "make test"
            print "\n*** Testing with:  " + cmd + " ..."
            (msg, cmdout) = runTestCommandPopen(cmd)
            print cmdout
    			
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
            cleanup(msg, startTime, sunCheckoutDir, tmpLogFile, logFileName)
            return 2
    
        # build successful - cleanup
        msg = "Success: Rev: " + svnrev + " | " + msg
        cleanup(msg, startTime, sunCheckoutDir, tmpLogFile, logFileName)
    
    # catch any other exception
    except:
        # send log file in progress...
        #print "Sending email from file: ", tmpLogFile
        msg = "ABORT: Build failed for rev: " + str(svnrev)
        cleanup(msg, startTime, sunCheckoutDir, tmpLogFile, logFileName)
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
def cleanup(msg, startTime, sunCheckoutDir, tmpLogFile, logFileName):

    # file spec of final log file in revision checkout directory
    revLogFile = os.path.join(sunCheckoutDir, logFileName)
    
    # print closing info
    print "\n" + msg
    print "For details see:  " + revLogFile
    
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
    cmdout = subprocess.check_output(cmd, shell=True)
    cmdoutList = cmdout.splitlines()
    # return 3rd from last line (% passed line), and full output.
    return (cmdoutList[-3], '\n'.join(cmdoutList))


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

	# me == the sender's email address
	me = "SUNDIALS.sunbuild@llnl.gov"
	# you == the recipient's email address
	you = "sundials-devs@llnl.gov"
	#you = "banks12@llnl.gov" # Eddy
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


