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
	svnRepo = "file:///usr/casc/sundials/svnrepo/trunk"
	sunBaseDir = "/usr/casc/sundials/devtest"
	sunBinDir = os.path.join(sunBaseDir, "bin")
	sunNightlyDir = os.path.join(sunBaseDir, "nightly")
	logFileName = "build.log"
	tmpLogFile = os.path.join(sunNightlyDir, logFileName)
	sys.stdout = Logger(tmpLogFile)
	
	print "SUNDIALS Automated Build/Test"
	startTime = datetime.datetime.now()
	print startTime.ctime()

	# parse command line options
	parser = OptionParser()
	parser.add_option("-r", "--rev", dest="rev", help="svn revision to checkout", metavar="number", default=0) # 0 -> latest

	(options, args) = parser.parse_args()
	
	
	# check for required options (none required at this time)

	try:	
		try:
			# determine revision to checkout
			svnrev = options.rev
			if (svnrev == 0):
				# get revision number from "Revision :  nnnn"
				cmd = "svn info " + svnRepo + " | grep Revision"
				svnrev = runCommandPopen(cmd).replace(' ','').split(':')[1].strip()
			#svnrev = 'rev' + svnrev
		
			print "svn revision number: " + svnrev
		
			# create directory for checkout
			sunCheckoutDir = os.path.join(sunNightlyDir, "rev"+svnrev)
			print "Creating checkout directory: " + sunCheckoutDir
			if not os.path.exists(sunCheckoutDir):
				os.makedirs(sunCheckoutDir)
			
			# if directory is not empty - exit
			if os.listdir(sunCheckoutDir) != []:
				msg = "Checkout revision: '" + sunCheckoutDir + "' already present."
				print msg
				raise msg
				
			# checkout desired revision
			sunSrcDir = os.path.join(sunCheckoutDir,"trunk")
			cmd = "svn checkout " + svnRepo + " " + sunSrcDir
			print "running: " + cmd + "..."
			cmdout = runCommandPopen(cmd)
			print cmdout
			
			# create build directory
			sunBuildDir = os.path.join(sunCheckoutDir, "build")
			print "Creating build directory: " + sunBuildDir
			if not os.path.exists(sunBuildDir):
				os.makedirs(sunBuildDir)
			
			# if directory is not empty - exit
			if os.listdir(sunBuildDir) != []:
				msg = "Error: build directory: '" + sunBuildDir + "' is not empty."
				print msg
				raise msg
				
			# change working directory to build dir
			os.chdir(sunBuildDir)
			
			# run CMake to configure
			cmd = "cmake -DMPI_ENABLE=ON " + sunSrcDir
			print "Configuring with:  " + cmd + " ..."
			cmdout = runCommandPopen(cmd)
			print cmdout
			
			# run make to build libs and executables
			cmd = "make"
			print "Building with:  " + cmd + " ..."
			cmdout = runCommandPopen(cmd)
			print cmdout
			
			# run 'make test' to test examples
			cmd = "make test"
			print "Testing with:  " + cmd + " ..."
			cmdout = runTestCommandPopen(cmd)
			print cmdout
			
			# cleanup nightly build directory - keep last 7
			
		
		# catch any exception
		except Exception, e:
			print "Execption: ", e
			# send log file in progress...
			print "Build log file: ", tmpLogFile
			sys.stdout.flush()
			sys.stdout.close()
			sendEmail(tmpLogFile)
			return 2

	# catch any exception
	except:
		print "Build Aborted!"
		# send log file in progress...
		#print "Sending email from file: ", tmpLogFile
		sys.stdout.flush()
		sys.stdout.close()
		sendEmail(tmpLogFile, "SUNDIALS build aborted for rev: " + svnrev)
		raise "ABORT"
	
	# Build Done - move and mail log file
	endTime = datetime.datetime.now()
	print "End time: ", endTime.ctime()
	elapsedTime = endTime - startTime
	print "Elapsed Time:", elapsedTime
	
	sys.stdout.flush()
	sys.stdout.close()

	# move the log file to the revision checkout directory
	revLogFile = os.path.join(sunCheckoutDir, logFileName)
	os.rename(tmpLogFile, revLogFile)
	
	# send email
	sendEmail(revLogFile, "New SUNDIALS build for svn rev: " + svnrev)

	# nothing else to do
	return 0

############################
# functions
############################
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
	p = subprocess.Popen(cmd, shell=True, \
	stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	(stdout, stderr) = p.communicate()
	errcode = p.returncode
	if p.returncode != 0:
		print ("Error: Command did not complete successfully. (error code: %s)" % p.returncode)
		print "Error/failure details:\n", stderr
		raise Exception(p.returncode)
	
	# return last line of output
	return ''.join(stdout.splitlines()[-1:])

#
# run external command 
#
def runTestCommandPopen(cmd):
	p = subprocess.Popen(cmd, shell=True, \
	stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	(stdout, stderr) = p.communicate()
	errcode = p.returncode
	if p.returncode != 0:
		print ("Error: Command did not complete successfully. (error code: %s)" % p.returncode)
		#print "Error/failure details:\n", stderr
		#raise Exception(p.returncode)
	
	# return last line of output  
	# TODO: capture only 'meaningful' last lines (may need to add a 'Completion Line' in testRunner to look for)
	return '\n'.join(stdout.splitlines()[-12:])


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
	you = "banks12@llnl.gov"
	#you = "sundials-devs@llnl.gov"
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


