#!/usr/bin/env python
# -------------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL 
# -------------------------------------------------------------------------------
# LLNS Copyright Start
# Copyright (c) 2014, Lawrence Livermore National Security
# This work was performed under the auspices of the U.S. Department 
# of Energy by Lawrence Livermore National Laboratory in part under 
# Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# LLNS Copyright End
# -------------------------------------------------------------------------------
# Send email notification if a SUNDIALS regression test status
# -------------------------------------------------------------------------------

def main():

    import argparse
    import sys, os

    parser = argparse.ArgumentParser(
        description='Send email notification based on regression test status',
        formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument('teststatus',type=str,
                        choices=['success','failure','unstable'],
                        help='Status of regression test')

    parser.add_argument('branchname',type=str,
                        help='Name of branch tested')

    parser.add_argument('testurl',type=str,
                        help='URL for viewing test results')

    # parse command line args 
    args = parser.parse_args()

    # name of regression test log file
    logfile = "suntest.log"

    # if log file exists add url, otherwise create log file
    if (os.path.isfile(logfile)):
        missinglogfile = False
        with open(logfile,"a") as log:
            log.write("View test output at:\n")
            log.write(args.testurl)
    else:
        print "Warning: log file not found"
        missinglogfile = True

        with open(logfile,"w") as log:
            log.write("Warning: log file not found\n")
            log.write("View test output at:\n")
            log.write(args.testurl)

    # determine notification recipient
    if ((args.branchname == 'master') or (args.branchname == 'develop')):
        # SUNDIALS developer list
        # recipient = "sundials-devs@llnl.gov"
        recipient = "gardner48@llnl.gov" # for testing
    else:
        # author of most recent commit
        cmd = "git log --format='%ae' -1"
        recipient = runCommand(cmd)

    # send notification if regression tests fail or log file not found
    if (args.teststatus != 'success'):

        subject = "FAILED: SUNDIALS "+args.branchname+" failed regression tests"
        print "Tests failed, sending notification to",recipient
        sendEmail(recipient, subject, logfile)

    elif (missinglogfile):

        subject = "ERROR: SUNDIALS "+args.branchname+" log file not found"
        print "Log file not found, sending notification to",recipient
        sendEmail(recipient, subject, logfile)

    else:
        print "Tests passed, no notifications sent"

#
# run external command 
#
def runCommand(cmd):

    import subprocess

    cmdout = subprocess.check_output(cmd, shell=True)

    return(cmdout)


#
# send email
#
def sendEmail(recipient, subject, message):

    import smtplib
    from email.MIMEText import MIMEText

    # Open a plain text file for reading. Assumed that
    # the text file contains only ASCII characters.
    fp = open(message, 'rb')
    
    # Create a text/plain message
    msg = MIMEText(fp.read())
    fp.close()

    # sender's email address
    sender = "SUNDIALS.suntest@llnl.gov"

    # email settings
    msg['Subject'] = subject
    msg['From'] = sender
    msg['To'] = recipient

    # Send the message via our own SMTP server, but don't include the
    # envelope header.
    s = smtplib.SMTP('smtp.llnl.gov')
    s.sendmail(sender, [recipient], msg.as_string())
    s.quit()


#
# just run the main routine
#
if __name__ == '__main__':
    main()
