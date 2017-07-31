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
# Send email notification if a SUNDIALS regression test fails
# -------------------------------------------------------------------------------

def main():

    import argparse
    import sys, os

    parser = argparse.ArgumentParser(
        description='Send email notification if regression test fails', 
        formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument('teststatus',type=int,
                        help='Status of regression tests, 0 = success')

    parser.add_argument('branchname',type=str,
                        help='Name of branch tested')

    parser.add_argument('testurl',type=str,
                        help='URL for viewing test results')

    # parse command line args 
    args = parser.parse_args()

    # add URL for test results to log file
    logfile = "suntest.log"
    with open(logfile, "a") as log:
        log.write("View output log at:\n")
        log.write(args.testurl)

    # send notification if any regression tests fail
    if (args.teststatus != 0):

        subject = "FAILED: SUNDIALS "+args.branchname+" branch failed regression tests"

        if ((args.branchname == 'master') or (args.branchname == 'develop')):      
            # SUNDIALS developer list
            recipient = "sundials-devs@llnl.gov"
        else:
            # author of last commit
            cmd = "git log --format='%ae' -1"
            recipient = runCommand(cmd)

        print "Tests failed, sending notification to",recipient

        sendEmail(recipient, subject, logfile)

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
