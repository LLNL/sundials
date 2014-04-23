#!/bin/csh


############################################################################
# $Revision$
# $Date$
############################################################################
#
# Filename: sundials_check_bash.sh
# Programmer: Aaron Collier @ LLNL
#
############################################################################


set CHANGE_FILE=$argv[1]
set USE_LOGIN_SHELL=$argv[2]

unalias cp

# look for BASH
set FIND_BASH=`which bash`
# determine local operating system
set OS_TYPE=`uname -s`

# AIX systems do not set PATH correctly when using SSH
if ( "$OS_TYPE" == "AIX" ) then
  if ( "$USE_LOGIN_SHELL" == "yes" ) then
    echo "#\!/usr/local/bin/bash --login" > temp.txt
  else if ( "$USE_LOGIN_SHELL" == "no" ) then
    echo "#\!/usr/local/bin/bash" > temp.txt
  endif
# other systems can use just search for BASH
else if ( "$OS_TYPE" != "AIX" ) then
  if ( "$USE_LOGIN_SHELL" == "yes" ) then
    echo "#\!$FIND_BASH --login" > temp.txt
  else if ( "$USE_LOGIN_SHELL" == "no" ) then
    echo "#\!$FIND_BASH" > temp.txt
  endif
endif
# update script
cat $CHANGE_FILE >> temp.txt
cp -f temp.txt $CHANGE_FILE
rm -f temp.txt

exit 0
