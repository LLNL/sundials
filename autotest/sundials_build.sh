#!/bin/bash --login


############################################################################
# $Revision: 1.1 $
# $Date: 2004-04-02 21:18:20 $
############################################################################
#
# Filename: sundials_build.sh
# Programmer: Aaron Collier @ LLNL
#
############################################################################


# get build/configure options
get_build_opts()
{
  while read IN_MACHINE_NAME IN_REMOTE_DIR IN_CONFIGURE_OPTS; do
    # get information for local system
    if [ "${IN_MACHINE_NAME}" = "#${SYSTEM_NAME}" ]; then
      # determine remote directory to use
      if [ "${IN_REMOTE_DIR}" = "default_dir" ]; then
        REMOTE_DIR="${LOCAL_DIR}"
      else
        REMOTE_DIR="${IN_REMOTE_DIR}"
      fi
      echo "remote dir: ${REMOTE_DIR}"
      # determine configure options to use to build software
      if [ "${IN_CONFIGURE_OPTS}" = "default_config" ]; then
        CONFIGURE_OPTS=""
      else
        CONFIGURE_OPTS="${IN_CONFIGURE_OPTS}"
      fi
      echo -e "configure opts: ${CONFIGURE_OPTS}\n"
      return
    fi
  done
}

# get system name
get_system_name()
{
  CONTINUE="yes"
  NUM_CHARS=`echo "${SYSTEM_NAME}" | wc -m`
  NUM_CHARS=$((${NUM_CHARS} - 1))

  # eliminate node number digit by digit
  while [ "${CONTINUE}" != "no" ]; do
    TEST=`echo "${SYSTEM_NAME}" | cut -c$((${NUM_CHARS}))`
    TEST=${TEST%[0-9]}
    if [ "${TEST}" = "" ]; then
      SYSTEM_NAME=${SYSTEM_NAME%[0-9]}
    else
      CONTINUE="no"
    fi
    NUM_CHARS=$((${NUM_CHARS} - 1))
  done

  echo "${SYSTEM_NAME}" > temp.txt
}

############################################################################
#
# check the input file
#
############################################################################

# if input file not readable (does not exist) then abort
TEMP_INPUT_FILE="$1"
if [ ! -r "$TEMP_INPUT_FILE" ]; then
  echo -e "\a\n\tERROR: input file either unreadable or does not exist\n"
  exit 0
# if input file is readable then source input file
else
  echo -ne "\nSourcing input file..."
  # determine system name
  MACHINE_TYPE=`uname -s`
  if [ "${MACHINE_TYPE}" = "SunOS" ]; then
    SYSTEM_NAME=`hostname`
  else
    SYSTEM_NAME=`hostname -s`
  fi
  STATUS=`echo "${SYSTEM_NAME}" | fgrep "tux"`
  if [ "${STATUS}" = "" ]; then
    # if CASC Linux system then keep node number
    STATUS=`echo "${LOCAL_MACHINE}" | fgrep "tux"`
    # if LC cluster system then remove node number from name
    # clusters often have a login node pool (multiple nodes)
    # want cluster name not node name
    if [ "${STATUS}" = "" ]; then
      echo "${SYSTEM_NAME}" | get_system_name
      SYSTEM_NAME=`cat temp.txt`
      rm -f temp.txt
    fi
  fi
  STATUS=`source "${TEMP_INPUT_FILE}" 2>&1 | fgrep "error"`
  if [ "${STATUS}" = "" ]; then
    source "${TEMP_INPUT_FILE}"
    echo -e "[OK]\n\tINPUT_FILE = ${TEMP_INPUT_FILE}\n"
  # if cannot source input file then abort
  else
    echo "[FAILED]"
    exit 0
  fi
  # get build options (arguments to supply to configure script)
  get_build_opts < "${TEMP_INPUT_FILE}"
fi

############################################################################
#
# build software
#
############################################################################

# if remote directory does not exist then abort
if [ ! -d "${REMOTE_DIR}" ]; then
  echo -e "\a\n\tERROR: directory does not exist\n"
  exit 0
else
  # set marker to allow motd header to be removed from log file (need a login shell)
  echo "START_HERE"
  cd "${PROJECT_NAME}"
  # create shell script to perform software build
  touch "${SYSTEM_NAME}-build.sh"
  echo "#\!/bin/sh" >> "${SYSTEM_NAME}-build.sh"
  echo "" >> "${SYSTEM_NAME}-build.sh"
  echo "echo \"#################################\"" >> "${SYSTEM_NAME}-build.sh"
  echo "echo \"#  Running configure script...  #\"" >> "${SYSTEM_NAME}-build.sh"
  echo "echo \"#################################\"" >> "${SYSTEM_NAME}-build.sh"
  echo "./configure ${CONFIGURE_OPTS}" >> "${SYSTEM_NAME}-build.sh"
  echo "echo \"#######################\"" >> "${SYSTEM_NAME}-build.sh"
  echo "echo \"#  Running 'make'...  #\"" >> "${SYSTEM_NAME}-build.sh"
  echo "echo \"#######################\"" >> "${SYSTEM_NAME}-build.sh"
  echo "make" >> "${SYSTEM_NAME}-build.sh"
  echo "echo \"###############################\"" >> "${SYSTEM_NAME}-build.sh"
  echo "echo \"#  Running 'make install'...  #\"" >> "${SYSTEM_NAME}-build.sh"
  echo "echo \"###############################\"" >> "${SYSTEM_NAME}-build.sh"
  echo "make install" >> "${SYSTEM_NAME}-build.sh"
  echo "echo \"################################\"" >> "${SYSTEM_NAME}-build.sh"
  echo "echo \"#  Running 'make examples'...  #\"" >> "${SYSTEM_NAME}-build.sh"
  echo "echo \"################################\"" >> "${SYSTEM_NAME}-build.sh"
  echo "make examples" >> "${SYSTEM_NAME}-build.sh"
  echo "" >> "${SYSTEM_NAME}-build.sh"
  echo "exit 0" >> "${SYSTEM_NAME}-build.sh"
  chmod +x "${SYSTEM_NAME}-build.sh"
  # execute software build script
  sh ./${SYSTEM_NAME}-build.sh
fi

exit 0
