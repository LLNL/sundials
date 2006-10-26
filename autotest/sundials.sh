#!/bin/bash


############################################################################
# $Revision: 1.10 $
# $Date: 2006-10-26 21:47:12 $
############################################################################
#
# Filename: sundials.sh
# Programmer(s): Aaron Collier @ LLNL
#
############################################################################


############################################################################
#
# environment variables
#
############################################################################

MAIL_DOMAIN="llnl.gov"
PROJECT_CVSROOT="/home/casc/repository"
REMOTE_COPY_CMD="/usr/bin/scp"
#REMOTE_COPY_ARGS="-oProtocol=1"
REMOTE_COPY_ARGS=""
REMOTE_LOGIN_CMD="/usr/bin/ssh"
#REMOTE_LOGIN_ARGS="-1"
REMOTE_LOGIN_ARGS=""
MAIL_CMD="/bin/mail"
ERROR_FLAG="no"
ERROR_DETECTED="no"

# keywords used when searching build log files for errors
ERROR_0="error"
ERROR_1="failed"
ERROR_2="warning"
ERROR_2="error"
ERROR_3="illegal"
ERROR_4="failure"
ERROR_5="cannot"
ERROR_6="not found"

echo -e "\nStarting..."

############################################################################
#
# check syntax (must include input file)
#
############################################################################

if [ $# -eq 0 ]; then
  echo -e "\a\n\tERROR: missing input file\n"
  exit 0
else
  INPUT_FILE="$1"
fi

############################################################################
#
# check the input file
#
############################################################################

if [ ! -r $INPUT_FILE ]; then
  echo -e "\a\n\tERROR: input file either unreadable or does not exist\n"
  exit 0
else
  echo -ne "\nSourcing input file..."
  # determine absolute path to input file if not given
  EXEC_NAME=`basename ${INPUT_FILE}`
  if [ "${EXEC_NAME}" = "${INPUT_FILE}" ]; then
    CURRENT_DIR=`pwd`
    INPUT_FILE="${CURRENT_DIR}/${INPUT_FILE}"
  fi
  # attempt to source input file to define variable values
  STATUS=`source "${INPUT_FILE}" 2>&1 | fgrep "error"`
  if [ "${STATUS}" = "" ]; then
    source "${INPUT_FILE}"
    echo -e "[OK]\n\tINPUT_FILE = ${INPUT_FILE}\n"
  else
    echo "[FAILED]"
    exit 0
  fi
fi

############################################################################
#
# check environment variabls
#
############################################################################

echo -e "Checking environment variables:\n"
echo -ne "\tChecking ${PROJECT_NAME}..."
# need a project name in order to proceed
if [ "${PROJECT_NAME}" = "" ]; then
  echo -e "[FAILED - none given]\n"
  exit 0
else
  echo -e "[OK]\n\t\tPROJECT_NAME = ${PROJECT_NAME}\n"
fi

# set proper variable values based upon info from input file (project name)
MAIN_SCRIPT="${PROJECT_NAME}.sh"
BUILD_SCRIPT="${PROJECT_NAME}_build.sh"
TEST_SCRIPT="${PROJECT_NAME}_test.sh"
TIMEOUT_SCRIPT="${PROJECT_NAME}_timeout.sh"
FIX_BASH="${PROJECT_NAME}_check_bash.sh"
KEYWORD_ERROR_FILE="${PROJECT_NAME}_errors.txt"
MAIL_FILE="${PROJECT_NAME}.info"

############################################################################

echo -ne "\tChecking LOCAL_USERNAME..."
# determine local user name
if [ "${LOCAL_USERNAME}" = "" ]; then
  # guess current user if none given and update input file
  LOCAL_USERNAME=`whoami`
  cat "${INPUT_FILE}" | sed -e "s/LOCAL_USERNAME=/LOCAL_USERNAME=${LOCAL_USERNAME}/" > temp.txt
  cp -f temp.txt "${INPUT_FILE}"
  rm -f temp.txt
  echo -e "[DEFALUT]\n\t\tLOCAL_USERNAME = ${LOCAL_USERNAME}\n"
else
  echo -e "[OK]\n\t\tLOCAL_USERNAME = ${LOCAL_USERNAME}\n"
fi

############################################################################

# update entry in input file
update_file_entry()
{
  touch temp.txt

  while read IN_LINE; do
    if [ "${IN_LINE}" != "" ]; then
      # update LOCAL_MACHINE entry
      if [ "${IN_LINE}" = "LOCAL_MACHINE=${OLD_NAME}" ]; then
        echo "LOCAL_MACHINE=${LOCAL_MACHINE}" >> temp.txt
      else
        echo "${IN_LINE}" >> temp.txt
      fi
    fi
  done

  # update input file
  cp -f temp.txt "${INPUT_FILE}"
  rm -f temp.txt
}

# get system name
get_system_name()
{
  CONTINUE="yes"
  NUM_CHARS=`echo "${LOCAL_MACHINE}" | wc -m`
  NUM_CHARS=$((${NUM_CHARS} - 1))

  # eliminate node number digit by digit
  while [ "${CONTINUE}" != "no" ]; do
    TEST=`echo "${LOCAL_MACHINE}" | cut -c$((${NUM_CHARS}))`
    TEST=${TEST%[0-9]}
    if [ "${TEST}" = "" ]; then
      LOCAL_MACHINE=${LOCAL_MACHINE%[0-9]}
    else
      CONTINUE="no"
    fi
    NUM_CHARS=$((${NUM_CHARS} - 1))
  done

  echo "${LOCAL_MACHINE}" > temp.txt
}

echo -ne "\tChecking LOCAL_MACHINE..."

# guess system name to verify entry in input file if given
# do NOT want FQDN so use appropriate flag based upon machine type
TEMP_LOCAL_MACHINE="${LOCAL_MACHINE}"
MACHINE_TYPE=`uname -s`
if [ "${MACHINE_TYPE}" = "SunOS" ]; then
  LOCAL_MACHINE=`hostname`
else
  LOCAL_MACHINE=`hostname -s`
fi
# if CASC Linux system then keep node number
STATUS=`echo "${LOCAL_MACHINE}" | fgrep "tux"`
# if LC cluster system then remove node number from name
# clusters often have a login node pool (multiple nodes)
# want cluster name not node name
if [ "${STATUS}" = "" ]; then
  echo "${LOCAL_MACHINE}" | get_system_name
  LOCAL_MACHINE=`cat temp.txt`
  rm -f temp.txt
fi
LOCAL_MACHINE_CHECK="${LOCAL_MACHINE}"
LOCAL_MACHINE="${TEMP_LOCAL_MACHINE}"

# get local system name
if [ "${LOCAL_MACHINE}" = "" ]; then
  LOCAL_MACHINE="${LOCAL_MACHINE_CHECK}"
  # update input file
  cat "${INPUT_FILE}" | sed -e "s/LOCAL_MACHINE=/LOCAL_MACHINE=${LOCAL_MACHINE}/" > temp.txt
  cp -f temp.txt "${INPUT_FILE}"
  rm -f temp.txt
  echo -e "[DEFAULT]\n\t\tLOCAL_MACHINE = ${LOCAL_MACHINE} (system running script)\n"
else
  OLD_NAME="${LOCAL_MACHINE}"
  # if CASC Linux system then keep node number
  STATUS=`echo "${LOCAL_MACHINE}" | fgrep "tux"`
  # if LC cluster system then remove node number from name
  # clusters often have a login node pool (multiple nodes)
  # want cluster name not node name
  if [ "${STATUS}" = "" ]; then
    echo "${LOCAL_MACHINE}" | get_system_name
    LOCAL_MACHINE=`cat temp.txt`
    rm -f temp.txt
  fi
  if [ "${LOCAL_MACHINE}" != "${LOCAL_MACHINE_CHECK}" ]; then
    LOCAL_MACHINE="${LOCAL_MACHINE_CHECK}"
    echo -e "[UPDATED]\n\t\tLOCAL_MACHINE = ${LOCAL_MACHINE}\n"
  elif [ "${LOCAL_MACHINE}" != "${OLD_NAME}" ]; then
    echo -e "[UPDATED]\n\t\tLOCAL_MACHINE = ${LOCAL_MACHINE}\n"
  else
    echo -e "[OK]\n\t\tLOCAL_MACHINE = ${LOCAL_MACHINE}\n"
  fi
  # update input file
  update_file_entry < "${INPUT_FILE}"
fi

############################################################################

echo -ne "\tChecking LOCAL_DIR..."
# check local work directory (does it exist)
# directory where a few autotest-related files are stored
if [ "${LOCAL_DIR}" != "" -a -d "${LOCAL_DIR}" ]; then
  echo -e "[OK]\n\t\tLOCAL_DIR = ${LOCAL_DIR}\n"
# would rather the user specify a work directory than just use current directory
elif [ "${LOCAL_DIR}" = "" ]; then
  echo -e "[FAILED - none given]\n"
  exit 0
# if the given directory does not exist then abort
else
  echo -e "[FAILED - given directory does not exist]\n"
  exit 0
fi

############################################################################

echo -ne "\tChecking REMOTE_USERNAME..."
# if no remote user name given then assume same as local user name
if [ "${REMOTE_USERNAME}" = "" ]; then
  REMOTE_USERNAME="${LOCAL_USERNAME}"
  # update input file
  cat "${INPUT_FILE}" | sed -e "s/REMOTE_USERNAME=/REMOTE_USERNAME=${REMOTE_USERNAME}/" > temp.txt
  cp -f temp.txt "${INPUT_FILE}"
  rm -f temp.txt
  echo -e "[DEFAULT]\n\t\tREMOTE_USERNAME = ${REMOTE_USERNAME}\n"
else
  echo -e "[OK]\n\t\tREMOTE_USERNAME = ${REMOTE_USERNAME}\n"
fi

############################################################################

echo -ne "\tChecking REMOTE_MACHINES..."
# if no remote machines given then assume only testing local system and update input file
if [ "${REMOTE_MACHINES}" = "" ]; then
  REMOTE_MACHINES="${LOCAL_MACHINE}"
  cat "${INPUT_FILE}" | sed -e "s/REMOTE_MACHINES=/REMOTE_MACHINES=${REMOTE_MACHINES}/" > temp.txt
  cp -f temp.txt "${INPUT_FILE}"
  rm -f temp.txt
  echo -e "[DEFAULT]\n\t\tREMOTE_MACHINES = ${REMOTE_MACHINES}\n"
else
  echo -e "[OK]\n\t\tREMOTE_MACHINES = ${REMOTE_MACHINES}\n"
fi

############################################################################

# parse system info from input file
check_config_opts()
{
  SERIAL_ONLY=""
  SYS_WITH_OPTS=""
  NUM_OK="0"

  echo -e "\tChecking CONFIGURE_OPTS for system(s)..."
  while read MACHINE_NAME REMOTE_DIR CONFIGURE_OPTS; do
    if [ "${MACHINE_NAME}" != "" -a "${CONFIGURE_OPTS}" != "" -a "${REMOTE_DIR}" != "" ]; then
      TEMP_MACHINE_NAME=`echo "${MACHINE_NAME}" | sed -e "s/#//"`
      STATUS=`echo "${REMOTE_MACHINES}" | fgrep "${TEMP_MACHINE_NAME}"`
      if [ ! "${STATUS}" = "" ]; then
        NUM_OK=$((${NUM_OK} + 1))
        # the default_config keyword results in a standard/default software build
        if [ "${CONFIGURE_OPTS}" = "default_config" ]; then
          echo -e "\t\t[DEFAULT]   SYSTEM = ${TEMP_MACHINE_NAME}\n"
        else
          echo -e "\t\t[OK]        SYSTEM = ${TEMP_MACHINE_NAME}\n"
        fi
	# determine if building software with MPI support
	STATUS_SERIAL=`echo "${CONFIGURE_OPTS}" | fgrep -e "--without-mpi"`
	if [ "${STATUS_SERIAL}" = "" ]; then
	  TEMP_SERIAL="no"
        else
          TEMP_SERIAL="yes"
        fi
        if [ "${SERIAL_ONLY}" = "" ]; then
          SERIAL_ONLY="${TEMP_SERIAL}"
        else
          SERIAL_ONLY="${SERIAL_ONLY} ${TEMP_SERIAL}"
        fi
        # keep track of which systems have config options
        if [ "${SYS_WITH_OPTS}" = "" ]; then
          SYS_WITH_OPTS="${TEMP_MACHINE_NAME}"
        else
          SYS_WITH_OPTS="${SYS_WITH_OPTS} ${TEMP_MACHINE_NAME}"
        fi
      else
        if [ $((${NUM_OK})) -eq 0 ]; then
          echo ""
          NUM_OK=$((${NUM_OK} + 1))
        fi
      fi
    # if not enough info given then abort will illegal format message
    else
      STATUS=`echo "${MACHINE_NAME}" | fgrep "#"`
      if [ ! "${STATUS}" = "" ]; then
        echo -e "\t\t[FAILED]"
        echo -e "\a\n\tERROR: illegal format\n"
        exit 0
      fi
    fi
  done

  echo -ne "\tChecking for necessary system configuration directives..."
  # checking if each system listed in REMOTE_MACHINES has a corresponding configuration directive
  NUM_SYS_GIVEN=`echo "${REMOTE_MACHINES}" | wc -w`
  MISSING_SYS_OPTS=""
  WITH_SYS_OPTS=""
  while [ $((${NUM_SYS_GIVEN})) -gt 0 ]; do
    TEMP_SYS_GIVEN=`echo "${REMOTE_MACHINES}" | cut -d' ' -f$((${NUM_SYS_GIVEN}))`
    STATUS=`echo "${SYS_WITH_OPTS}" | fgrep "${TEMP_SYS_GIVEN}"`
    # system has a configuration directive
    if [ "${STATUS}" != "" ]; then
      if [ "${WITH_SYS_OPTS}" = "" ]; then
        WITH_SYS_OPTS="${TEMP_SYS_GIVEN}"
      else
        WITH_SYS_OPTS="${WITH_SYS_OPTS} ${TEMP_SYS_GIVEN}"
      fi
    # system does not have a configuration directive
    else
      if [ "${MISSING_SYS_OPTS}" = "" ]; then
        MISSING_SYS_OPTS="${TEMP_SYS_GIVEN}"
      else
        MISSING_SYS_OPTS="${MISSING_SYS_OPTS} ${TEMP_SYS_GIVEN}"
      fi
    fi
    NUM_SYS_GIVEN=$((${NUM_SYS_GIVEN} - 1))
  done
  # if at least one system does not have a configuration directive, then update REMOTE_MACHINES
  if [ "${MISSING_SYS_OPTS}" != "" ]; then
    # if no system has a configuration directive then abort
    if [ "${WITH_SYS_OPTS}" = "" ]; then
      echo -e "[FAILED - missing system configuration directive(s)]\n"
      exit 0
    # if at least one system has a configuration directive then continue
    else
      echo -e "[FAILED - missing system configuration directive(s)]\n"
      echo -e "\t\tSkipping the following system(s): ${MISSING_SYS_OPTS}\n"
    fi
  # if each systems given have a corresponding configuration directive then continue
  else
    echo -e "[OK]\n"
  fi

  echo "${SERIAL_ONLY}" > t2.txt
  echo "${WITH_SYS_OPTS}" > t3.txt
}

check_config_opts < "${INPUT_FILE}"
SERIAL_ONLY=`cat t2.txt`
REMOTE_MACHINES=`cat t3.txt`
rm -f t2.txt t3.txt

############################################################################

# update directory entry in input file
# sed command chokes on slashes in directory path names
update_dir_entry()
{
  touch temp.txt

  while read IN_LINE; do
    if [ "${IN_LINE}" != "" ]; then
      # update LOG_DIR entry
      if [ "${IN_LINE}" = "LOG_DIR=" -a "${UPDATE_DIR_ENTRY_FLAG}" = "log" ]; then
        echo "LOG_DIR=${LOG_DIR}" >> temp.txt
      # update SCRIPT_DIR entry
      elif [ "${IN_LINE}" = "SCRIPT_DIR=" -a "${UPDATE_DIR_ENTRY_FLAG}" = "script" ]; then
        echo "SCRIPT_DIR=${SCRIPT_DIR}" >> temp.txt
      # update MPI_OPT_DIR entry
      elif [ "${IN_LINE}" = "MPI_OPT_DIR=" -a "${UPDATE_DIR_ENTRY_FLAG}" = "mpi" ]; then
        echo "MPI_OPT_DIR=${MPI_OPT_DIR}" >> temp.txt
      else
        echo "${IN_LINE}" >> temp.txt
      fi
    fi
  done

  # update input file
  cp -f temp.txt "${INPUT_FILE}"
  rm -f temp.txt
}

echo -ne "\tChecking MPI_OPT_DIR..."
# check MPI_OPT_DIR (directory containing optional MPI configuration files)
if [ "${MPI_OPT_DIR}" = "" ]; then
  # by default MPI_OPT_DIR is directory containing INPUT_FILE
  MPI_OPT_DIR=`dirname "${INPUT_FILE}"`
  # update input file
  UPDATE_DIR_ENTRY_FLAG="mpi"
  update_dir_entry < "${INPUT_FILE}"
  echo -e "[DEFAULT]\n\t\tMPI_OPT_DIR = ${MPI_OPT_DIR}\n"
else
  # check if given directory exists
  if [ -d "${MPI_OPT_DIR}" ]; then
    echo -e "[OK]\n\t\tMPI_OPT_DIR = ${MPI_OPT_DIR}\n"
  # if directory does not exist then abort
  else
    echo -e "[ERROR - directory does not exist]\n"
    exit 0
  fi
fi

############################################################################

# check MPI options
OPT_EXISTS="no"
OPT_DIR="${MPI_OPT_DIR}"
echo -e "\tChecking MPI options..."
NUM_MACHINES=`echo "${REMOTE_MACHINES}" | wc -w`
while [ $((${NUM_MACHINES})) -gt 0 ]; do
  TEMP_ONE=`echo "${REMOTE_MACHINES}" | cut -d' ' -f$((${NUM_MACHINES}))`
  TEMP_TWO=`echo "${SERIAL_ONLY}" | cut -d' ' -f$((${NUM_MACHINES}))`
  UPDATE_FLAG="no"
  # if MPI support enabled then check MPI input file
  if [ "${TEMP_TWO}" = "no" ]; then
    OPT_EXISTS="yes"
    # if MPI input file exists then source
    if [ -r "${OPT_DIR}"/"${TEMP_ONE}".opt ]; then
      OPT_INPUT_FILE="${OPT_DIR}/${TEMP_ONE}.opt"
      OPT_EXEC_NAME=`basename ${OPT_INPUT_FILE}`
      if [ "${OPT_EXEC_NAME}" = "${OPT_INPUT_FILE}" ]; then
        CURRENT_DIR=`pwd`
	OPT_INPUT_FILE="${CURRENT_DIR}/${OPT_INPUT_FILE}"
      fi
      STATUS_OPT=`source "${OPT_INPUT_FILE}" 2>&1 | fgrep "error"`
      if [ "${STATUS_OPT}" = "" ]; then
        source "${OPT_INPUT_FILE}"
	# set missing values to "default" keyword
	if [ "${MPI_VERSION}" = "" ]; then
          MPI_VERSION="default"
	  UPDATE_FLAG="yes"
	fi
	if [ "${MPI_DIR}" = "default" ]; then
	  MPI_DIR=""
	  UPDATE_FLAG="yes"
	fi
	if [ "${MPI_COMMAND}" = "" ]; then
	  MPI_COMMAND="default"
	  UPDATE_FLAG="yes"
	fi
	# update MPI input file if necessary
	if [ "${UPDATE_FLAGE}" = "yes" ]; then
	  echo "MPI_VERSION=\"${MPI_VERSION}\"" > t2.txt
	  echo "MPI_DIR=\"${MPI_DIR}\"" >> t2.txt
	  echo "MPI_COMMAND=\"${MPI_COMMAND}\"" >> t2.txt
	  cp -f t2.txt "${OPT_INPUT_FILE}"
	  rm -f t2.txt
	fi
      # if cannot source input file then generate default file
      else
        rm -f "${OPT_INPUT_FILE}"
	touch "${OPT_INPUT_FILE}"
        echo "MPI_VERSION=\"default\"" > "${OPT_INPUT_FILE}"
        echo "MPI_DIR=\"\"" >> "${OPT_INPUT_FILE}"
        echo "MPI_COMMAND=\"default\"" >> "${OPT_INPUT_FILE}"
      fi
    # if MPI input file does not exist then create with default values
    else
      touch "${OPT_DIR}/${TEMP_ONE}.opt"
      OPT_INPUT_FILE="${OPT_DIR}/${TEMP_ONE}.opt"
      echo "MPI_VERSION=\"default\"" > "${OPT_INPUT_FILE}"
      echo "MPI_DIR=\"\"" >> "${OPT_INPUT_FILE}"
      echo "MPI_COMMAND=\"default\"" >> "${OPT_INPUT_FILE}"
    fi
    echo -e "\t\t[PARALLEL]   SYSTEM = ${TEMP_ONE}\n"
  # not building software with MPI support
  else
    echo -e "\t\t[SERIAL]     SYSTEM = ${TEMP_ONE}\n"
  fi
  NUM_MACHINES=$((${NUM_MACHINES} - 1))
done

############################################################################

echo -ne "\tChecking MPI_HOST_FILE_DIR..."
# check MPI host file directory
if [ "${MPI_HOST_FILE_DIR}" != "" ]; then
  if [ -d "${MPI_HOST_FILE_DIR}" ]; then
    echo -e "[OK]\n\t\tMPI_HOST_FILE_DIR = ${MPI_HOST_FILE_DIR}\n"
  else
    echo -e "[ERROR - directory does not exist]\n"
    exit 0
  fi
fi

############################################################################

# check if a remote directory was given (does not confirm existence)
# the remote directory is where the software is built and installed
check_remote_dir()
{
  echo -e "\tChecking REMOTE_DIR for system(s)..."
  while read MACHINE_NAME REMOTE_DIR CONFIGURE_OPTS; do
    if [ "${MACHINE_NAME}" != "" -a "${CONFIGURE_OPTS}" != "" -a "${REMOTE_DIR}" != "" ]; then
      TEMP_MACHINE_NAME=`echo "${MACHINE_NAME}" | sed -e "s/#//"`
      STATUS=`echo "${REMOTE_MACHINES}" | fgrep "${TEMP_MACHINE_NAME}"`
      if [ ! "${STATUS}" = "" ]; then
        # use of the default_dir keyword implies REMOTE_DIR = LOCAL_DIR
        if [ "${REMOTE_DIR}" = "default_dir" ]; then
          echo -e "\t\t[DEFAULT]   REMOTE_DIR = ${LOCAL_DIR} (local directory)"
          echo -e "\t\t            SYSTEM = ${TEMP_MACHINE_NAME}\n"
        else
          echo -e "\t\t[OK]        REMOTE_DIR = ${REMOTE_DIR}"
          echo -e "\t\t            SYSTEM = ${TEMP_MACHINE_NAME}\n"
        fi
      fi
    fi
  done
}

check_remote_dir < "${INPUT_FILE}"

############################################################################

echo -ne "\tChecking MAIL_RECIPIENTS..."
# check e-mail list (addresses where error messages are sent)
if [ "${MAIL_RECIPIENTS}" = "" ]; then
  # if none given then assume local user should be sent e-mail containing error messages
  MAIL_RECIPIENTS="${LOCAL_USERNAME}@${MAIL_DOMAIN}"
  # update input file
  cat "${INPUT_FILE}" | sed -e "s/MAIL_RECIPIENTS=/MAIL_RECIPIENTS=${MAIL_RECIPIENTS}/" > temp.txt
  cp -f temp.txt "${INPUT_FILE}"
  rm -f temp.txt
  echo -e "[DEFAULT]\n\t\tMAIL_RECIPIENTS = ${MAIL_RECIPIENTS}\n"
# make certain each entry in list is an actual e-mail address
else
  TEMP_MAIL_RECIPIENTS=""
  NUM_PEOPLE=`echo "${MAIL_RECIPIENTS}" | wc -w`
  while [ $((${NUM_PEOPLE})) -gt 0 ]; do
    TEMP=`echo "${MAIL_RECIPIENTS}" | cut -d' ' -f$((${NUM_PEOPLE}))`
    # if only user name given then assume e-mail address is actually user_name@MAIL_DOMAIN
    STATUS=`echo "${TEMP}" | fgrep "@"`
    if [ "${STATUS}" = "" ]; then
      TEMP="${TEMP}@${MAIL_DOMAIN}"
    fi
    if [ "${TEMP_MAIL_RECIPIENTS}" = "" ]; then
      TEMP_MAIL_RECIPIENTS="${TEMP}"
    else
      TEMP_MAIL_RECIPIENTS="${TEMP_MAIL_RECIPIENTS} ${TEMP}"
    fi
    NUM_PEOPLE=$((${NUM_PEOPLE} - 1))
  done
  # update input file
  cat "${INPUT_FILE}" | sed -e "s/MAIL_RECIPIENTS=\"${MAIL_RECIPIENTS}\"/MAIL_RECIPIENTS=\"${TEMP_MAIL_RECIPIENTS}\"/" > temp.txt
  cp -f temp.txt "${INPUT_FILE}"
  rm -f temp.txt
  MAIL_RECIPIENTS="${TEMP_MAIL_RECIPIENTS}"
  echo -e "[OK]\n\t\tMAIL_RECIPIENTS = ${MAIL_RECIPIENTS}\n"
fi

############################################################################

echo -ne "\tChecking AUTOTEST_ADMIN..."
if [ "${AUTOTEST_ADMIN}" = "" ]; then
  # if none given then assume local user should be sent main script output in e-mail form
  AUTOTEST_ADMIN="${LOCAL_USERNAME}@${MAIL_DOMAIN}"
  # update input file
  cat "${INPUT_FILE}" | sed -e "s/AUTOTEST_ADMIN=/AUTOTEST_ADMIN=${AUTOTEST_ADMIN}/" > temp.txt
  cp -f temp.txt "${INPUT_FILE}"
  rm -f temp.txt
  echo -e "[DEFAULT]\n\t\tAUTOTEST_ADMIN = ${AUTOTEST_ADMIN}\n"
# appdend default mail domain to user name (if none given) to form e-mail address
else
  STATUS=`echo "${AUTOTEST_ADMIN}" | fgrep "@"`
  if [ "${STATUS}" = "" ]; then
    TEMP_AUTOTEST_ADMIN="${AUTOTEST_ADMIN}"
    AUTOTEST_ADMIN="${AUTOTEST_ADMIN}@${MAIL_DOMAIN}"
    # update input file
    cat "${INPUT_FILE}" | sed -e "s/AUTOTEST_ADMIN=${TEMP_AUTOTEST_ADMIN}/AUTOTEST_ADMIN=${AUTOTEST_ADMIN}/" > temp.txt
    cp -f temp.txt "${INPUT_FILE}"
    rm -f temp.txt
    echo -e "[OK]\n\t\tAUTOTEST_ADMIN = ${AUTOTEST_ADMIN}\n"
  else
    echo -e "[OK]\n\t\tAUTOTEST_ADMIN = ${AUTOTEST_ADMIN}\n"
  fi
fi

############################################################################

# check maximum running time
# running time is measured in units of minutes
echo -ne "\tChecking MAX_TIME..."
# must give maximum running time to prevent infinite loops from causing problems
if [ "${MAX_TIME}" = "" ]; then
  echo "[ERROR - MAX_TIME not specified]"
  exit 0
# maximum running time must be positive
else
  if [ $((${MAX_TIME})) -lt 1 ]; then
    echo "[ERROR - MAX_TIME must be >= 1]"
    exit 0
  else
    echo -e "[OK]\n\t\tMAX_TIME = ${MAX_TIME} minutes\n"
  fi
fi

############################################################################

echo -ne "\tChecking LOG_DIR..."
# check log directory (does it exist)
# directory where log files/archives are stored on local system
if [ "${LOG_DIR}" != "" ]; then
  if [ -d "${LOG_DIR}" ]; then
    echo -e "[OK]\n\t\tLOG_DIR = ${LOG_DIR}\n"
  else
    mkdir "${LOG_DIR}"
    echo -e "[OK]\n\t\tLOG_DIR = ${LOG_DIR}\n"
  fi
elif [ "${LOG_DIR}" = "" ]; then
  # if not given assume LOG_DIR = LOCAL_DIR
  LOG_DIR="${LOCAL_DIR}/log"
  UPDATE_DIR_ENTRY_FLAG="log"
  # update input file
  update_dir_entry < "${INPUT_FILE}"
  echo -e "[DEFAULT]\n\t\tLOG_DIR = ${LOG_DIR}\n"
fi

############################################################################

echo -ne "\tChecking NUM_LOG_FILES..."
# check how may log files to keep
# if number < ZERO then abort
if [ $((${NUM_LOG_FILES})) -lt 0 ]; then
  echo -e "[ERROR - cannot be negative]\n"
  exit 0
# if no value given or value is zero then set to default value
# really only need the later case
elif [ "${NUM_LOG_FILES}" = "" -o $((${NUM_LOG_FILES})) -eq 0 ]; then
  if [ "${NUM_LOG_FILES}" = "" ]; then
    NUM_LOG_FILES="1"
    cat "${INPUT_FILE}" | sed -e "s/NUM_LOG_FILES=/NUM_LOG_FILES=${NUM_LOG_FILES}/" > temp.txt
  else
    NUM_LOG_FILES="1"
    cat "${INPUT_FILE}" | sed -e "s/NUM_LOG_FILES=0/NUM_LOG_FILES=${NUM_LOG_FILES}/" > temp.txt
  fi
  cp -f temp.txt "${INPUT_FILE}"
  rm -f temp.txt
  echo -e "[DEFAULT]\n\t\tNUM_LOG_FILES = ${NUM_LOG_FILES}\n"
# valid number was given
elif [ $((${NUM_LOG_FILES})) -gt 10 ]; then
  echo -e "[ERROR - cannot be larger than 10]\n"
  exit 0
else
  echo -e "[OK]\n\t\tNUM_LOG_FILES = ${NUM_LOG_FILES}\n"
fi

############################################################################

echo -ne "\tChecking SCRIPT_DIR..."
# check script directory
if [ "${SCRIPT_DIR}" = "" ]; then
  # check if script directory if current working directory
  SCRIPT_DIR=`pwd`
  # if main autotest script is in current working directory then guess is correct
  if [ -x "${SCRIPT_DIR}/${MAIN_SCRIPT}" ]; then
    UPDATE_DIR_ENTRY_FLAG="script"
    # update input file
    update_dir_entry < "${INPUT_FILE}"
    echo -e "[DEFAULT]\n\t\tSCRIPT_DIR = ${SCRIPT_DIR}\n"
  # if could not guess script directory then abort
  else
    echo -e "[ERROR - unable to guess directory]\n"
    exit 0
  fi
elif [ "${SCRIPT_DIR}" != "" ]; then
  # check if given script directory actually exists
  if [ -d "${SCRIPT_DIR}" ]; then
    # check if given script directory actually contains the autotest scripts
    if [ -x "${SCRIPT_DIR}/${MAIN_SCRIPT}" ]; then
      echo -e "[OK]\n\t\tSCRIPT_DIR = ${SCRIPT_DIR}\n"
    # if not then abort
    else
      echo -e "[ERROR - cannot find scripts]\n"
      exit 0
    fi
  # if given script directory does not exist then abort
  else
    echo -e "[ERROR - directory does not exist]\n"
    exit 0
  fi
fi

############################################################################

# check if should use CVS to retrieve project files
if [ "${SOURCE_DIR}" = "" ]; then
  USE_CVS="yes"
else
  echo -ne "\tChecking SOURCE_DIR..."
  # cannot use scratch directory
  if [ "${SOURCE_DIR}" = "${LOCAL_DIR}/${LOCAL_MACHINE}/${PROJECT_NAME}" ]; then
    echo -e "[ERROR - cannot use scratch directory]\n"
    exit 0
  # if given directory exists then tar contents of directory for distribution
  elif [ -d "${SOURCE_DIR}" ]; then
    USE_CVS="no"
    echo -e "[OK]\n\t\tSOURCE_DIR = ${SOURCE_DIR}\n"
  # if directory does not exist then abort
  else
    echo -e "[ERROR - directory does not exist]\n"
    exit 0
  fi
fi

############################################################################
#
# copy files to local and remote machine(s)
#
############################################################################

# remove old work directory if necessary
if [ -d "${LOCAL_DIR}/${LOCAL_MACHINE}" ]; then
  rm -Rf "${LOCAL_DIR}/${LOCAL_MACHINE}"
fi
mkdir "${LOCAL_DIR}/${LOCAL_MACHINE}"

if [ "${USE_CVS}" = "yes" ]; then
  echo -n "Checking ${PROJECT_NAME} out of CVS repository..."
  # get project files out of CVS repository
  # the destination directory is the local work directory
  cd "${LOCAL_DIR}/${LOCAL_MACHINE}"
  STATUS=`cvs -d "${PROJECT_CVSROOT}" export -D now "${PROJECT_NAME}" 2>&1 | fgrep "cannot"`
  # if CVS checkout failed then abort
  if [ ! "${STATUS}" = "" ]; then
    echo "[FAILED]"
    echo -e "\a\n\tERROR: CVS checkout failed\n"
    exit 0
  else
    cvs -d "${PROJECT_CVSROOT}" export -D now "${PROJECT_NAME}" 1> /dev/null 2> /dev/null
    echo -e "[DONE]\n"
  fi
# use given source tree
else
  echo -n "Using ${SOURCE_DIR} as source directory..."
  cp -pR "${SOURCE_DIR}" "${LOCAL_DIR}/${LOCAL_MACHINE}"
  cd "${LOCAL_DIR}/${LOCAL_MACHINE}"
  echo -e "[DONE]\n"
fi

############################################################################

echo -n "Moving files to appropriate directory on local system..."
# copy all necessary scripts/files to the local work directory on local system
cp -f "${SCRIPT_DIR}/${BUILD_SCRIPT}" "${LOCAL_DIR}/${LOCAL_MACHINE}"
cp -f "${SCRIPT_DIR}/${TEST_SCRIPT}" "${LOCAL_DIR}/${LOCAL_MACHINE}"
cp -f "${SCRIPT_DIR}/${FIX_BASH}" "${LOCAL_DIR}/${LOCAL_MACHINE}"
cp -f "${SCRIPT_DIR}/${TIMEOUT_SCRIPT}" "${LOCAL_DIR}/${LOCAL_MACHINE}"
cp -f "${INPUT_FILE}" "${LOCAL_DIR}/${LOCAL_MACHINE}"
# copy MPI option file(s) only if used
if [ "${OPT_EXISTS}" = "yes" ]; then
  cp -f "${OPT_DIR}"/*.opt "${LOCAL_DIR}/${LOCAL_MACHINE}"
fi
if [ "${MPI_HOST_FILE_DIR}" != "" ]; then
  cp -f "${MPI_HOST_FILE_DIR}"/*.host "${LOCAL_DIR}/${LOCAL_MACHINE}"
fi
echo -e "[DONE]\n"

############################################################################

# determine remote directory
find_remote_dir()
{
  while read IN_MACHINE_NAME IN_REMOTE_DIR IN_CONFIGURE_OPTS; do
    if [ "${IN_MACHINE_NAME}" != "" -a "${IN_CONFIGURE_OPTS}" != "" -a "${IN_REMOTE_DIR}" != "" ]; then
      TEMP_MACHINE_NAME=`echo "${IN_MACHINE_NAME}" | sed -e "s/#//"`
      FIND_STATUS=`echo "${REMOTE_MACHINES}" | fgrep "${TEMP_MACHINE_NAME}"`
      if [ "${FIND_STATUS}" != "" -a "${TEMP_MACHINE_NAME}" = "${TEMP}" ]; then
        if [ "${IN_REMOTE_DIR}" = "default_dir" ]; then
          TEMP_REMOTE_DIR="${LOCAL_DIR}"
        else
          TEMP_REMOTE_DIR="${IN_REMOTE_DIR}"
        fi
        echo "${TEMP_REMOTE_DIR}" > temp-dir.txt
        return
      fi
    fi
  done 
}

echo "Copying files to remote system(s):"
# copy project files to remote systems
LIST_REMOTE_DIRS=""
REMOVE_MACHINES=""
STATUS=`echo "${REMOTE_MACHINES}" | fgrep "${LOCAL_MACHINE}"`
NUM_MACHINES=`echo "${REMOTE_MACHINES}" | wc -w`
# delete old project tar file
if [ "${STATUS}" = "" -o $((${NUM_MACHINES})) -gt 1 ]; then
  if [ -e "${PROJECT_NAME}-test.tar.gz" ]; then
    rm -f "${PROJECT_NAME}-test.tar.gz"
  fi
  # create project tar ball (compressed)
  TAR_STATUS=`tar -cf "${PROJECT_NAME}"-test.tar "${PROJECT_NAME}" "${BUILD_SCRIPT}" "${TEST_SCRIPT}" "${FIX_BASH}" "${TIMEOUT_SCRIPT}" "${EXEC_NAME}" *.opt *.host 2>&1 | fgrep -i "error"`
  # if tar command fails then abort
  if [ ! "${TAR_STATUS}" = "" ]; then
    echo -e "\a\n\tERROR: tar utility failed\n"
    exit 0
  fi
  tar -cf "${PROJECT_NAME}"-test.tar "${PROJECT_NAME}" "${BUILD_SCRIPT}" "${TEST_SCRIPT}" "${FIX_BASH}" "${TIMEOUT_SCRIPT}" "${FORTRAN_SCRIPT}" "${EXEC_NAME}" *.opt *.host 1> /dev/null 2> /dev/null
  gzip "${PROJECT_NAME}"-test.tar 1> /dev/null 2> /dev/null
  while [ $((${NUM_MACHINES})) -gt 0 ]; do
    TEMP=`echo "${REMOTE_MACHINES}" | cut -d' ' -f$((${NUM_MACHINES}))`
    if [ ! "${TEMP}" = "${LOCAL_MACHINE}" ]; then
      # determine remote directory associated with current system name (check input file)
      find_remote_dir < "${INPUT_FILE}"
      TEMP_REMOTE_DIR=`cat temp-dir.txt`
      rm -f temp-dir.txt
      # check for duplicate destination/remote directories to avoid resending project tar ball
      DIR_STATUS=`echo "${LIST_REMOTE_DIRS}" | fgrep "${TEMP_REMOTE_DIR}"`
      # if remote system then use the remote copy command
      if [ "${DIR_STATUS}" = "" -a "${TEMP_REMOTE_DIR}" != "${LOCAL_DIR}" ]; then
        LIST_REMOTE_DIRS="${LIST_REMOTE_DIRS} ${TEMP_REMOTE_DIR}"
        echo -ne "\tCopying file(s) to ${TEMP}..."
        CPY_STATUS=`"${REMOTE_COPY_CMD}" "${REMOTE_COPY_ARGS}" "${PROJECT_NAME}-test.tar.gz" "${REMOTE_USERNAME}"@"${TEMP}":. 2>&1 | fgrep "not"`
        # if copy failed then flag system name for removal from list
        if [ ! "${CPY_STATUS}" = "" ]; then
          echo -e "[FAILED - unknown host]\n"
          if [ ! "${REMOVE_MACHINES}" = "" ]; then
            REMOVE_MACHINES="${REMOVE_MACHINES} ${TEMP}"
          else
            REMOVE_MACHINES="${TEMP}"
          fi
        else
          "${REMOTE_COPY_CMD}" "${REMOTE_COPY_ARGS}" "${PROJECT_NAME}-test.tar.gz" "${REMOTE_USERNAME}"@"${TEMP}":. 1> /dev/null 2> /dev/null
          echo -e "[DONE]\n"
        fi
      # if local system then use normal copy command rather than remote copy command
      else
        if [ "${TEMP_REMOTE_DIR}" = "${LOCAL_DIR}" ]; then
          cp -f "${LOCAL_DIR}/${LOCAL_MACHINE}/${PROJECT_NAME}-test.tar.gz" ~/
        fi
        echo -e "\tCopying file(s) to ${TEMP}...[SKIPPING - mounts shared directory]\n"
      fi
    fi
    NUM_MACHINES=$((${NUM_MACHINES}-1))
  done
# if only testing local system then just continue
else
  echo -e "\t[SKIPPING] (only local system given)\n"
fi

############################################################################

# remove systems from list for which the remote copy failed
if [ ! "${REMOVE_MACHINES}" = "" ]; then
  NUM_MACHINES=`echo "${REMOTE_MACHINES}" | wc -w`
  NUM_BAD_MACHINES=`echo "${REMOVE_MACHINES}" | wc -w`
  # if all systems were unreachable then abort
  if [ $((${NUM_BAD_MACHINES})) -eq $((${NUM_MACHINES})) ]; then
    echo -e "\a\n\tERROR: all systems are unreachable\n"
    exit 0
  fi
  # edit the system list if a system failed the previous step
  if [ $((${NUM_BAD_MACHINES})) -gt 0 ]; then
    TEMP=""
    while [ $((${NUM_MACHINES})) -gt 0 ]; do
      TEMP_REMOTE=`echo "${REMOTE_MACHINES}" | cut -d' ' -f$((${NUM_MACHINES}))`
      NUM_BAD_MACHINES=`echo "${REMOVE_MACHINES}" | wc -w`
      STATUS="good"
      while [ $((${NUM_BAD_MACHINES})) -gt 0 ]; do
        TEMP_BAD=`echo "${REMOVE_MACHINES}" | cut -d' ' -f$((${NUM_BAD_MACHINES}))`
        if [ "${TEMP_REMOTE}" = "${TEMP_BAD}" ]; then
          STATUS="bad"
        fi
        NUM_BAD_MACHINES=$((${NUM_BAD_MACHINES}-1))
      done
      if [ "${STATUS}" = "good" ]; then
        if [ ! "${TEMP}" = "" ]; then
          TEMP="${TEMP} ${TEMP_REMOTE}"
        else
          TEMP="${TEMP_REMOTE}"
        fi
      fi
      NUM_MACHINES=$((${NUM_MACHINES}-1))
    done
    echo -e "Skipping the following system(s): ${REMOVE_MACHINES}\n"
    REMOTE_MACHINES="${TEMP}"
  fi
fi

############################################################################
#
# specify log directory
#
############################################################################

# initialize scratch directory (LOG_DIR/examples)
if [ -d "${LOG_DIR}" ]; then
  rm -Rf "${LOG_DIR}/examples"
  mkdir "${LOG_DIR}/examples"
# make LOG_DIR and subdirectories if do not exist (LOG_DIR/log is created later)
elif [ ! -d "${LOG_DIR}" ]; then
  mkdir "${LOG_DIR}"
  mkdir "${LOG_DIR}/examples"
fi

############################################################################

# initialize mail file (contains error messages)
if [ -f "${LOG_DIR}/${MAIL_FILE}" ]; then
  rm -f "${LOG_DIR}/${MAIL_FILE}"
fi
touch "${LOG_DIR}/${MAIL_FILE}"
TIME_STAMP=`date`
echo "##################################" >> "${LOG_DIR}/${MAIL_FILE}"
echo "#" >> "${LOG_DIR}/${MAIL_FILE}"
echo "#  AUTOTEST INFORMATION" >> "${LOG_DIR}/${MAIL_FILE}"
echo "#" >> "${LOG_DIR}/${MAIL_FILE}"
echo "#  Tested System(s): ${REMOTE_MACHINES}" >> "${LOG_DIR}/${MAIL_FILE}"
echo "#" >> "${LOG_DIR}/${MAIL_FILE}"
echo "#  Date: ${TIME_STAMP}" >> "${LOG_DIR}/${MAIL_FILE}"
echo "#" >> "${LOG_DIR}/${MAIL_FILE}"
echo "##################################" >> "${LOG_DIR}/${MAIL_FILE}"
echo "" >> "${LOG_DIR}/${MAIL_FILE}"
echo "" >> "${LOG_DIR}/${MAIL_FILE}"

############################################################################
#
# build software
#
############################################################################

echo "Building ${PROJECT_NAME} on system(s):"
# buile software on each system
NUM_MACHINES=`echo "${REMOTE_MACHINES}" | wc -w`
while [ $((${NUM_MACHINES})) -gt 0 ]; do
  TEMP_REMOTE=`echo "${REMOTE_MACHINES}" | cut -d' ' -f$((${NUM_MACHINES}))`
  echo -ne "\tBuilding ${PROJECT_NAME} on ${TEMP_REMOTE}..."
  # if remote system then execute build commands via remote login command
  if [ ! "${TEMP_REMOTE}" = "${LOCAL_MACHINE}" ]; then
    TEMP="${TEMP_REMOTE}"
    # determine remote directory (check input file)
    find_remote_dir < "${INPUT_FILE}"
    REMOTE_DIR=`cat temp-dir.txt`
    rm -f temp-dir.txt
    # determine absolute path to tar command on remote system
    TAR_EXEC=`${REMOTE_LOGIN_CMD} ${REMOTE_LOGIN_ARGS} -l ${REMOTE_USERNAME} ${TEMP_REMOTE} " which tar "`
    # remove old installation directory
    # uncompress and extract contents of project tar ball
    ${REMOTE_LOGIN_CMD} ${REMOTE_LOGIN_ARGS} -l ${REMOTE_USERNAME} ${TEMP_REMOTE} \
    " rm -Rf ${REMOTE_DIR}/${TEMP_REMOTE} && \
    mkdir ${REMOTE_DIR}/${TEMP_REMOTE} && \
    cp -f ${PROJECT_NAME}-test.tar.gz ${REMOTE_DIR}/${TEMP_REMOTE} && \
    cd ${REMOTE_DIR}/${TEMP_REMOTE} && \
    gunzip -f ${PROJECT_NAME}-test.tar.gz && \
    ${TAR_EXEC} -xf ${PROJECT_NAME}-test.tar && \
    ./${FIX_BASH} ${BUILD_SCRIPT} yes && \
    ./${BUILD_SCRIPT} ${EXEC_NAME} " &> "${LOG_DIR}/${TEMP_REMOTE}"-build.log
    echo -e "[DONE]\n"
  # if local system then just change to appropriate directory and build software
  else
    cd "${LOCAL_DIR}/${LOCAL_MACHINE}" && \
    ./"${FIX_BASH}" "${BUILD_SCRIPT}" no && \
    ./"${BUILD_SCRIPT}" "${EXEC_NAME}" &> "${LOG_DIR}/${LOCAL_MACHINE}"-build.log
    echo -e "[DONE]\n"
  fi
  NUM_MACHINES=$((${NUM_MACHINES}-1))
done

############################################################################
#
# search for build/compile error messages
#
############################################################################

# remove motd header from log file before performing keyword search
remove_motd_header()
{
  while read IN_LINES IN_FILE; do
    if [ "${IN_LINES}" != "" -a "${IN_FILE}" != "" ]; then
      echo "TOTAL_LINES=${IN_LINES}" > temp.txt
      return
    fi
  done
}

BUILD_ERROR_LIST=""
echo "Checking build log file(s) for possible error messages:"
# check software build log files for error messages
NUM_MACHINES=`echo "${REMOTE_MACHINES}" | wc -w`
while [ $((${NUM_MACHINES})) -gt 0 ]; do
  TEMP_REMOTE=`echo "${REMOTE_MACHINES}" | cut -d' ' -f$((${NUM_MACHINES}))`
  echo -ne "\tChecking build log file(s) from ${TEMP_REMOTE}..."
  # remove motd header to avoid erroneous error messages
  # keyword search for error terms could return a flase-positive otherwise
  SKIP_LINES=`fgrep -n "START_HERE" "${LOG_DIR}/${TEMP_REMOTE}"-build.log | cut -d':' -f1`
  TOTAL_LINES=`wc -l "${LOG_DIR}/${TEMP_REMOTE}"-build.log`
  echo "${TOTAL_LINES}" | remove_motd_header
  source temp.txt
  SKIP_LINES=$((${TOTAL_LINES} - ${SKIP_LINES}))
  rm -f temp.txt
  MACHINE_TYPE=`uname -s`
  if [ "${MACHINE_TYPE}" != "SunOS" ]; then
    tail -n $((${SKIP_LINES})) "${LOG_DIR}/${TEMP_REMOTE}"-build.log > temp.txt
  else
    tail -$((${SKIP_LINES})) "${LOG_DIR}/${TEMP_REMOTE}"-build.log > temp.txt
  fi
  cp -f temp.txt "${LOG_DIR}/${TEMP_REMOTE}"-build.log
  rm -f temp.txt
  # search for keywords that may indicate an error occurred during software build procedure
  # case-insensitive keyword search (-i flag)
  STATUS_0=`fgrep -i "${ERROR_0}" "${LOG_DIR}/${TEMP_REMOTE}"-build.log`
  STATUS_1=`fgrep -i "${ERROR_1}" "${LOG_DIR}/${TEMP_REMOTE}"-build.log`
  # Thunder displays a trivial error that can safely be ignored
  if [ "${LOCAL_MACHINE}" = "thunder" ]; then
    :
  else
    STATUS_2=`fgrep -i "${ERROR_2}" "${LOG_DIR}/${TEMP_REMOTE}"-build.log`
  fi
  STATUS_3=`fgrep -i "${ERROR_3}" "${LOG_DIR}/${TEMP_REMOTE}"-build.log`
  STATUS_4=`fgrep -i "${ERROR_4}" "${LOG_DIR}/${TEMP_REMOTE}"-build.log`
  STATUS_5=`fgrep -i "${ERROR_5}" "${LOG_DIR}/${TEMP_REMOTE}"-build.log`
  STATUS_6=`fgrep -i "${ERROR_6}" "${LOG_DIR}/${TEMP_REMOTE}"-build.log`
  # if an error was detected then set error flag
  if [ "${STATUS_0}" != "" -o "${STATUS_1}" != "" -o "${STATUS_2}" != "" -o "${STATUS_3}" != "" -o "${STATUS_4}" != "" -o "${STATUS_5}" != "" -o "${STATUS_6}" != "" ]; then
    ERROR_FLAG="yes"
    echo -e "[ERROR]\n"
    # flag system if build failed
    if [ "${BUILD_ERROR_LIST}" = "" ]; then
      BUILD_ERROR_LIST="${TEMP_REMOTE}"
    else
      BUILD_ERROR_LIST="${BUILD_ERROR_LIST} ${TEMP_REMOTE}"
    fi
    # add ERROR tag to log files containing errors
    if [ -r "${LOG_DIR}/${TEMP_REMOTE}"-build.log ]; then
      mv -f "${LOG_DIR}/${TEMP_REMOTE}"-build.log "${LOG_DIR}/${TEMP_REMOTE}"-build.log-ERROR
      # add informative header and log file to e-mail message
      echo "##################################" >> "${LOG_DIR}/${MAIL_FILE}"
      echo "#" >> "${LOG_DIR}/${MAIL_FILE}"
      echo "#  BUILD ERROR" >> "${LOG_DIR}/${MAIL_FILE}"
      echo "#" >> "${LOG_DIR}/${MAIL_FILE}"
      echo "#  SYSTEM: ${TEMP_REMOTE}" >> "${LOG_DIR}/${MAIL_FILE}"
      echo "#  LOG FILE: ${TEMP_REMOTE}-build.log" >> "${LOG_DIR}/${MAIL_FILE}"
      echo "#" >> "${LOG_DIR}/${MAIL_FILE}"
      echo "##################################" >> "${LOG_DIR}/${MAIL_FILE}"
      cat "${LOG_DIR}/${TEMP_REMOTE}"-build.log-ERROR >> "${LOG_DIR}/${MAIL_FILE}"
    fi
  else
    echo -e "[DONE]\n"
  fi
  NUM_MACHINES=$((${NUM_MACHINES}-1))
done

############################################################################

# remove systems from list which failed to build software
if [ "${ERROR_FLAG}" = "yes" ]; then
  echo -e "Skipping the following system(s) because of possible build error(s): ${BUILD_ERROR_LIST}\n"
  NUM_MACHINES=`echo "${REMOTE_MACHINES}" | wc -w`
  NUM_BAD_MACHINES=`echo "${BUILD_ERROR_LIST}" | wc -w`
  # if all systems failed to build software then abort
  if [ $((${NUM_BAD_MACHINES})) -eq $((${NUM_MACHINES})) ]; then
    echo -e "\a\n\tERROR: all systems failed to build ${PROJECT_NAME}\n"
    # e-mail error message to appropriate people
    echo -n "E-mailing exit status to ${MAIL_RECIPIENTS} before exiting..."
    MAIL_SUBJECT="Autotest Message (${PROJECT_NAME}): Error"
    "${MAIL_CMD}" -s "${MAIL_SUBJECT}" "${MAIL_RECIPIENTS}" < "${LOG_DIR}/${MAIL_FILE}"
    echo -e "[DONE]\n"
    exit 0
  fi
  # if not all systems failed to successfully complete software build procedure then edit system list
  if [ $((${NUM_BAD_MACHINES})) -gt 0 ]; then
    TEMP=""
    while [ $((${NUM_MACHINES})) -gt 0 ]; do
      TEMP_REMOTE=`echo "${REMOTE_MACHINES}" | cut -d' ' -f$((${NUM_MACHINES}))`
      NUM_BAD_MACHINES=`echo "${BUILD_ERROR_LIST}" | wc -w`
      STATUS="good"
      while [ $((${NUM_BAD_MACHINES})) -gt 0 ]; do
        TEMP_BAD=`echo "${BUILD_ERROR_LIST}" | cut -d' ' -f$((${NUM_BAD_MACHINES}))`
        if [ "${TEMP_REMOTE}" = "${TEMP_BAD}" ]; then
          STATUS="bad"
        fi
        NUM_BAD_MACHINES=$((${NUM_BAD_MACHINES}-1))
      done
      if [ "${STATUS}" = "good" ]; then
        if [ "${TEMP}" = "" ]; then
          TEMP="${TEMP_REMOTE}"
        else
          TEMP="${TEMP} ${TEMP_REMOTE}"
        fi
      fi
      NUM_MACHINES=$((${NUM_MACHINES}-1))
    done
    REMOTE_MACHINES="${TEMP}"
  fi
fi

############################################################################

# create compressed tar ball containing all build log files
echo -n "Zipping build log files..."
# if at least one system failed to build the software then create a separate tar ball containing only build log
# files with errors
if [ "${ERROR_FLAG}" = "yes" ]; then
  ( cd "${LOG_DIR}" && tar -cf "${PROJECT_NAME}"-build-log-ERROR.tar *-build.log-ERROR && gzip "${PROJECT_NAME}"-build-log-ERROR.tar )
  # reset error flag
  ERROR_FLAG="no"
  ERROR_DETECTED="yes"
fi
# create a tar ball containing build log files for successful builds
( cd "${LOG_DIR}" && tar -cf "${PROJECT_NAME}"-build-log.tar *-build.log && gzip "${PROJECT_NAME}"-build-log.tar )
# delete build log files to save space
( cd "${LOG_DIR}" && rm -f *.log *.log-ERROR )
echo -e "[DONE]\n"

############################################################################
#
# run test code (examples)
#
############################################################################

# run examples on systems
echo "Testing ${PROJECT_NAME} on system(s):"
LIST_BAD_TESTS=""
NUM_MACHINES=`echo "${REMOTE_MACHINES}" | wc -w`
while [ $((${NUM_MACHINES})) -gt 0 ]; do
  TEMP_REMOTE=`echo "${REMOTE_MACHINES}" | cut -d' ' -f $((${NUM_MACHINES}))`
  echo -ne "\tTesting ${PROJECT_NAME} on ${TEMP_REMOTE}..."
  # if remote system then execute test commands via remote login command
  if [ ! "${TEMP_REMOTE}" = "${LOCAL_MACHINE}" ]; then
    TEMP="${TEMP_REMOTE}"
    # determine remote directory (check input file)
    find_remote_dir < "${INPUT_FILE}"
    REMOTE_DIR=`cat temp-dir.txt`
    rm -f temp-dir.txt

    NODE_TEMP_REMOTE="${TEMP_REMOTE}"

    # begin software test on remote node
    ${REMOTE_LOGIN_CMD} ${REMOTE_LOGIN_ARGS} -l ${REMOTE_USERNAME} ${NODE_TEMP_REMOTE} \
    " cd ${REMOTE_DIR}/${TEMP_REMOTE} && \
    ./${FIX_BASH} ${TIMEOUT_SCRIPT} no && \
    ./${FIX_BASH} ${TEST_SCRIPT} no && \
    ./${TEST_SCRIPT} ${REMOTE_DIR}/${TEMP_REMOTE}/${PROJECT_NAME} ${TIMEOUT_SCRIPT} ${REMOTE_USERNAME} ${PROJECT_NAME} ${MAX_TIME} ${REMOTE_LOGIN_CMD} ${REMOTE_LOGIN_ARGS}" &> "${LOG_DIR}/${TEMP_REMOTE}"-test.log
    echo -e "[DONE]\n"
  # commands necessary to test software on local system (remote command execution not required)
  else
    cd "${LOCAL_DIR}/${LOCAL_MACHINE}" && \
    ./"${FIX_BASH}" "${TIMEOUT_SCRIPT}" no && \
    ./"${FIX_BASH}" "${TEST_SCRIPT}" no && \
    ./"${TEST_SCRIPT}" "${LOCAL_DIR}/${LOCAL_MACHINE}/${PROJECT_NAME}" "${TIMEOUT_SCRIPT}" "${LOCAL_USERNAME}" "${PROJECT_NAME}" "${MAX_TIME}" "${REMOTE_LOGIN_CMD}" "${REMOTE_LOGIN_ARGS}" &> "${LOG_DIR}/${LOCAL_MACHINE}"-test.log
    echo -e "[DONE]\n"
  fi

  if [ ! "${TEMP_REMOTE}" = "${LOCAL_MACHINE}" ]; then
    # retieve test log files from remote system
    # actual output from the example programs is checked later
    echo -ne "\tRetrieving test log file(s) from ${TEMP_REMOTE}..."
    STATUS=`"${REMOTE_COPY_CMD}" "${REMOTE_COPY_ARGS}" "${REMOTE_USERNAME}"@"${TEMP_REMOTE}":"${REMOTE_DIR}/${TEMP_REMOTE}/${PROJECT_NAME}/${TEMP_REMOTE}"-examples-logs.tar.gz ${LOG_DIR}/examples/ 2>&1`
    CHECK_STATUS_0=`echo "${STATUS}" | fgrep " not"`
    CHECK_STATUS_1=`echo "${STATUS}" | fgrep "No such file or directory"`
    CHECK_STATUS_2=`echo "${STATUS}" | fgrep -i "warning: Executing scp1"`
    # Check if we did get the log file
    if test -f ${LOG_DIR}/examples/${TEMP_REMOTE}-examples-logs.tar.gz ; then
      IGNORE_RET_ERR="yes"
    else
      IGNORE_RET_ERR="no"
    fi
    if [ "${IGNORE_RET_ERR}" = "no" ]; then
      # Retry once if test log retrieval fails
      if [ "${CHECK_STATUS_0}" != "" -o "${CHECK_STATUS_1}" != "" -a "${CHECK_STATUS_2}" = "" ]; then
        sleep 30
        echo -e "[FAILED - retrying]"
        STATUS=""
        CHECK_STATUS_0=""
        CHECK_STATUS_1=""
        CHECK_STATUS_2=""
        echo -ne "\tRetrieving test log file(s) from ${TEMP_REMOTE}..."
        STATUS=`"${REMOTE_COPY_CMD}" "${REMOTE_COPY_ARGS}" "${REMOTE_USERNAME}"@"${TEMP_REMOTE}":"${REMOTE_DIR}/${TEMP_REMOTE}/${PROJECT_NAME}/${TEMP_REMOTE}"-examples-logs.tar.gz ${LOG_DIR}/examples/ 2>&1`
        CHECK_STATUS_0=`echo "${STATUS}" | fgrep " not"`
        CHECK_STATUS_1=`echo "${STATUS}" | fgrep "No such file or directory"`
        CHECK_STATUS_2=`echo "${STATUS}" | fgrep -i "warning: Executing scp1"`
      fi
      # if test log file cannot be retrieved then flag system name
      if [ "${CHECK_STATUS_0}" != "" -o "${CHECK_STATUS_1}" != "" -a "${CHECK_STATUS_2}" = "" ]; then
        if [ "${LIST_BAD_TESTS}" = "" ]; then
          LIST_BAD_TESTS="${TEMP_REMOTE}"
        else
          LIST_BAD_TESTS="${LIST_BAD_TESTS} ${TEMP_REMOTE}"
        fi
        # indicate if file retieval failed because system name could not be resolved
        if [ "${CHECK_STATUS_0}" != "" ]; then
          echo -e "[FAILED - unknown host]\n"
          echo -e "\tRetrieving test log file(s) from ${TEMP_REMOTE}...[FAILED - unknown host]" >> "${LOG_DIR}/${TEMP_REMOTE}"-test.log
        # indicate if file retrieval failed because test log file simply does not exist
        else
          echo -e "[FAILED - no log files found]\n"
          echo -e "\tRetrieving test log file(s) from ${TEMP_REMOTE}...[FAILED - no log files found]" >> "${LOG_DIR}/${TEMP_REMOTE}"-test.log
        fi
      # get test log files (remote copy)
      else
        "${REMOTE_COPY_CMD}" "${REMOTE_COPY_ARGS}" "${REMOTE_USERNAME}"@"${TEMP_REMOTE}":"${REMOTE_DIR}/${TEMP_REMOTE}/${PROJECT_NAME}/${TEMP_REMOTE}"-examples-logs.tar.gz "${LOG_DIR}"/examples 1> /dev/null 2> /dev/null
        echo -e "[DONE]\n"
      fi
    # get test log files (remote copy)
    else
      "${REMOTE_COPY_CMD}" "${REMOTE_COPY_ARGS}" "${REMOTE_USERNAME}"@"${TEMP_REMOTE}":"${REMOTE_DIR}/${TEMP_REMOTE}/${PROJECT_NAME}/${TEMP_REMOTE}"-examples-logs.tar.gz "${LOG_DIR}"/examples 1> /dev/null 2> /dev/null
      echo -e "[DONE]\n"
    fi
  # if local system just move the test log file to the appropriate directory
  else
    echo -ne "\tMoving test log file(s) to appropriate directory on local system..."
    STATUS=`cp -f "${LOCAL_DIR}/${LOCAL_MACHINE}/${PROJECT_NAME}/${LOCAL_MACHINE}"-examples-logs.tar.gz "${LOG_DIR}/examples" 2>&1`
    CHECK_STATUS_0=`echo "${STATUS}" | fgrep "No such file or directory"`
    # if test log files could not be moved (because does not exist) then flag system name
    if [ ! "${CHECK_STATUS_0}" = "" ]; then
      if [ "${LIST_BAD_TESTS}" = "" ]; then
        LIST_BAD_TESTS="${LOCAL_MACHINE}"
      else
        LIST_BAD_TESTS="${LIST_BAD_TESTS} ${LOCAL_MACHINE}"
      fi
      echo -e "[FAILED - no log files found]\n"
      echo -e "\tMoving test log file(s) to appropriate directory on local system...[FAILED - no log files found]" >> "${LOG_DIR}/${TEMP_REMOTE}"-test.log
    else
      cp -f "${LOCAL_DIR}/${LOCAL_MACHINE}/${PROJECT_NAME}/${LOCAL_MACHINE}"-examples-logs.tar.gz "${LOG_DIR}"/examples
      echo -e "[DONE]\n"
    fi
  fi
  
  NUM_MACHINES=$((${NUM_MACHINES} - 1))
done

############################################################################

# only check available test log files (those that could be retrieved)
if [ ! "${LIST_BAD_TESTS}" = "" ]; then
  echo -e "Skipping analysis of logs from the following system(s) because of file retrieval error(s): ${LIST_BAD_TESTS}\n"
  NUM_MACHINES=`echo "${REMOTE_MACHINES}" | wc -w`
  NUM_BAD_MACHINES=`echo "${LIST_BAD_TESTS}" | wc -w`
  # if errors were indicated on every system then e-mail error message and abort
  # stated another way, if no log files could be retrieved then abort because there is nothing to check
  if [ $((${NUM_BAD_MACHINES})) -eq $((${NUM_MACHINES})) ]; then
    echo -e "\a\n\tERROR: all systems failed to run tests for ${PROJECT_NAME}\n"
    echo -n "E-mailing exit status to ${MAIL_RECIPIENTS} before exiting..."
    MAIL_SUBJECT="Autotest Message (${PROJECT_NAME}): Error"
    echo -e "\nERROR: all systems failed to run tests for ${PROJECT_NAME}\n" >> "${LOG_DIR}/${MAIL_FILE}"
    "${MAIL_CMD}" -s "${MAIL_SUBJECT}" "${MAIL_RECIPIENTS}" < "${LOG_DIR}/${MAIL_FILE}"
    echo -e "[DONE]\n"
    exit 0
  fi
  # if log files could not be retrieved from at least one system then edit system list
  if [ $((${NUM_BAD_MACHINES})) -gt 0 ]; then
    TEMP=""
    while [ $((${NUM_MACHINES})) -gt 0 ]; do
      TEMP_REMOTE=`echo "${REMOTE_MACHINES}" | cut -d' ' -f$((${NUM_MACHINES}))`
      NUM_BAD_MACHINES=`echo "${LIST_BAD_TESTS}" | wc -w`
      STATUS="good"
      while [ $((${NUM_BAD_MACHINES})) -gt 0 ]; do
        TEMP_BAD=`echo "${LIST_BAD_TESTS}" | cut -d' ' -f$((${NUM_BAD_MACHINES}))`
        if [ "${TEMP_REMOTE}" = "${TEMP_BAD}" ]; then
          STATUS="bad"
        fi
        NUM_BAD_MACHINES=$((${NUM_BAD_MACHINES}-1))
      done
      if [ "${STATUS}" = "good" ]; then
        if [ "${TEMP}" = "" ]; then
          TEMP="${TEMP_REMOTE}"
        else
          TEMP="${TEMP} ${TEMP_REMOTE}"
        fi
      fi
      NUM_MACHINES=$((${NUM_MACHINES}-1))
    done
    REMOTE_MACHINES="${TEMP}"
  fi
fi

############################################################################

echo -n "Extracting test log files..."
# extract available test log files for analysis
NUM_MACHINES=`echo "${REMOTE_MACHINES}" | wc -w`
while [ $((${NUM_MACHINES})) -gt 0 ]; do
  TEMP_MACHINE=`echo "${REMOTE_MACHINES}" | cut -d' ' -f$((${NUM_MACHINES}))`
  ( cd "${LOG_DIR}/examples" && gunzip "${TEMP_MACHINE}"-examples-logs.tar.gz && tar -xf "${TEMP_MACHINE}"-examples-logs.tar )
  NUM_MACHINES=$((${NUM_MACHINES} - 1))
done
# move test log files to directory where example log files are kept
mv -f "${LOG_DIR}"/*.log "${LOG_DIR}/examples"
echo -e "[DONE]\n"

############################################################################
#
# search for error messages
#
############################################################################

# tag log files containing error messages with ERROR keyword
tag_error_files()
{
  while read TEMP_ERROR_FILE; do
    if [ ! "${TEMP_ERROR_FILE}" = "" ]; then
      if [ -r "${TEMP_ERROR_FILE}" ]; then
	# append informative header and log file to e-mail message
        TAG_SYSTEM=`basename "${TEMP_ERROR_FILE}" | cut -d'-' -f1`
        TAG_MODULE=`basename "${TEMP_ERROR_FILE}" | cut -d'-' -f2`
        TAG_FILE=`basename "${TEMP_ERROR_FILE}" | cut -d '-' -f3 | cut -d'.' -f1`
        mv -f "${TEMP_ERROR_FILE}" "${TEMP_ERROR_FILE}"-ERROR
        echo "#####################################" >> "${LOG_DIR}/${MAIL_FILE}"
        echo "#" >> "${LOG_DIR}/${MAIL_FILE}"
        echo "#  TEST ERROR" >> "${LOG_DIR}/${MAIL_FILE}"
        echo "#" >> "${LOG_DIR}/${MAIL_FILE}"
        echo "#  SYSTEM: ${TAG_SYSTEM}" >> "${LOG_DIR}/${MAIL_FILE}"
        echo "#  MODULE: ${TAG_MODULE}" >> "${LOG_DIR}/${MAIL_FILE}"
        echo "#  FILE: ${TAG_FILE}" >> "${LOG_DIR}/${MAIL_FILE}"
        echo "#  LOG FILE: ${TEMP_ERROR_FILE}" >> "${LOG_DIR}/${MAIL_FILE}"
        echo "#" >> "${LOG_DIR}/${MAIL_FILE}"
        echo "#####################################" >> "${LOG_DIR}/${MAIL_FILE}"
        cat "${TEMP_ERROR_FILE}"-ERROR >> "${LOG_DIR}/${MAIL_FILE}"
      fi
    fi
  done
}

# check log files for error messages
check_error_list()
{
  KEYWORD_FLAG="no"
  ALL_STATUS=""
  # perform case-insensitive keyword search using terms found in file indicated by KEYWORD_ERROR_FILE variable
  TEMP_STATUS=`fgrep -i -f "${SCRIPT_DIR}/${KEYWORD_ERROR_FILE}" "${LOG_DIR}/examples/${TEMP_REMOTE}"-*.log`
  if [ ! "${TEMP_STATUS}" = "" ]; then
    KEYWORD_FLAG="yes"
    # create list of all log files associated with current system which contain error messages
    if [ ! "${ALL_STATUS}" = "" ]; then
      ALL_STATUS="${ALL_STATUS}
${TEMP_STATUS}"
    else
      ALL_STATUS="${TEMP_STATUS}"
    fi
  fi
  echo "${KEYWORD_FLAG}" > keyword-temp.txt
}

echo "Checking test log file(s) for possible error messages:"
TEST_ERROR_LIST=""
NUM_MACHINES=`echo "${REMOTE_MACHINES}" | wc -w`
# check log files for error messages
while [ $((${NUM_MACHINES})) -gt 0 ]; do
  TEMP_REMOTE=`echo "${REMOTE_MACHINES}" | cut -d' ' -f$((${NUM_MACHINES}))`
  echo -ne "\tChecking test log file(s) from ${TEMP_REMOTE}..."
  # check log files associated with system TEMP_REMOTE
  check_error_list
  STATUS=`cat keyword-temp.txt`
  rm -f keyword-temp.txt
  if [ "${STATUS}" = "yes" ]; then
    ERROR_FLAG="yes"
    echo -e "[ERROR]\n"
    # if system log file contains an error message then update list of systems with errors
    if [ "${TEST_ERROR_LIST}" = "" ]; then
      TEST_ERROR_LIST="${TEMP_REMOTE}"
    else
      TEST_ERROR_LIST="${TEST_ERROR_LIST} ${TEMP_REMOTE}"
    fi
    # tag each log file containing an error message
    TEST_FILE_LIST=`echo "${ALL_STATUS}" | cut -d':' -f1`
    echo "${TEST_FILE_LIST}" | tag_error_files
  else
    echo -e "[DONE]\n"
  fi
  NUM_MACHINES=$((${NUM_MACHINES} - 1))
done
if [ ! "${TEST_ERROR_LIST}" = "" ]; then
  echo -e "Possible test error(s) on the following system(s): ${TEST_ERROR_LIST}\n"
fi
NUM_BAD_MACHINES=`echo "${TEST_ERROR_LIST}" | wc -w`
NUM_MACHINES=`echo "${REMOTE_MACHINES}" | wc -w`

############################################################################

# tar and compress all log files
if [ "${ERROR_FLAG}" = "no" ]; then
  echo "No errors were detected :-)" >> "${LOG_DIR}/${MAIL_FILE}"
fi
echo -n "Zipping test log files..."
( cd "${LOG_DIR}"/examples && rm -f *.tar )
# tar and compress files with errors separately
if [ "${ERROR_FLAG}" = "yes" ]; then
  ( cd "${LOG_DIR}"/examples && tar -cf "${PROJECT_NAME}"-test-logs-ERROR.tar *.log-ERROR && gzip "${PROJECT_NAME}"-test-logs-ERROR.tar && mv -f "${PROJECT_NAME}"-test-logs-ERROR.tar.gz ../ )
  ERROR_DETECTED="yes"
fi
# tar and compress normal log files
( cd "${LOG_DIR}"/examples && tar -cf "${PROJECT_NAME}"-test-logs.tar *.log && gzip "${PROJECT_NAME}"-test-logs.tar && mv -f "${PROJECT_NAME}"-test-logs.tar.gz ../ )
echo -e "[DONE]\n"

rm -Rf "${LOG_DIR}"/examples

############################################################################
#
# rotate log files
#
############################################################################

# archive all log files and miscellaneous files related to current autotest run
( cd "${LOG_DIR}" && tar -cf "${PROJECT_NAME}"-logs.tar *.tar.gz "${PROJECT_NAME}".info )
# delete duplicate copies to save space
( cd "${LOG_DIR}" && rm -f *.tar.gz )
# create archive directory if does not exist
if [ ! -d "${LOG_DIR}/archive" ]; then
  mkdir "${LOG_DIR}/archive"
fi
# rotate log files
LAST_FILE=`ls -t "${LOG_DIR}"/archive/"${PROJECT_NAME}"-logs.tar-[0-9] 2> /dev/null | head -n1`
if [ ! "${LAST_FILE}" = "" ]; then
  LAST_FILE=`basename "${LAST_FILE}"`
  LAST_FILE=`echo "${LAST_FILE}" | cut -d'-' -f3`
else
  LAST_FILE=$((${NUM_LOG_FILES} - 1))
fi
LOG_FILE=$((${LAST_FILE} + 1))
LOG_FILE=$((${LOG_FILE} % ${NUM_LOG_FILES}))
rm -f "${LOG_DIR}"/archive/"${PROJECT_NAME}"-logs.tar-"${LOG_FILE}"
mv -f "${LOG_DIR}"/"${PROJECT_NAME}"-logs.tar "${LOG_DIR}"/archive/"${PROJECT_NAME}"-logs.tar-"${LOG_FILE}"
echo -e "Log file: ${LOG_DIR}/archive/${PROJECT_NAME}-logs.tar-${LOG_FILE}\n"

############################################################################
#
# e-mail results
#
############################################################################

# send e-mail only if errors were detected
if [ "${ERROR_DETECTED}" = "yes" ]; then
  echo -n "E-mailing error message(s) to: ${MAIL_RECIPIENTS}..."
  MAIL_SUBJECT="Autotest Message (${PROJECT_NAME}): Error"
  STATUS=`"${MAIL_CMD}" -s "${MAIL_SUBJECT}" "${MAIL_RECIPIENTS}" < "${LOG_DIR}/${MAIL_FILE}"`
  if [ "${STATUS}" = "" ]; then
    rm -f "${LOG_DIR}/${MAIL_FILE}"
    echo -e "[DONE]\n"
  else
    mv -f "${LOG_DIR}/${MAIL_FILE}" "${LOG_DIR}/archive/${MAIL_FILE}-${LOG_FILE}"
    echo -e "[FAILED]\n"
    echo -e "Error Message: ${LOG_DIR}/archive/${MAIL_FILE}-${LOG_FILE}\n"
  fi
else
  echo -n "E-mailing test completion notification to: ${MAIL_RECIPIENTS}..."
  MAIL_SUBJECT="Autotest Message (${PROJECT_NAME}): OK"
  STATUS=`"${MAIL_CMD}" -s "${MAIL_SUBJECT}" "${MAIL_RECIPIENTS}" < "${LOG_DIR}/${MAIL_FILE}"`
  if [ "${STATUS}" = "" ]; then
    rm -f "${LOG_DIR}/${MAIL_FILE}"
    echo -e "[DONE]\n"
  else
    mv -f "${LOG_DIR}/${MAIL_FILE}" "${LOG_DIR}/archive/${MAIL_FILE}-${LOG_FILE}"
    echo -e "[FAILED]\n"
    echo -e "Status Message: ${LOG_DIR}/archive/${MAIL_FILE}-${LOG_FILE}\n"
  fi
fi

echo -e "Finished :-)\n"

exit 0
