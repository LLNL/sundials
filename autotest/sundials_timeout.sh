#!/bin/bash


############################################################################
# $Revision$
# $Date$
############################################################################
#
# Filename: sundials_timeout.sh
# Programmer: Aaron Collier @ LLNL
#
############################################################################


# get system name
get_system_name()
{
  CONTINUE="yes"
  NUM_CHARS=`echo "${LOCAL_SYSTEM}" | wc -m`
  NUM_CHARS=$((${NUM_CHARS} - 1))

  # eliminate node number digit by digit
  while [ "${CONTINUE}" != "no" ]; do
    TEST=`echo "${LOCAL_SYSTEM}" | cut -c$((${NUM_CHARS}))`
    TEST=${TEST%[0-9]}
    if [ "${TEST}" = "" ]; then
      LOCAL_SYSTEM=${LOCAL_SYSTEM%[0-9]}
    else
      CONTINUE="no"
    fi
    NUM_CHARS=$((${NUM_CHARS} - 1))
  done

  echo "${LOCAL_SYSTEM}" > temp.txt
}

# truncate program names to account for Solaris weirdness
fix_sun()
{
  TEMP_LIST=""

  while read IN_LINE; do
    if [ "${IN_LINE}" != "" ]; then
      IN_LINE=`echo "${IN_LINE}" | cut -c1-8`
      if [ "${TEMP_LIST}" = "" ]; then
        TEMP_LIST="${IN_LINE}"
      else
        TEMP_LIST="${TEMP_LIST}
${IN_LINE}"
      fi
    fi
  done

  rm -f "${IN_LIST}"
  echo "${TEMP_LIST}" > "${IN_LIST}"
}

# get IPC (Inter-process Communication) resource id information
get_ipc_id()
{
  IN_FILE="$1"
  TEMP_ID=""

  while read IN_KEY IN_ID IN_OWNER IN_PERMS IN_BYTES IN_NATTCH IN_STATUS; do
    if [ "${MACHINE_TYPE}" = "Linux" ]; then
      if [ "${IN_OWNER}" = "${IN_REMOTE_USER}" ]; then
        if [ "${TEMP_ID}" = "" ]; then
          TEMP_ID="${IN_ID}"
        else
          TEMP_ID="${TEMP_ID} ${IN_ID}"
        fi
      fi
    else
      if [ "${IN_BYTES}" = "${IN_REMOTE_USER}" ]; then
        if [ "${TEMP_ID}" = "" ]; then
          TEMP_ID="${IN_ID}"
        else
          TEMP_ID="${TEMP_ID} ${IN_ID}"
        fi
      fi
    fi
  done

  echo "${TEMP_ID}" > "${IN_FILE}"
}

# release appropriate IPC resources
remove_ipc_id()
{
  BEFORE_ID=`cat $1`
  AFTER_ID=`cat $2`

  if [ "${AFTER_ID}" = "" ]; then
    return
  fi

  if [ "${BEFORE_ID}" = "" ]; then
    BEFORE_ID="blank"
  fi

  if [ "${MACHINE_TYPE}" != "Linux" ]; then
    if [ "${IPC_TYPE}" = "shm" ]; then
      IPC_TYPE="-m"
    elif [ "${IPC_TYPE}" = "sem" ]; then
      IPC_TYPE="-s"
    else
      IPC_TYPE="-q"
    fi
  fi

  NUM_ID=`echo "${AFTER_ID}" | wc -w`
  while [ $((${NUM_ID})) -gt 0 ]; do
    CURRENT_ID=`echo "${AFTER_ID}" | cut -d' ' -f$((${NUM_ID}))`
    STATUS=`echo "${BEFORE_ID}" | fgrep "${CURRENT_ID}"`
    if [ "${STATUS}" = "" ]; then
      ipcrm ${IPC_TYPE} ${CURRENT_ID} 1> /dev/null 2> /dev/null
    fi
  NUM_ID=$((${NUM_ID} - 1))
  done
}

# get process id
get_proc_pid()
{
  TEMP_PID_LIST=""

  while read IN_PID IN_IGNORE; do
    if [ "${IN_PID}" != "" -a "${IN_IGNORE}" != "" ]; then
      if [ "${TEMP_PID_LIST}" = "" ]; then
        TEMP_PID_LIST="${IN_PID}"
      else
        TEMP_PID_LIST="${TEMP_PID_LIST} ${IN_PID}"
      fi
    fi
  done

  echo "${TEMP_PID_LIST}" > temp-pid.txt
}

# determine total running time for given example
# strictly speaking, finds maximum running time for all processes in given group (all threads/sub-processes related to
# a particular process)
find_max_time()
{
  MAX_TIME=""

  while read IN_PROC_LINE; do
    if [ "${IN_PROC_LINE}" != "" ]; then
      if [ "${MACHINE_TYPE}" = "Linux" ]; then
        IN_MIN=`echo "${IN_PROC_LINE}" | cut -d':' -f2`
      else
        NUM_CHAR_SKIP="1"
        CURRENT_CHAR=""
        while [ "${CURRENT_CHAR}" != ":" ]; do
          CURRENT_CHAR=`echo "${IN_PROC_LINE}" | cut -c$((${NUM_CHAR_SKIP}))`
          NUM_CHAR_SKIP=$((${NUM_CHAR_SKIP} + 1))
        done
        IN_MIN=`echo "${IN_PROC_LINE}" | cut -c$((${NUM_CHAR_SKIP} - 2)),$((${NUM_CHAR_SKIP} - 3))`
      fi
      if [ "${MAX_TIME}" = "" ]; then
        MAX_TIME="${IN_MIN}"
      elif [ $((${IN_MIN})) -gt $((${MAX_TIME})) ]; then
        MAX_TIME="${IN_MIN}" 
      fi
    fi
  done

  echo "${MAX_TIME}" > temp-time.txt
}

# kill all user processes on remote nodes and release all reserved IPC resources
# only call this shell function if MPICH-MPI is being used
kill_node_procs()
{
  while read IN_NODE; do
    if [ "${IN_NODE}" != "" -a "${IN_NODE}" != "${LOCAL_SYSTEM}" ]; then
      "${IN_REMOTE_LOGIN_CMD}" "${IN_REMOTE_LOGIN_ARGS}" -l "${IN_REMOTE_USER}" "${IN_NODE}" " kill -9 -1 ; ${IN_MPI_DIR}/sbin/cleanipcs " 1> /dev/null 2> /dev/null
    fi
  done
}

# systematically check all processes related to examples
# determine if maximum allowable running time has been exceeded,
# or if process has entered zombie state (meaning total running time has
# remained unchanged for certain period of time)
check_procs()
{
  while read PROC_LIST; do
    if [ "${PROC_LIST}" != "" ]; then

      # Linux systems...
      if [ "${MACHINE_TYPE}" = "Linux" ]; then
        PROC_STATUS=`${PS_CMD} 2> /dev/null`
        PROC_STATUS=`echo "${PROC_STATUS}" | fgrep "${PROC_LIST}"`
        echo "${PS_IGNORE_KEYWORDS}" > temp-linux-list.txt
	# find current total running time for all example-related processes
        echo "${PROC_STATUS}" | find_max_time
        NUM_MINUTES=`cat temp-time.txt`
        rm -f temp-time.txt
	# if zombie process (determined by number of consecutive passes reporting total running time as zero minutes) then
        # set total running time to maximum allowable running time to force process termination
	# checks if process in question is related to currently running example (remember sequential execution)
	# process considered zombie if total running time remains zero for (max_run_time*2) passes
	# this odd behavior has only been seen on Linux systems running early versions of MPICH-MPI
        if [ $((${NUM_MINUTES})) -eq 0 ]; then
          if [ "${LAST_PROG}" = "" ]; then
            LAST_PROG=`ls -1 -rt ${IN_PROJECT_TITLE}/${LOCAL_SYSTEM}-*-*.log | tail -n1 | cut -d'-' -f3 | cut -d'.' -f1`
          else
            LAST_PROG="${CURR_PROG}"
          fi
          CURR_PROG=`ls -1 -rt ${IN_PROJECT_TITLE}/${LOCAL_SYSTEM}-*-*.log | tail -n1 | cut -d'-' -f3 | cut -d'.' -f1`
          if [ "${CURR_PROG}" = "${LAST_PROG}" ]; then
            if [ "${PROC_LIST}" = "${CURR_PROG}" ]; then
              if [ "${NUM_PASSES}" = "" ]; then
                NUM_PASSES="1"
              else
                if [ $((${NUM_PASSES})) -ge $((${IN_SLEEP} * 2)) ]; then
                  NUM_MINUTES=$((${IN_SLEEP}))
                  NUM_PASSES=""
                else
                  NUM_PASSES=$((${NUM_PASSES} + 1))
                fi
              fi
            fi
          else
            NUM_PASSES=""
          fi
        fi
        echo "${PROC_LIST}" >> temp-linux-list.txt
        if [ "${NUM_MINUTES}" != "" ]; then
	  PROC_STATUS=`${PS_CMD} 2> /dev/null | fgrep -f temp-linux-list.txt`
        fi
        rm -f temp-linux-list.txt
        sync

      # Compaq systems...
      elif [ "${MACHINE_TYPE}" = "OSF1" ]; then
        PROC_STATUS=`${PS_CMD} 2> /dev/null`
        PROC_STATUS=`echo "${PROC_STATUS}" | fgrep "${PROC_LIST}"`
        TEMP_PROC_STATUS=`echo "${PROC_STATUS}" | fgrep -v -e "${PS_IGNORE_KEYWORDS}"`
        echo "${TEMP_PROC_STATUS}" | find_max_time
        NUM_MINUTES=`cat temp-time.txt`
        rm -f temp-time.txt
        if [ "${NUM_MINUTES}" != "" ]; then
          PROC_STATUS=`${PS_CMD} 2> /dev/null`
          TEMP_PS_IGNORE_KEYWORDS="${PS_IGNORE_KEYWORDS}
${PROC_LIST}"
          PROC_STATUS=`echo "${PROC_STATUS}" | fgrep -e "${TEMP_PS_IGNORE_KEYWORDS}"`
        fi

      # SUN systems...
      elif [ "${MACHINE_TYPE}" = "SunOS" ]; then
        PROC_STATUS=`${PS_CMD} 2> /dev/null`
        PROC_STATUS=`echo "${PROC_STATUS}" | fgrep "${PROC_LIST}"`
        echo "${PS_IGNORE_KEYWORDS}" > temp-sun-list.txt
        echo "${PROC_STATUS}" | find_max_time
        NUM_MINUTES=`cat temp-time.txt`
        rm -f temp-time.txt
        echo "${PROC_LIST}" >> temp-sun-list.txt
        if [ "${NUM_MINUTES}" != "" ]; then
          PROC_STATUS=`${PS_CMD} 2> /dev/null | fgrep -f temp-sun-list.txt`
        fi
        rm -f temp-sun-list.txt
        sync
      fi

      if [ "${NUM_MINUTES}" != "" ]; then
        # kill process if has exceeded maximum allowable running time
        if [ $((${NUM_MINUTES})) -ge $((${IN_SLEEP})) ]; then
          # get process id's of all related processes
          echo "${PROC_STATUS}" | get_proc_pid
          KILL_LIST=`cat temp-pid.txt`
          rm -f temp-pid.txt
	  # issue warning message if exceeded maximum running time
          if [ "${WARNING_MESSAGE}" = "yes" ]; then
            echo "WARNING: ${PROC_LIST} exceeded maximum allotted execution time - ${IN_SLEEP} minute(s)"
	  # issue message if process killed during post-script clean-up
          elif [ "${WARNING_MESSAGE}" = "no" ]; then
            echo "NOTE: lingering ${PROC_LIST} processes were found"
          fi
          if [ "${MACHINE_TYPE}" != "OSF1" ]; then
            sync
          fi
	  # kill relevant processes
          kill -9 ${KILL_LIST} 1> /dev/null 2> /dev/null
        fi
      fi

    fi

  done
}

IN_LIST="$1"
IN_SLEEP="$2"
IN_SCRIPT="$3"
IN_USE_LAM_MPI="$4"
IN_REMOTE_USER="$5"
IN_PROJECT_TITLE="$6"
IN_REMOTE_LOGIN_CMD="$7"
IN_REMOTE_LOGIN_ARGS="$8"
IN_MPI_DIR="$9"

# determine machine type (operating system)
MACHINE_TYPE=`uname -s`
# determine which flags should be used with the ps command
# certain keywords are ignored later when determining total running time for a given process
if [ "${MACHINE_TYPE}" = "SunOS" ]; then
  PS_CMD="ps -U ${IN_REMOTE_USER}"
  PS_IGNORE_KEYWORDS="mpirun"
elif [ "${MACHINE_TYPE}" = "OSF1" ]; then
  PS_CMD="ps -U ${IN_REMOTE_USER}"
  PS_IGNORE_KEYWORDS="dmpirun
mpihost"
elif [ "${MACHINE_TYPE}" = "Linux" ]; then
  PS_CMD="ps --no-headers --User ${IN_REMOTE_USER}"
  PS_IGNORE_KEYWORDS="mpirun"
fi

# get name of local system
if [ "${MACHINE_TYPE}" = "SunOS" ]; then
  LOCAL_SYSTEM=`hostname`
else
  LOCAL_SYSTEM=`hostname -s`
fi
# if CASC Linux system then keep node number
STATUS=`echo "${LOCAL_SYSTEM}" | fgrep "tux"`
# if LC cluster system then remove node number from name
# clusters often have a login node pool (multiple nodes)
# want cluster name not node name
if [ "${STATUS}" = "" ]; then
  echo "${LOCAL_SYSTEM}" | get_system_name
  LOCAL_SYSTEM=`cat temp.txt`
  rm -f temp.txt
fi

# if not using LAM-MPI then IPC resources may not be released if a program aborts (returns with non-zero exit status)
# MPICH-MPI has a facility to clean IPC's, but will delete all IPC's associated with the current user
# this script will do slightly better - will take snapshots of reserved IPC resources before and after the test
# script has been executed and will compare the two snapshots
# any reserved IPC resources appearing in the post-test snapshot but not appearing in the pre-test snapshot will
# be freed
if [ "${IN_USE_LAM_MPI}" = "no" ]; then
  ipcs -m > temp-shm.txt
  ipcs -s > temp-sem.txt
  ipcs -q > temp-msg.txt
  get_ipc_id shm.var-before < temp-shm.txt
  get_ipc_id sem.var-before < temp-sem.txt
  get_ipc_id msg.var-before < temp-msg.txt
  rm -f temp-shm.txt temp-sem.txt temp-msg.txt
  if [ "${MACHINE_TYPE}" != "OSF1" ]; then
    sync
  fi
fi

WARNING_MESSAGE="yes"
NUM_PASSES=""
LAST_PROG=""

# Solaris systems limit process information displayed by the ps command to only eight characters, so
# all of the example executable names need to be truncated accordingly
# if not, then a search for an executable name of length greater than eight characters in the ps command
# output will always fail and the process will not be terminated if it exceeds the maximum allowable running time
if [ "${MACHINE_TYPE}" = "SunOS" ]; then
  fix_sun < "${IN_LIST}"
fi

# wait until at least one example process has started before beginning zombie process loop
# the file named script_running.info is created when the example program execution script generated by the test script
# begins execution
while [ ! -e script_running.info ]; do
  sleep 2
done

# start monitoring processes (check for zombies)
STATUS=`${PS_CMD} 2> /dev/null`
STATUS_0=`echo "${STATUS}" | fgrep "${IN_SCRIPT}"`
STATUS_1=`echo "${STATUS}" | fgrep "${PS_IGNORE_KEYWORDS}"`
while [ "${STATUS_0}" != "" -o "${STATUS_1}" != "" ]; do
  sleep 30
  check_procs < "${IN_LIST}"
  STATUS=`${PS_CMD} 2> /dev/null`
  if [ "${MACHINE_TYPE}" = "OSF1" ]; then
    STATUS_0=`echo "${STATUS}" | fgrep "${IN_SCRIPT}" | fgrep -v "${IN_LIST}"`
    STATUS_1=`echo "${STATUS}" | fgrep -e "${PS_IGNORE_KEYWORDS}"`
  else
    STATUS_0=`echo "${STATUS}" | fgrep "${IN_SCRIPT}"`
    STATUS_1=`echo "${STATUS}" | fgrep "${PS_IGNORE_KEYWORDS}"`
  fi
done

sleep 10

# clean-up all remaining processes (defunct processes) and issue a warning message rather than
# an error message
if [ "${IN_USE_LAM_MPI}" = "no" ]; then
  WARNING_MESSAGE="no"
  sleep 25
  # changing maximum running time to zero minutes causes all remaining processes to be terminated
  IN_SLEEP="0"
  check_procs < "${IN_LIST}"
fi

rm -f "${IN_LIST}"

# check if using MPICH-MPI
if [ -x "${IN_MPI_DIR}"/sbin/cleanipcs ]; then
  USE_MPICH_MPI="yes"
else
  USE_MPICH_MPI="no"
fi

# if not using LAM-MPI then release IPC resources
if [ "${IN_USE_LAM_MPI}" = "no" ]; then
  if [ "${USE_MPICH_MPI}" = "yes" ]; then
    if [ -r "${LOCAL_SYSTEM}.host" ]; then
      cat "${LOCAL_SYSTEM}".host | fgrep -v -e "#" | cut -d'.' -f1 > t3.txt
      kill_node_procs < t3.txt
      rm -f t3.txt
    fi
  fi
  ipcs -m > temp-shm.txt
  ipcs -s > temp-sem.txt
  ipcs -q > temp-msg.txt
  get_ipc_id shm.var-after < temp-shm.txt
  get_ipc_id sem.var-after < temp-sem.txt
  get_ipc_id msg.var-after < temp-msg.txt
  rm -f temp-shm.txt temp-sem.txt temp-msg.txt
  if [ "${MACHINE_TYPE}" != "OSF1" ]; then
    sync
  fi
  IPC_TYPE="shm"
  remove_ipc_id shm.var-before shm.var-after
  IPC_TYPE="sem"
  remove_ipc_id sem.var-before sem.var-after
  IPC_TYPE="msg"
  remove_ipc_id msg.var-before msg.var-after
  rm -f *.var-*
  if [ "${MACHINE_TYPE}" != "OSF1" ]; then
    sync
  fi
fi

exit 0
