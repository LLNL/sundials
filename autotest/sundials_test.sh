#!/bin/bash


############################################################################
# $Revision: 1.6 $
# $Date: 2006-07-25 01:06:28 $
############################################################################
#
# Filename: sundials_test.sh
# Programmer: Aaron Collier @ LLNL
#
############################################################################


# clean-up list of compiled examples by removing object file extension (*.o)
get_example_list()
{
  TEMP_EXAMPLE_LIST=""

  while read IN_EXAMPLE_FILE; do
    if [ "${IN_EXAMPLE_FILE}" != "" ]; then
      # remove *.o file extension
      IN_EXAMPLE_FILE=`basename "${IN_EXAMPLE_FILE}" .o`
      # remove '-updated' suffix from names of FCMIX examples
      TEMP_A=`echo "${IN_EXAMPLE_FILE}" | fgrep "updated"`
      if [ "${TEMP_A}" != "" ]; then
        IN_EXAMPLE_FILE=`echo "${IN_EXAMPLE_FILE}" | cut -d'-' -f1`
      fi
      # update list of executables
      if [ "${TEMP_EXAMPLE_LIST}" = "" ]; then
        TEMP_EXAMPLE_LIST="${IN_EXAMPLE_FILE}"
      else
        TEMP_EXAMPLE_LIST="${TEMP_EXAMPLE_LIST} ${IN_EXAMPLE_FILE}"
      fi
    fi
  done

  # return modified program list
  echo "${TEMP_EXAMPLE_LIST}" > temp.txt
}

# generate script to execute all examples found
run_examples()
{
  EXAMPLES_DIR="$1"
  PSUB_BATCH_LIST=""

  # generate list of all compiled examples by searching for object files
  EXAMPLE_FILES=`ls -1 ${EXAMPLES_DIR}/*.o 2>&1`
  STATUS=`echo "${EXAMPLE_FILES}" | fgrep -i "no"`
  # continue if current directory contains compiled examples
  if [ "${EXAMPLE_FILES}" != "" -a "${STATUS}" = "" ]; then
    # remove object file extension (*.o) to derive executable name
    echo "${EXAMPLE_FILES}" | get_example_list
    EXAMPLE_FILES=`cat temp.txt`
    rm -f temp.txt
    NUM_EXAMPLES=`echo "${EXAMPLE_FILES}" | wc -w`
    while [ $((${NUM_EXAMPLES})) -gt 0 ]; do
      TEMP_EXAMPLE_FILE=`echo "${EXAMPLE_FILES}" | cut -d' ' -f$((${NUM_EXAMPLES}))`
      TEMP_EXAMPLE_FILE=`basename "${TEMP_EXAMPLE_FILE}" .o`
      EXEC_TYPE=`basename "${EXAMPLES_DIR}"`
      # determine if serial or parallel example
      if [ "${EXEC_TYPE}" = "serial" -o "${EXEC_TYPE}" = "fcmix_serial" ]; then
        STATUS="OK"
	# if not using PSUB then generate execution script rather than job script
        if [ "${USE_PSUB}" = "no" ]; then
	  # examples with sensitivity options must be given command-line options to run properly
	  # execute each such example with all possible combinations of options
          if [ "${TEMP_MODULE}" = "cvodes" ]; then
            if [ "${TEMP_EXAMPLE_FILE}" = "cvfdx" -o "${TEMP_EXAMPLE_FILE}" = "cvfkx" -o "${TEMP_EXAMPLE_FILE}" = "cvfnx" -o "${TEMP_EXAMPLE_FILE}" = "pvfkt" ]; then
              echo "echo \"Running ${EXAMPLES_DIR}/${TEMP_EXAMPLE_FILE} -nosensi...\"" >> "${SCRIPT_NAME}"
              echo "${TEMP_EXAMPLE_FILE}" >> "${KILL_INFO}"
              echo "cd ${EXAMPLES_DIR} && ./${TEMP_EXAMPLE_FILE} -nosensi &> a.log" >> "${SCRIPT_NAME}"
              echo "echo \"Running ${EXAMPLES_DIR}/${TEMP_EXAMPLE_FILE} -sensi sim t...\"" >> "${SCRIPT_NAME}"
              echo "cd ${EXAMPLES_DIR} && ./${TEMP_EXAMPLE_FILE} -sensi sim t  &> b.log" >> "${SCRIPT_NAME}"
              echo "echo \"Running ${EXAMPLES_DIR}/${TEMP_EXAMPLE_FILE} -sensi sim f...\"" >> "${SCRIPT_NAME}"
              echo "cd ${EXAMPLES_DIR} && ./${TEMP_EXAMPLE_FILE} -sensi sim f &> c.log" >> "${SCRIPT_NAME}"
              echo "echo \"Running ${EXAMPLES_DIR}/${TEMP_EXAMPLE_FILE} -sensi stg t...\"" >> "${SCRIPT_NAME}"
              echo "cd ${EXAMPLES_DIR} && ./${TEMP_EXAMPLE_FILE} -sensi stg t &> d.log" >> "${SCRIPT_NAME}"
              echo "echo \"Running ${EXAMPLES_DIR}/${TEMP_EXAMPLE_FILE} -sensi stg f...\"" >> "${SCRIPT_NAME}"
              echo "cd ${EXAMPLES_DIR} && ./${TEMP_EXAMPLE_FILE} -sensi stg f &> e.log" >> "${SCRIPT_NAME}"
              echo "echo \"Running ${EXAMPLES_DIR}/${TEMP_EXAMPLE_FILE} -sensi stg1 t...\"" >> "${SCRIPT_NAME}"
              echo "cd ${EXAMPLES_DIR} && ./${TEMP_EXAMPLE_FILE} -sensi stg1 t &> f.log" >> "${SCRIPT_NAME}"
              echo "echo \"Running ${EXAMPLES_DIR}/${TEMP_EXAMPLE_FILE} -sensi stg1 f...\"" >> "${SCRIPT_NAME}"
              echo "cd ${EXAMPLES_DIR} && ./${TEMP_EXAMPLE_FILE} -sensi stg1 f &> g.log" >> "${SCRIPT_NAME}"
              if [ "${MACHINE_TYPE}" != "OSF1" ]; then
                echo "sync" >> "${SCRIPT_NAME}"
              fi
	      # combine separate output files (each corresponding to a single combination of options) into a single output file
              echo "cd ${EXAMPLES_DIR} && cat a.log b.log c.log d.log e.log f.log g.log > ${BASE_DIR}/${LOCAL_MACHINE}-${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.log" >> "${SCRIPT_NAME}"
	    # allow for normal example programs (meaing no command-line options) in subdirectories containing examples with sensitivity options
            else
              echo "echo \"Running ${EXAMPLES_DIR}/${TEMP_EXAMPLE_FILE}...\"" >> "${SCRIPT_NAME}"
              echo "${TEMP_EXAMPLE_FILE}" >> "${KILL_INFO}"
              echo "cd ${EXAMPLES_DIR} && ./${TEMP_EXAMPLE_FILE} &> ${BASE_DIR}/${LOCAL_MACHINE}-${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.log" >> "${SCRIPT_NAME}"
            fi
	  # if normal example program then just add execution command to script file
          else
            echo "echo \"Running ${EXAMPLES_DIR}/${TEMP_EXAMPLE_FILE}...\"" >> "${SCRIPT_NAME}"
            echo "${TEMP_EXAMPLE_FILE}" >> "${KILL_INFO}"
            echo "cd ${EXAMPLES_DIR} && ./${TEMP_EXAMPLE_FILE} &> ${BASE_DIR}/${LOCAL_MACHINE}-${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.log" >> "${SCRIPT_NAME}"
          fi
# ----------------------------
	# if using PSUB then generate job submission script
	# each example program has a separate job script which is called by the master job submission script
        else
          echo "echo \"Submitting ${EXAMPLES_DIR}/${TEMP_EXAMPLE_FILE}...\"" >> "${SCRIPT_NAME}"
          echo "psub ${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job" >> "${SCRIPT_NAME}"
          echo "sleep 5" >> "${SCRIPT_NAME}"
          touch "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
	  # examples with sensitivity options must be given command-line options to run properly
	  # execute each such example with all possible combinations of options
          if [ "${TEMP_MODULE}" = "cvodes" ]; then
            if [ "${TEMP_EXAMPLE_FILE}" = "cvfdx" -o "${TEMP_EXAMPLE_FILE}" = "cvfkx" -o "${TEMP_EXAMPLE_FILE}" = "cvfnx" -o "${TEMP_EXAMPLE_FILE}" = "pvfkt" ]; then
              if [ "${LOCAL_MACHINE}" = "thunder" ]; then
                echo "#PSUB -ln 1" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
              else
                echo "#PSUB -ln 1 -g 1" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
              fi
              echo "#PSUB -c pbatch,${LOCAL_MACHINE}" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
              echo "#PSUB -s /usr/local/bin/bash" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
	      PSUB_WAIT_TIME=$((${WAIT_TIME} * 7))
              echo "#PSUB -tM ${PSUB_WAIT_TIME}" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
              if [ "${LOCAL_MACHINE}" = "thunder" ]; then
                echo "#PSUB -b casc" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
              fi
              echo "" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
              echo "cd ${EXAMPLES_DIR} && ./${TEMP_EXAMPLE_FILE} -nosensi &> a.log" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
              echo "cd ${EXAMPLES_DIR} && ./${TEMP_EXAMPLE_FILE} -sensi sim t &> b.log" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
              echo "cd ${EXAMPLES_DIR} && ./${TEMP_EXAMPLE_FILE} -sensi sim f &> c.log" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
              echo "cd ${EXAMPLES_DIR} && ./${TEMP_EXAMPLE_FILE} -sensi stg t &> d.log" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
              echo "cd ${EXAMPLES_DIR} && ./${TEMP_EXAMPLE_FILE} -sensi stg f &> e.log" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
              echo "cd ${EXAMPLES_DIR} && ./${TEMP_EXAMPLE_FILE} -sensi stg1 t &> f.log" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
              echo "cd ${EXAMPLES_DIR} && ./${TEMP_EXAMPLE_FILE} -sensi stg1 f &> g.log" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
              echo "" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
              echo "sync" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
	      # combine output files
              echo "cd ${EXAMPLES_DIR} && cat a.log b.log c.log d.log e.log f.log g.log > ${BASE_DIR}/${LOCAL_MACHINE}-${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.log" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
              echo "" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
              echo "exit 0" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
	    # allow for normal example programs (meaning no command-line options) in subdirectories containing examples with sensitivity options
            else
              if [ "${LOCAL_MACHINE}" = "thunder" ]; then
                echo "#PSUB -ln 1" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
              else
                echo "#PSUB -ln 1 -g 1" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
              fi
              echo "#PSUB -c pbatch,${LOCAL_MACHINE}" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
              echo "#PSUB -s /usr/local/bin/bash" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
	      PSUB_WAIT_TIME=$((${WAIT_TIME}))
              echo "#PSUB -tM ${PSUB_WAIT_TIME}" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
              if [ "${LOCAL_MACHINE}" = "thunder" ]; then
                echo "#PSUB -b casc" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
              fi
              echo "" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
              echo "cd ${EXAMPLES_DIR} && ./${TEMP_EXAMPLE_FILE} &> ${BASE_DIR}/${LOCAL_MACHINE}-${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.log" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
              echo "" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
              echo "exit 0" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
            fi
	  # if normal example program then just add execution commands to job script
          else
            if [ "${LOCAL_MACHINE}" = "thunder" ]; then
              echo "#PSUB -ln 1" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
            else
              echo "#PSUB -ln 1 -g 1" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
            fi
            echo "#PSUB -c pbatch,${LOCAL_MACHINE}" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
            echo "#PSUB -s /usr/local/bin/bash" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
	    PSUB_WAIT_TIME="${WAIT_TIME}"
            echo "#PSUB -tM ${WAIT_TIME}" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
            if [ "${LOCAL_MACHINE}" = "thunder" ]; then
              echo "#PSUB -b casc" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
            fi
            echo "" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
            echo "cd ${EXAMPLES_DIR} && ./${TEMP_EXAMPLE_FILE} &> ${BASE_DIR}/${LOCAL_MACHINE}-${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.log" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
            echo "" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
            echo "exit 0" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
          fi
          chmod +x "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
        fi
# ----------------------------
      # parallel examples have different execution commands
      elif [ "${EXEC_TYPE}" = "parallel" -o "${EXEC_TYPE}" = "fcmix_parallel" ]; then
        STATUS="OK"
	# if not using PSUB then generate executable script containing program execution commands
        if [ "${USE_PSUB}" = "no" ]; then
	  # example programs with sensitivity options require command-line arguments
          if [ "${TEMP_MODULE}" = "cvodes" ]; then
            if [ "${TEMP_EXAMPLE_FILE}" = "pvfkx" -o "${TEMP_EXAMPLE_FILE}" = "pvfnx" -o "${TEMP_EXAMPLE_FILE}" = "pvfkt" ]; then
              echo "echo \"Running ${EXAMPLES_DIR}/${TEMP_EXAMPLE_FILE} -nosensi...\"" >> "${SCRIPT_NAME}"
              echo "echo \"${TEMP_EXAMPLE_FILE}\" > temp-running.chk" >> "${SCRIPT_NAME}"
              echo "${TEMP_EXAMPLE_FILE}" >> "${KILL_INFO}"
              echo "cd ${EXAMPLES_DIR} && ${MPI_CMD} ${TEMP_EXAMPLE_FILE} -nosensi &> a.log" >> "${SCRIPT_NAME}"
              echo "echo \"Running ${EXAMPLES_DIR}/${TEMP_EXAMPLE_FILE} -sensi sim t...\"" >> "${SCRIPT_NAME}"
              echo "cd ${EXAMPLES_DIR} && ${MPI_CMD} ${TEMP_EXAMPLE_FILE} -sensi sim t  &> b.log" >> "${SCRIPT_NAME}"
              echo "echo \"Running ${EXAMPLES_DIR}/${TEMP_EXAMPLE_FILE} -sensi sim f...\"" >> "${SCRIPT_NAME}"
              echo "cd ${EXAMPLES_DIR} && ${MPI_CMD} ${TEMP_EXAMPLE_FILE} -sensi sim f &> c.log" >> "${SCRIPT_NAME}"
              echo "echo \"Running ${EXAMPLES_DIR}/${TEMP_EXAMPLE_FILE} -sensi stg t...\"" >> "${SCRIPT_NAME}"
              echo "cd ${EXAMPLES_DIR} && ${MPI_CMD} ${TEMP_EXAMPLE_FILE} -sensi stg t &> d.log" >> "${SCRIPT_NAME}"
              echo "echo \"Running ${EXAMPLES_DIR}/${TEMP_EXAMPLE_FILE} -sensi stg f...\"" >> "${SCRIPT_NAME}"
              echo "cd ${EXAMPLES_DIR} && ${MPI_CMD} ${TEMP_EXAMPLE_FILE} -sensi stg f &> e.log" >> "${SCRIPT_NAME}"
              echo "echo \"Running ${EXAMPLES_DIR}/${TEMP_EXAMPLE_FILE} -sensi stg1 t...\"" >> "${SCRIPT_NAME}"
              echo "cd ${EXAMPLES_DIR} && ${MPI_CMD} ${TEMP_EXAMPLE_FILE} -sensi stg1 t &> f.log" >> "${SCRIPT_NAME}"
              echo "echo \"Running ${EXAMPLES_DIR}/${TEMP_EXAMPLE_FILE} -sensi stg1 f...\"" >> "${SCRIPT_NAME}"
              echo "cd ${EXAMPLES_DIR} && ${MPI_CMD} ${TEMP_EXAMPLE_FILE} -sensi stg1 f &> g.log" >> "${SCRIPT_NAME}"
              if [ "${MACHINE_TYPE}" != "OSF1" ]; then
                echo "sync" >> "${SCRIPT_NAME}"
              fi
	      # combine output files
              echo "cd ${EXAMPLES_DIR} && cat a.log b.log c.log d.log e.log f.log g.log > ${BASE_DIR}/${LOCAL_MACHINE}-${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.log" >> "${SCRIPT_NAME}"
	    # if just a normal example in a subdirectory containing sensitivity examples then use normal execution command
            else
              echo "echo \"Running ${EXAMPLES_DIR}/${TEMP_EXAMPLE_FILE}...\"" >> "${SCRIPT_NAME}"
              echo "${TEMP_EXAMPLE_FILE}" >> "${KILL_INFO}"
              echo "cd ${EXAMPLES_DIR} && ${MPI_CMD} ${TEMP_EXAMPLE_FILE} &> ${BASE_DIR}/${LOCAL_MACHINE}-${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.log" >> "${SCRIPT_NAME}"
            fi
	  # if normal example then just use normal execution command
          else
            echo "echo \"Running ${EXAMPLES_DIR}/${TEMP_EXAMPLE_FILE}...\"" >> "${SCRIPT_NAME}"
            echo "${TEMP_EXAMPLE_FILE}" >> "${KILL_INFO}"
            echo "cd ${EXAMPLES_DIR} && ${MPI_CMD} ${TEMP_EXAMPLE_FILE} &> ${BASE_DIR}/${LOCAL_MACHINE}-${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.log" >> "${SCRIPT_NAME}"
          fi
# ----------------------------
	# if using PSUB then generate job scripts and master job submission script
        else
          echo "echo \"Submitting ${EXAMPLES_DIR}/${TEMP_EXAMPLE_FILE}...\"" >> "${SCRIPT_NAME}"
          echo "psub ${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job" >> "${SCRIPT_NAME}"
          echo "sleep 5" >> "${SCRIPT_NAME}"
          touch "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
	  # recognize sensitivity examples which require command-line options
          if [ "${TEMP_MODULE}" = "cvodes" ]; then
            if [ "${TEMP_EXAMPLE_FILE}" = "pvfkx" -o "${TEMP_EXAMPLE_FILE}" = "pvfnx" -o "${TEMP_EXAMPLE_FILE}" = "pvfkt" ]; then
              if [ "${LOCAL_MACHINE}" = "thunder" ]; then
                echo "#PSUB -ln 1" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
              else
                echo "#PSUB -ln 1 -g 4" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
              fi
              echo "#PSUB -c pbatch,${LOCAL_MACHINE}" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
              echo "#PSUB -s /usr/local/bin/bash" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
	      PSUB_WAIT_TIME=$((${WAIT_TIME} * 7))
              echo "#PSUB -tM ${PSUB_WAIT_TIME}" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
              if [ "${LOCAL_MACHINE}" = "thunder" ]; then
                echo "#PSUB -b casc" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
              fi
              echo "" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
              if [ "${LOCAL_MACHINE}" = "thunder" ]; then
                echo "cd ${EXAMPLES_DIR} && srun -n4 ${TEMP_EXAMPLE_FILE} -nosensi &> a.log" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
                echo "cd ${EXAMPLES_DIR} && srun -n4 ${TEMP_EXAMPLE_FILE} -sensi sim t &> b.log" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
                echo "cd ${EXAMPLES_DIR} && srun -n4 ${TEMP_EXAMPLE_FILE} -sensi sim f &> c.log" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
                echo "cd ${EXAMPLES_DIR} && srun -n4 ${TEMP_EXAMPLE_FILE} -sensi stg t &> d.log" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
                echo "cd ${EXAMPLES_DIR} && srun -n4 ${TEMP_EXAMPLE_FILE} -sensi stg f &> e.log" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
                echo "cd ${EXAMPLES_DIR} && srun -n4 ${TEMP_EXAMPLE_FILE} -sensi stg1 t &> f.log" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
                echo "cd ${EXAMPLES_DIR} && srun -n4 ${TEMP_EXAMPLE_FILE} -sensi stg1 f &> g.log" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
              else
                echo "cd ${EXAMPLES_DIR} && ./${TEMP_EXAMPLE_FILE} -nosensi &> a.log" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
                echo "cd ${EXAMPLES_DIR} && ./${TEMP_EXAMPLE_FILE} -sensi sim t &> b.log" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
                echo "cd ${EXAMPLES_DIR} && ./${TEMP_EXAMPLE_FILE} -sensi sim f &> c.log" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
                echo "cd ${EXAMPLES_DIR} && ./${TEMP_EXAMPLE_FILE} -sensi stg t &> d.log" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
                echo "cd ${EXAMPLES_DIR} && ./${TEMP_EXAMPLE_FILE} -sensi stg f &> e.log" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
                echo "cd ${EXAMPLES_DIR} && ./${TEMP_EXAMPLE_FILE} -sensi stg1 t &> f.log" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
                echo "cd ${EXAMPLES_DIR} && ./${TEMP_EXAMPLE_FILE} -sensi stg1 f &> g.log" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
              fi
              echo "" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
              echo "sync" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
              echo "cd ${EXAMPLES_DIR} && cat a.log b.log c.log d.log e.log f.log g.log > ${BASE_DIR}/${LOCAL_MACHINE}-${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.log" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
              echo "" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
              echo "exit 0" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
	    # normal examples do not require command-line options
            else
              if [ "${LOCAL_MACHINE}" = "thunder" ]; then
                echo "#PSUB -ln 1" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
              else
                echo "#PSUB -ln 1 -g 4" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
              fi
              echo "#PSUB -c pbatch,${LOCAL_MACHINE}" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
              echo "#PSUB -s /usr/local/bin/bash" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
              PSUB_WAIT_TIME="${WAIT_TIME}"
              echo "#PSUB -tM ${PSUB_WAIT_TIME}" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
              if [ "${LOCAL_MACHINE}" = "thunder" ]; then
                echo "#PSUB -b casc" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
              fi
              echo "" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
              if [ "${LOCAL_MACHINE}" = "thunder" ]; then
                echo "cd ${EXAMPLES_DIR} && srun -n4 ${TEMP_EXAMPLE_FILE} &> ${BASE_DIR}/${LOCAL_MACHINE}-${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.log" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
              else
                echo "cd ${EXAMPLES_DIR} && ./${TEMP_EXAMPLE_FILE} &> ${BASE_DIR}/${LOCAL_MACHINE}-${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.log" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
              fi
              echo "" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
              echo "exit 0" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
            fi
	  # normal examples do not require command-line options
          else
            if [ "${LOCAL_MACHINE}" = "thunder" ]; then
              echo "#PSUB -ln 1" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
            else
              echo "#PSUB -ln 1 -g 4" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
            fi
            echo "#PSUB -c pbatch,${LOCAL_MACHINE}" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
            echo "#PSUB -s /usr/local/bin/bash" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
	    PSUB_WAIT_TIME="${WAIT_TIME}"
            echo "#PSUB -tM ${PSUB_WAIT_TIME}" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
            if [ "${LOCAL_MACHINE}" = "thunder" ]; then
              echo "#PSUB -b casc" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
            fi
            echo "" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
            if [ "${LOCAL_MACHINE}" = "thunder" ]; then
              echo "cd ${EXAMPLES_DIR} && srun -n4 ${TEMP_EXAMPLE_FILE} &> ${BASE_DIR}/${LOCAL_MACHINE}-${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.log" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
            else
              echo "cd ${EXAMPLES_DIR} && ./${TEMP_EXAMPLE_FILE} &> ${BASE_DIR}/${LOCAL_MACHINE}-${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.log" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
            fi
            echo "" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
            echo "exit 0" >> "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
          fi
          chmod +x "${TEMP_MODULE}-${TEMP_EXAMPLE_FILE}.job"
        fi
# ----------------------------
      else
        STATUS="FAILED"
      fi
      NUM_EXAMPLES=$((${NUM_EXAMPLES}-1))
    done
  else
    STATUS="SKIPPING - no examples found"
  fi
}

# remove the motd header generated when starting a BASH login shell
remove_junk_header()
{
  while read IN_LINES IN_FILE; do
    if [ "${IN_LINES}" != "" -a "${IN_FILE}" != "" ]; then
      echo "TOTAL_LINES=${IN_LINES}" > junk-temp.txt
      return
    fi
  done
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

MODULE_LIST="cvode cvodes ida idas kinsol"
PROBLEM_FLAG="yes"
BASE_DIR="$1"
KILL_SCRIPT="$2"
REMOTE_USER="$3"
PROJECT_TITLE="$4"
WAIT_TIME="$5"
IN_REMOTE_CMD="$6"
IN_REMOTE_ARGS="$7"
IN_FIX_BASH="$8"

# get local system name
MACHINE_TYPE=`uname -s`
if [ "${MACHINE_TYPE}" = "SunOS" ]; then
  LOCAL_MACHINE=`hostname`
else
  LOCAL_MACHINE=`hostname -s`
fi

# default is to not use PSUB
USE_PSUB="no"

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

# read MPI input file to determine MPI options
MPI_OPT="no"
if [ -r "${LOCAL_MACHINE}".opt ]; then
  MPI_OPT="yes"
  OPT_INPUT_FILE="${LOCAL_MACHINE}.opt"
  OPT_EXEC_NAME=`basename ${OPT_INPUT_FILE}`
  if [ "${OPT_EXEC_NAME}" = "${OPT_INPUT_FILE}" ]; then
    CURRENT_DIR=`pwd`
    OPT_INPUT_FILE="${CURRENT_DIR}/${OPT_INPUT_FILE}"
  fi
  STATUS_OPT=`source "${OPT_INPUT_FILE}" 2>&1 | fgrep "error"`
  if [ "${STATUS_OPT}" = "" ]; then
    source "${OPT_INPUT_FILE}"
  # if file exists but cannot be sourced then abort
  else
    echo "FILE_ERROR: unable to source MPI input file - ${OPT_INPUT_FILE}"
    exit 0
  fi
fi

# determine appropriate execution command for parallel examples given contents of MPI options file for local system
# different systems have different defaults
if [ "${LOCAL_MACHINE}" = "gps" ]; then
  if [ "${MPI_OPT}" = "no" ]; then
    PARALLEL_EXEC_CMD="dmpirun -np 4"
  elif [ "${MPI_VERSION}" = "default" -a "${MPI_COMMAND}" = "default" ]; then
    PARALLEL_EXEC_CMD="dmpirun -np 4"
  elif [ "${MPI_VERSION}" = "default" -a "${MPI_COMMAND}" != "default" ]; then
    PARALLEL_EXEC_CMD="${MPI_COMMAND}"
  else
    if [ "${MPI_DIR}" != "" -a "${MPI_COMMAND}" != "" ]; then
      PARALLEL_EXEC_CMD="${MPI_DIR}/bin/${MPI_COMMAND}"
    elif [ "${MPI_DIR}" = "" -a "${MPI_COMMAND}" != "" ]; then
      PARALLEL_EXEC_CMD="${MPI_COMMAND}"
    fi
  fi

elif [ "${LOCAL_MACHINE}" = "tc2k" ]; then
  if [ "${MPI_OPT}" = "no" ]; then
    PARALLEL_EXEC_CMD="prun -n4 -ppdebug"
  elif [ "${MPI_VERSION}" = "default" -a "${MPI_COMMAND}" = "default" ]; then
    PARALLEL_EXEC_CMD="prun -n4 -ppdebug"
  elif [ "${MPI_VERSION}" = "default" -a "${MPI_COMMAND}" != "default" ]; then
    PARALLEL_EXEC_CMD="${MPI_COMMAND}"
  else
    if [ "${MPI_DIR}" != "" -a "${MPI_COMMAND}" != "" ]; then
      PARALLEL_EXEC_CMD="${MPI_DIR}/bin/${MPI_COMMAND}"
    elif [ "${MPI_DIR}" = "" -a "${MPI_COMMAND}" != "" ]; then
      PARALLEL_EXEC_CMD="${MPI_COMMAND}"
    fi
  fi

# force use of PSUB on blue and frost
elif [ "${LOCAL_MACHINE}" = "blue" -o "${LOCAL_MACHINE}" = "frost" -o "${LOCAL_MACHINE}" = "thunder" ]; then
  PARALLEL_EXEC_CMD=""
  USE_PSUB="yes"

elif [ "${LOCAL_MACHINE}" = "mcr" -o "${LOCAL_MACHINE}" = "ilx" -o "${LOCAL_MACHINE}" = "pengra" -o "${LOCAL_MACHINE}" = "alc" -o "${LOCAL_MACHINE}" = "pvc" ]; then
  if [ "${MPI_OPT}" = "no" ]; then
    PARALLEL_EXEC_CMD="srun -n4"
  elif [ "${MPI_VERSION}" = "default" -a "${MPI_COMMAND}" = "default" ]; then
    PARALLEL_EXEC_CMD="srun -n4"
  elif [ "${MPI_VERSION}" = "default" -a "${MPI_COMMAND}" != "default" ]; then
    PARALLEL_EXEC_CMD="${MPI_COMMAND}"
  else
    if [ "${MPI_DIR}" != "" -a "${MPI_COMMAND}" != "" ]; then
      PARALLEL_EXEC_CMD="${MPI_DIR}/bin/${MPI_COMMAND}"
    elif [ "${MPI_DIR}" = "" -a "${MPI_COMMAND}" != "" ]; then
      PARALLEL_EXEC_CMD="${MPI_COMMAND}"
    fi
  fi

elif [ "${MACHINE_TYPE}" = "Linux" -o "${MACHINE_TYPE}" = "SunOS" ]; then
  if [ "${MPI_OPT}" = "no" ]; then
    PARALLEL_EXEC_CMD="mpirun -np 4"
  elif [ "${MPI_VERSION}" = "default" -a "${MPI_COMMAND}" = "default" ]; then
    PARALLEL_EXEC_CMD="mpirun -np 4"
  elif [ "${MPI_VERSION}" = "default" -a "${MPI_COMMAND}" != "default" ]; then
    PARALLEL_EXEC_CMD="${MPI_COMMAND}"
  else
    if [ "${MPI_DIR}" != "" -a "${MPI_COMMAND}" != "" ]; then
      PARALLEL_EXEC_CMD="${MPI_DIR}/bin/${MPI_COMMAND}"
    elif [ "${MPI_DIR}" = "" -a "${MPI_COMMAND}" != "" ]; then
      PARALLEL_EXEC_CMD="${MPI_COMMAND}"
    fi
  fi
fi

MPI_CMD="${PARALLEL_EXEC_CMD}"

# determine if using LAM-MPI
# if so then can use lamboot and lamclean
if [ "${USE_PSUB}" = "no" ]; then
  PARALLEL_EXEC_CMD=`echo "${MPI_CMD}" | cut -d' ' -f1`
  STATUS_MPIRUN=`basename "${PARALLEL_EXEC_CMD}"`
  if [ "${STATUS_MPIRUN}" != "mpirun" ]; then
    USE_LAM_MPI="no"
  else
    if [ "${MPI_DIR}" != "" ]; then
      MPIRUN_DIR="${MPI_DIR}/bin"
    else
      FIND_MPIRUN=`which "${PARALLEL_EXEC_CMD}"`
      MPIRUN_DIR=`dirname "${FIND_MPIRUN}"`
    fi
    if [ -x "${MPIRUN_DIR}"/lamboot ]; then
      USE_LAM_MPI="yes"
    else
      USE_LAM_MPI="no"
    fi
  fi
elif [ "${USE_PSUB}" = "yes" ]; then
  USE_LAM_MPI="no"
fi

# begin writing program execution script
TEMP_SCRIPT_NAME=`echo "${LOCAL_MACHINE}" | cut -c1-5`
SCRIPT_NAME="${TEMP_SCRIPT_NAME}.sh"
touch "${SCRIPT_NAME}"
# must use a login shell
if [ "${MACHINE_TYPE}" = "Linux" ]; then
  echo "#!/bin/sh --login" > "${SCRIPT_NAME}"
else
  echo "#!/bin/bash --login" > "${SCRIPT_NAME}"
fi
echo "" >> "${SCRIPT_NAME}"
echo "touch script_running.info" >> "${SCRIPT_NAME}"
echo "" >> "${SCRIPT_NAME}"
echo "echo \"START_HERE\"" >> "${SCRIPT_NAME}"

# if necessary, update executable search path
# also, add LAM-specific commands to script if using LAM-MPI
if [ "${MPI_OPT}" = "yes" -a "${USE_PSUB}" = "no" ]; then
  if [ "${MPI_VERSION}" != "default" -a "${MPI_DIR}" != "" ]; then
    echo "export PATH=${MPI_DIR}:\${PATH}" >> "${SCRIPT_NAME}"
  fi
  if [ "${USE_LAM_MPI}" = "yes" ]; then
    echo "lamboot 1> /dev/null 2> /dev/null" >> "${SCRIPT_NAME}"
    echo "" >> "${SCRIPT_NAME}"
  fi
fi

chmod +x "${SCRIPT_NAME}"

# if not using PSUB then generate information file for zombie script
if [ "${USE_PSUB}" = "no" ]; then
  KILL_INFO=`basename "${KILL_SCRIPT}" .sh`
  KILL_INFO="${KILL_INFO}.txt"
  touch "${KILL_INFO}"
fi

# make certain script has correct path to BASH shell interpreter
if [ "${MACHINE_TYPE}" = "Linux" ]; then
  :
else
  ./${IN_FIX_BASH} "${SCRIPT_NAME}" yes
fi

# search all relevant directories for examples and generate script
NUM_MODULES=`echo "${MODULE_LIST}" | wc -w`
while [ $((${NUM_MODULES})) -gt 0 ]; do
  TEMP_MODULE=`echo "${MODULE_LIST}" | cut -d' ' -f$((${NUM_MODULES}))`
  echo -ne "\nChecking module ${TEMP_MODULE}..."
  if [ -d "${BASE_DIR}/examples/${TEMP_MODULE}" ]; then
    MODULE_DIR="${BASE_DIR}/examples/${TEMP_MODULE}"
    echo "[CHECKING]"

    if [ -d "${MODULE_DIR}/serial" ]; then
      echo -ne "\tChecking for serial examples..."
      run_examples ${MODULE_DIR}/serial
      echo "[${STATUS}]"
    fi

    if [ -d "${BASE_DIR}/test_examples/${TEMP_MODULE}/serial" ]; then
      MODULE_TEST_DIR="${BASE_DIR}/test_examples/${TEMP_MODULE}"
      echo -ne "\tChecking for serial dev examples..."
      run_examples ${MODULE_TEST_DIR}/serial
      echo "[${STATUS}]"
    fi

    if [ -d "${MODULE_DIR}/parallel" ]; then
      echo -ne "\tChecking for parallel examples..."
      run_examples ${MODULE_DIR}/parallel
      echo "[${STATUS}]"
    fi

    if [ -d "${BASE_DIR}/test_examples/${TEMP_MODULE}/parallel" ]; then
      MODULE_TEST_DIR="${BASE_DIR}/test_examples/${TEMP_MODULE}"
      echo -ne "\tChecking for parallel dev examples..."
      run_examples ${MODULE_TEST_DIR}/parallel
      echo "[${STATUS}]"
    fi

    if [ -d "${MODULE_DIR}/fcmix_serial" ]; then
      echo -ne "\tChecking serial fcmix examples..."
      run_examples ${MODULE_DIR}/fcmix_serial
      echo "[${STATUS}]"
    fi

    if [ -d "${MODULE_DIR}/fcmix_parallel" ]; then
      echo -ne "\tChecking parallel fcmix examples..."
      run_examples ${MODULE_DIR}/fcmix_parallel
      echo "[${STATUS}]"
    fi

  else
    echo "[SKIPPING]"
  fi
  NUM_MODULES=$((${NUM_MODULES}-1))
done

# add wait loop to script if using PSUB
# script will not return until all queued examples have finished
if [ "${USE_PSUB}" = "yes" ]; then
  echo "" >> "${SCRIPT_NAME}"
  echo "JOB_USER=\`whoami\`" >> "${SCRIPT_NAME}"
  echo "STILL_RUNNING=\`pstat -u \${JOB_USER} -m ${LOCAL_MACHINE}\`" >> "${SCRIPT_NAME}"
  echo "while [ \"\${STILL_RUNNING}\" != \"\" ]; do" >> "${SCRIPT_NAME}"
  echo "  sleep 300" >> "${SCRIPT_NAME}"
  echo "  STILL_RUNNING=\`pstat -u \${JOB_USER} -m ${LOCAL_MACHINE}\`" >> "${SCRIPT_NAME}"
  echo "done" >> "${SCRIPT_NAME}"
# if not using PSUB but using LAM-MPI then add LAM-specific commands at end of script to perform clean-up
elif [ "${USE_PSUB}" = "no" ]; then
  if [ "${MPI_VERSION}" = "lam" -o "${USE_LAM_MPI}" = "yes" ]; then
    echo "" >> "${SCRIPT_NAME}"
    echo "sleep 10" >> "${SCRIPT_NAME}"
    echo "sync" >> "${SCRIPT_NAME}"
    echo "lamclean 1> /dev/null 2> /dev/null" >> "${SCRIPT_NAME}"
    echo "lamhalt 1> /dev/null 2> /dev/null" >> "${SCRIPT_NAME}"
  fi
fi
echo "" >> "${SCRIPT_NAME}"
if [ "${MACHINE_TYPE}" != "OSF1" ]; then
  echo "sync" >> "${SCRIPT_NAME}"
fi
# add commands to script to put all log files into a compressed tar ball for easy retrieval by main system
echo "cd ${BASE_DIR}" >> "${SCRIPT_NAME}"
echo "mv -f config.log config-log" >> "${SCRIPT_NAME}"
echo "tar -cf ${LOCAL_MACHINE}-examples-logs.tar *.log" >> "${SCRIPT_NAME}"
echo "gzip ${LOCAL_MACHINE}-examples-logs.tar" >> "${SCRIPT_NAME}"
echo "" >> "${SCRIPT_NAME}"
echo "exit 0" >> "${SCRIPT_NAME}"

echo -e "\nExecuting automatically generated script..."
if [ "${MACHINE_TYPE}" != "OSF1" ]; then
  sync
fi
# only execute zombie script if not using PSUB
if [ "${USE_PSUB}" = "no" ]; then
  ./${KILL_SCRIPT} "${KILL_INFO}" "${WAIT_TIME}" "${SCRIPT_NAME}" "${USE_LAM_MPI}" "${REMOTE_USER}" "${PROJECT_TITLE}" "${IN_REMOTE_CMD}" "${IN_REMOTE_ARGS}" "${MPI_DIR}" 1> kill.output 2> kill.output &
  sleep 10
fi
# execute script and direct output (stdio and stderr) to specified file
./${SCRIPT_NAME} &> temp-output.txt
# remove motd header from log file
SKIP_LINES=`fgrep -n "START_HERE" temp-output.txt | cut -d':' -f1`
TOTAL_LINES=`wc -l temp-output.txt`
echo "${TOTAL_LINES}" | remove_junk_header
source junk-temp.txt
SKIP_LINES=$((${TOTAL_LINES} - ${SKIP_LINES}))
rm -f junk-temp.txt
if [ "${MACHINE_TYPE}" = "SunOS" ]; then
  tail -$((${SKIP_LINES})) temp-output.txt > junk-temp.txt
else
  tail -n $((${SKIP_LINES})) temp-output.txt > junk-temp.txt
fi
cp -f junk-temp.txt temp-output.txt
rm -f junk-temp.txt

# make script output available to main system via standard otuput (stdio)
cat temp-output.txt
if [ "${USE_PSUB}" = "no" ]; then
  sleep 120
  cat kill.output
  rm -f kill.output
fi
rm -f temp-output.txt
rm -f script_running.info
echo -e "[DONE]\n"

exit 0
