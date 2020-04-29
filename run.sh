#!/bin/bash

BASENAME="${0##*/}"
log () {
  echo "${BASENAME}: ${1}"
}

# Set child by default switch to main if on main node container
NODE_TYPE="child"
if [ "${AWS_BATCH_JOB_MAIN_NODE_INDEX}" == "${AWS_BATCH_JOB_NODE_INDEX}" ]; then
  NODE_TYPE="main"
fi

# If host is child node, then exit
if [ "${NODE_TYPE}" = "child" ]; then
    log "Immediately terminated since this is a child node"
    exit 0
fi

sleep 2 # Avoid completing main node before exiting child node

# Detect the number of local and physical cores
NUM_LOGICAL_CORES=$(lscpu -p | egrep -v '^#' | wc -l)
NUM_PHYSICAL_CORES=$(lscpu -p | egrep -v '^#' | sort -u -t, -k 2,4 | wc -l)
log "Host has ${NUM_LOGICAL_CORES} logical and ${NUM_PHYSICAL_CORES} physical cores"

# Computing Environment of Parallel Track in SAT 2020 Competition:
#  m4.16xlarge, which has 64 virtual cores and 256GB of memory (32 physical cores)
# In this branch, we use 64 threads
NUM_THREADS=64

# If NUM_THREADS is not defined, we use NUM_PHYSICAL_CORES as NUM_THREADS
if [ -z "${NUM_THREADS}" ]; then
    NUM_THREADS=${NUM_PHYSICAL_CORES}
fi
log "Solver will invoke ${NUM_THREADS} threads"

log "Downloading problem from s3://${S3_BKT}/${COMP_S3_PROBLEM_PATH}"
if [[ "${COMP_S3_PROBLEM_PATH}" == *".xz" ]];
then
  aws s3 cp s3://${S3_BKT}/${COMP_S3_PROBLEM_PATH} supervised-scripts/test.cnf.xz
  unxz supervised-scripts/test.cnf.xz
else
  aws s3 cp s3://${S3_BKT}/${COMP_S3_PROBLEM_PATH} supervised-scripts/test.cnf
fi

## Main
COMMAND="/manyglucose-4.1-60/parallel/manyglucose-4.1-60 -verb=0 -model -real-time-lim=5000 -nthreads=${NUM_THREADS} supervised-scripts/test.cnf"
log "Invoking solver: ${COMMAND}"
(time ${COMMAND}) &> test.log
RET_CODE=$?

if [ -z "${COMP_S3_RESULT_PATH}" ]; then
    cat test.log
else
    cat test.log
    aws s3 cp test.log s3://${S3_BKT}/${COMP_S3_RESULT_PATH}
fi

case "${RET_CODE}" in
    "10"         ) EXIT_CODE=0 ;; # SATISFIABLE
    "20" | "138" ) EXIT_CODE=0 ;; # UNSATISFIABLE
    *            ) EXIT_CODE=${RET_CODE} ;;
esac

exit ${EXIT_CODE}
