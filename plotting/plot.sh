#!/bin/bash

# plotting command
BINARY=Rscript

# sanity check, can only accept True or False
if [ -z "$1" ] || [ -z "$2" ];  then
  echo "Usage: $(basename $0) [LOG_SCALE] (True or False) [END_TIME] (in seconds)"
  exit 1
fi

LOG_SCALE=${1:0:1}
END_TIME=$2

if [ "$LOG_SCALE" = "T" -o "$LOG_SCALE" = "F" ]; then
  $BINARY plotBranch.R "$LOG_SCALE" "$END_TIME" # before branch "right"
  $BINARY plotAfterBranch.R "$LOG_SCALE" "$END_TIME" # after branch
  $BINARY plotAfterBranch2.R "$LOG_SCALE" "$END_TIME" # after branch2
  $BINARY plotAfterBranch3.R "$LOG_SCALE" "$END_TIME" # after branch3
  $BINARY plotDest.R "$LOG_SCALE" "$END_TIME" # soma
  $BINARY plotSource.R "$LOG_SCALE" "$END_TIME" # neurite tip
fi
