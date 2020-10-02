#!/bin/bash

for file in syn_simu*.sh; do
  echo "Submtting >> ${file} << now ..."
  qsub "$file"
done
