#!/bin/bash

#prep inputs
input=$1
output=$2

#setup env
#source /home/zrechav/testdev_sh.sh
eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/setup.sh`
export HDF5_USE_FILE_LOCKING='FALSE'
#run job                                              

/cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/RHEL_7_x86_64/metaprojects/combo/V01-01-01/env-shell.sh python /home/zrechav/SelfVeto_Correlation_Tables/scripts/compiling/OutputFileCompiler.py -i ${input} -o ${output} 
