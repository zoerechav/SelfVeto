#!/bin/bash

#prep inputs
infile=$1
#setup env
eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/setup.sh`
source /home/zrechav/.venvs/icetray_version_lock_venv/bin/activate

export HDF5_USE_FILE_LOCKING='FALSE'
echo ${infile}
echo ${outfile}
/home/zrechav/i3/icetray/build/env-shell.sh python /home/zrechav/SelfVeto_Correlation_Tables/scripts/picklemaking/picklemaker.py -i ${infile}   
