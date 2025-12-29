import os
import sys
import glob
import numpy as np

dataset = '22803'
subdir = '0000000-0019999'

def get_Grid_job(infile):
    
    job_name = 'make_pickle '
    lines =[
        'JOB ' + job_name + ' /home/zrechav/SelfVeto_Correlation_Tables/scripts/picklemaking/picklemaking.sub',
        'VARS ' + job_name + ' infile="' + str(infile) + '" ',
        'Retry ' + job_name + ' 2',
        ]
    return lines

input_path = '/data/ana/Diffuse/DNNCascades_Diffuse/Corsika_Temp/combined/'+ dataset + '_' + subdir + '_correlated.hdf5'
#output_path = '/data/user/zrechav/compiled_hdf5s/airshowered_corsika/' + dataset + '/correlation/'+ dataset + '_' + subdir + '_correlated.hdf5'

input_files = input_path


dag_lines = []
dag_lines.extend(
get_Grid_job(
    infile=input_files,

    )
)
outfile_name='/home/zrechav/SelfVeto_Correlation_Tables/scripts/picklemaking/dags/' + dataset + '_' + subdir + '_picklemaker.dag'
#outfile_name='/home/zrechav/SelfVeto_Correlation_Tables/scripts/compiling/dags/' + dataset + '_' + subdir + '_weight_airshowered_corsika.dag'
with open(outfile_name, 'w') as f:
    f.write('\n'.join(dag_lines))

