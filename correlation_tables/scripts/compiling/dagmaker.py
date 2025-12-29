import os
import sys
import glob
import numpy as np

dataset = '22803'
subdir = '0019000-0019999'

def get_Grid_job(infile,outfile):
    
    job_name = 'make_sampled_corsika_hdf'
    lines =[
        'JOB ' + job_name + ' /home/zrechav/SelfVeto_Correlation_Tables/scripts/compiling/make_sampled_corsika_hdf.sub',
        'VARS ' + job_name + ' infile="' + infile + '" outfile="' + outfile +'" ',
        'Retry ' + job_name + ' 2',
        ]
    return lines


input_path = '/data/user/zrechav/airshowered_corsika/' + dataset + '/random_test/'
#output_path = '/data/ana/Diffuse/DNNCascades_Diffuse/Corsika_Temp/' + dataset + '_' + subdir + '_unweighted.hdf5'
output_path = '/data/ana/Diffuse/DNNCascades_Diffuse/Corsika_Temp/random_test' + dataset + '_' + subdir + '_unweighted.hdf5'

# def get_Grid_job(infile,outfile):
    
#     job_name = 'weight_airshowered_corsika'
#     lines =[
#         'JOB ' + job_name + ' /home/zrechav/SelfVeto_Correlation_Tables/scripts/compiling/weight_airshowered_corsika.sub',
#         'VARS ' + job_name + ' infile="' + infile + '" outfile="' + outfile +'" ',
#         'Retry ' + job_name + ' 2',
#         ]
#     return lines



# input_path = '/data/user/zrechav/compiled_hdf5s/airshowered_corsika/' + dataset + '/' + dataset + '_' + subdir + '_unweighted.hdf5'
# output_path = '/data/user/zrechav/compiled_hdf5s/airshowered_corsika/' + dataset + '/' + dataset + '_' + subdir + '_weighted.hdf5'

dag_lines = []
dag_lines.extend(
get_Grid_job(
    infile=input_path,
    outfile = output_path,

    )
)
outfile_name='/home/zrechav/SelfVeto_Correlation_Tables/scripts/compiling/dags/' + dataset + '_' + subdir + '_make_sampled_corsika_hdf_random_test.dag'
# outfile_name='/home/zrechav/SelfVeto_Correlation_Tables/scripts/compiling/dags/' + dataset + '_' + subdir + '_weight_airshowered_corsika.dag'
with open(outfile_name, 'w') as f:
    f.write('\n'.join(dag_lines))

