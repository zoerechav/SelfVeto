print('hey bro')
import os
import sys
import glob
import numpy as np

##universal variables

dataset = '22803'
#subdir = '0000000-0000999'
#short_subdir = '00000-00999'
nfiles=4
Compiled_path='/data/user/zrechav/compiled_hdf5s/airshowered_corsika/correlation/'

input_base_path='/data/user/zrechav/compiled_hdf5s/airshowered_corsika/'  + dataset +'/' 

def get_Grid_job(infile,outfile):  
    job_name = 'COMILE_' 
    lines =[
        'JOB ' + job_name+"_Compile" + ' /home/zrechav/SelfVeto_Correlation_Tables/scripts/compiling/OutputFileCompiler_NPX.sub',
        'VARS ' + job_name+"_Compile" + ' infile="' + infile + '" outfile="' + outfile ,
        'Retry ' + job_name+"_Compile" + ' 2']
    return lines            

def get_Grid_jobs(Compiled_path=None,input_path=None):    
    raw_infile_prefix=input_path+'*_weighted.hdf5' 
    lines = []  
    infiles=raw_infile_prefix
    filenum=str(nfiles)   
    CompiledFile=Compiled_path#+dataset 
    print(CompiledFile)
    lines.extend(get_Grid_job(infiles,CompiledFile))
    return lines

def dag_lines(input_path,compiled_path):
    dag_lines = []
    dag_lines.extend(
    get_Grid_jobs(
        input_path=input_base_path,
        Compiled_path=compiled_path
        )
    )
    outfile_name='/home/zrechav/SelfVeto_Correlation_Tables/scripts/compiling/dags/GRAND_COMPILE_' +dataset + '.dag' #outputcompile_vertex_flag_NuGen_bfr1_' + dataset + '.dag'
    #outfile_name='/home/zrechav/test/output_compiler_MuonGun_'+name+'.dag'
    print(outfile_name)
    print(dag_lines)
    with open(outfile_name, 'w') as f:
        f.write('\n'.join(dag_lines))
    dag_lines = []


dag_lines(input_base_path, Compiled_path)
    
print('peace out bro')