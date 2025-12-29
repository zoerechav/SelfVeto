import os
import sys
import glob
import numpy as np
#/data/sim/IceCube/2020/filtered/level2/CORSIKA-in-ice/20904/0000000-0000999/
dataset = '22803'
subdir = '0000000-0000999'
#short_subdir = '00000-00999'
#tag = 'not_contained'
def get_Grid_job(infile,filenum,outfile,gcd, nugen,muongun,corsika,burn):
    if nugen==1:
        job_name = 'propogate_' +filenum
        lines =[
            'JOB ' + job_name + ' /home/zrechav/SelfVeto_Correlation_Tables/scripts/air_shower_reader/AirShowerReader_NPX.sub',
            'VARS ' + job_name + ' infile="' + infile + '" outfile="' + outfile +'" filenum="'+str(filenum)+'"',
            'Retry ' + job_name + ' 2',
            ]
        return lines


    
def get_Grid_jobs(input_path=None,output_path=None, gcd=None,nugen=1,muongun=0,corsika=0,burn=0,propagate=0):
    raw_infiles=sorted(glob.glob(input_path+ '/*.i3.zst'))[:nfiles]
    infile_numbers = [infile.split('/')[-1].split('.')[2] for infile in raw_infiles]
    print(infile_numbers[0])
    lines = []
    for (infile, filenum) in zip(raw_infiles,infile_numbers):
        
        #outfile=output_path+"vertex_split_" + dataset + '_' + filenum
        outfile = output_path + 'airshowered_corsika_' + dataset + '_' + filenum
        
        #print(infile)
        lines.extend(get_Grid_job(infile,filenum,outfile,gcd,nugen,muongun,corsika,burn))
    return lines


nugen=1
muongun=0

corsika=0
burn =0
nfiles=10
dag_lines = []
if nugen ==1:
    gcd = '/cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_2020.Run134142.Pass2_V0.i3.gz'
    input_path = '/data/sim/IceCube/2023/generated/CORSIKA-in-ice/'+dataset +'/' + subdir + '/generated/'
    #input_path = '/data/user/zrechav/propagated_corsika/20904/' + subdir + '/'
    #input_path = '/data/user/zrechav/propagated_corsika/22784/0000000-0000999/'
    #input_path = '/data/user/zrechav/propagated_corsika/' + dataset + '/' + subdir + '/'
    #input_path = '/data/user/zrechav/test_corsika/'
    #input_path = '/data/sim/IceCube/2023/generated/CORSIKA-in-ice/' + dataset + '/' + subdir + '/generated/'
    output_path = '/data/user/zrechav/airshowered_corsika/22803/random_test/'
    dag_lines.extend(
    get_Grid_jobs(
        input_path=input_path,
        output_path = output_path,
        gcd=gcd,
        nugen=1,
        muongun=0,corsika=0,propagate=0,burn=0
        )
    )
outfile_name='/home/zrechav/SelfVeto_Correlation_Tables/scripts/air_shower_reader/dags/random_test_' + dataset + '_' + subdir + '.dag'
with open(outfile_name, 'w') as f:
    f.write('\n'.join(dag_lines))

