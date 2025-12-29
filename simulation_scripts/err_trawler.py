import os
import sys
import glob
import numpy as np

#errfiles=sorted(glob.glob('/scratch/zrechav/simulation_scripts/21537_step_1_snowstorm_propagation/logs/*.err'))
errfiles = sorted(glob.glob('/scratch/vbasu/simulation_scripts/2201*_step_1_snowstorm_propagation/logs/*.err'))
print(len(errfiles))
badfile = open('/home/zrechav/Vedant_bad_files.txt', 'a')
for file in errfiles:
    with open(file,'r') as f:
        errname = file.split('/')[-1]
        for line in f.readlines():
            if 'Maximum number of photons exceeded' in line:
                #print(line)
                badfile.write(errname + '\n')
    f.close()      
badfile.close()