#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/icetray-start
#METAPROJECT /home/zrechav/i3/icetray/build
import numpy as np
import math
import os.path
from os import path
from icecube.tableio import I3TableWriter
#from icecube.hdfwriter import I3HDFTableService, I3HDFWriter
import tables
import h5py
import numpy as np
import math
import os.path
from os import path
from icecube.tableio import I3TableWriter
#from icecube.hdfwriter import I3HDFTableService, I3HDFWriter
import tables
import glob
import yaml
import h5py
import pandas as pd
from icecube.dataclasses import I3Double, I3Particle, I3Direction, I3Position, I3VectorI3Particle, I3Constants, I3VectorOMKey
import traceback
import simweights




with open('/home/zrechav/SelfVeto_Correlation_Tables/scripts/config.yaml', 'r') as yaml_file:
    config = yaml.safe_load(yaml_file)


Muon_Energy_L1, Muon_Energy_L2, Muon_Energy_L3, Muon_Energy_L4, Muon_Energy_L5 = np.array([]), np.array([]), np.array([]), np.array([]), np.array([])

Muon_Depth_L1, Muon_Depth_L2, Muon_Depth_L3, Muon_Depth_L4, Muon_Depth_L5 = np.array([]), np.array([]), np.array([]), np.array([]), np.array([])

Muon_Zenith_L1, Muon_Zenith_L2, Muon_Zenith_L3, Muon_Zenith_L4, Muon_Zenith_L5, = np.array([]), np.array([]), np.array([]), np.array([]), np.array([])

Muon_Neutrino_separation_L1, Muon_Neutrino_separation_L2, Muon_Neutrino_separation_L3, Muon_Neutrino_separation_L4, Muon_Neutrino_separation_L5 = np.array([]), np.array([]), np.array([]), np.array([]), np.array([])

Muon_Shower_separation_L1, Muon_Shower_separation_L2, Muon_Shower_separation_L3, Muon_Shower_separation_L4, Muon_Shower_separation_L5 = np.array([]), np.array([]), np.array([]), np.array([]), np.array([])

MuonMultiplicity= np.array([])



Max_Deposited_Energy = np.array([])     
Energy_Shower_Neutrino = np.array([])   
Zenith_Shower_Neutrino = np.array([])  
Depth_Shower_Neutrino = np.array([])    


Total_Muon_energy = np.array([])

Flavour_Shower_Neutrino = np.array([])
GaisserH4a_weight = np.array([])
PrimaryMass, PrimaryEnergyperNucleon= np.array([]),np.array([])

default_filename = config['corsika_sample']
nfiles = eval(config['nfiles'])

import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", type=str, default = default_filename, dest = "INPUT", help = "input i3 file directory")

args = parser.parse_args()

print(args)

#filename = path(args.INPUT)
filename = config['corsika_sample']

#filename = '/data/user/zrechav/compiled_hdf5s/airshowered_corsika/correlation/22803_0000000-0009999_correlated.hdf5' 
print(filename)


def calc_weight(inp = filename,numfiles = nfiles):
    with pd.HDFStore(inp, "r") as hdffile:
        weighter = simweights.CorsikaWeighter(hdffile, nfiles=numfiles)
        flux = simweights.GaisserH4a()
        weights = weighter.get_weights(flux)  
    return weights
    

if (path.exists(filename)):
    try:
        print('I EXIST')
        #GaisserH4a_weight = calc_weight()
        #print(GaisserH4a_weight)
        
        
        hdf=h5py.File(filename, 'r')
        
        ##weighting info?
        
        print('lets begin')
########ENERGY, ZENITH, AND WEIGHTING INFORMATION#######
        #GaisserH4a_weight = np.append(GaisserH4a_weight,np.asarray(hdf.get('weights')['value']))
        Zenith_Shower_Neutrino = np.append(Zenith_Shower_Neutrino,np.asarray((hdf.get('shower_neutrino_zenith')['value'])))
        print('skittles')
        Flavour_Shower_Neutrino = np.append(Flavour_Shower_Neutrino,np.asarray(hdf.get('shower_neutrino_type')['value']))
        print('cheez its')
        Energy_Shower_Neutrino = np.append(Energy_Shower_Neutrino,np.asarray(hdf.get('shower_neutrino_energy')['value']))
        print('i live')
        Depth_Shower_Neutrino = np.append(Depth_Shower_Neutrino,np.asarray(hdf.get('shower_neutrino_depth')['value']))
        print('oh ima goin')
        

#########MUON INFORMATION############
        MuonMultiplicity = np.append(MuonMultiplicity,np.asarray(hdf.get('MuonMultiplicity')['value']))
        print('making progress')
        Muon_Energy_L1 = np.append(Muon_Energy_L1, np.asarray(hdf.get('Muon_Energy_L1')['value']))
        print('making progress')
        Muon_Energy_L2 = np.append(Muon_Energy_L2, np.asarray(hdf.get('Muon_Energy_L2')['value']))
        print('making progress')
        Muon_Energy_L3 = np.append(Muon_Energy_L3, np.asarray(hdf.get('Muon_Energy_L3')['value']))
        print('L3 muon energy : ',len(Muon_Energy_L3))
        print('making progress')
        Muon_Energy_L4 = np.append(Muon_Energy_L4, np.asarray(hdf.get('Muon_Energy_L4')['value']))
        #Muon_Energy_L5 = np.append(Muon_Energy_L5, np.asarray(hdf.get('Muon_Energy_L5')['value']))
        print('making progress')
        Muon_Depth_L1 = np.append(Muon_Depth_L1, np.asarray(hdf.get('Muon_L1_Depth')['value']))
        print('UAU')
        Muon_Depth_L2 = np.append(Muon_Depth_L2, np.asarray(hdf.get('Muon_L2_Depth')['value']))
        print('UAU')
        Muon_Depth_L3 = np.append(Muon_Depth_L3, np.asarray(hdf.get('Muon_L3_Depth')['value']))
        print('UAU')
        Muon_Depth_L4 = np.append(Muon_Depth_L4, np.asarray(hdf.get('Muon_L4_Depth')['value']))
        #Muon_Depth_L5 = np.append(Muon_Depth_L5, np.asarray(hdf.get('Muon_L5_Depth')['value']))
        print('UAU')
        Muon_Zenith_L1=np.append(Muon_Zenith_L1,np.asarray(hdf.get('Muon_L1')['zenith']))
        print('yeet')
        Muon_Zenith_L2=np.append(Muon_Zenith_L2,np.asarray(hdf.get('Muon_L2')['zenith']))
        print('yeet')
        Muon_Zenith_L3=np.append(Muon_Zenith_L3,np.asarray(hdf.get('Muon_L3')['zenith']))
        print('yeet')
        Muon_Zenith_L4=np.append(Muon_Zenith_L4,np.asarray(hdf.get('Muon_L4')['zenith']))

        print('almost therrrrr')

        Total_Muon_energy = np.append(Total_Muon_energy,np.asarray(hdf.get('Total_Muon_Energy')['value']))
        print('I FINISHED LOADING MY ARRAYS WITH HDF INFO')
        hdf.close()
    except Exception as e:
        print(filename+" Is Faulty")
        print(e)
        print(f"Occurred on line: {traceback.extract_tb(e.__traceback__)[-1][1]}")

        hdf.close()
        
MuonEnergyThreshMask=[energy>10. for energy in Total_Muon_energy]
print('MuonEnergyThreshMask : ', len(MuonEnergyThreshMask))
##MASKING AND SAVING INFORMATION TO SAVE TIME LATER ON                   
np.save(config['Zenith_Shower_Neutrino'], Zenith_Shower_Neutrino[MuonEnergyThreshMask])
np.save(config['Flavour_Shower_Neutrino'], Flavour_Shower_Neutrino[MuonEnergyThreshMask])
np.save(config['Energy_Shower_Neutrino'], Energy_Shower_Neutrino[MuonEnergyThreshMask])
np.save(config['Depth_Shower_Neutrino'], Depth_Shower_Neutrino[MuonEnergyThreshMask])
np.save(config['MuonMultiplicity'], MuonMultiplicity[MuonEnergyThreshMask])
np.save(config['Muon_Energy_L1'], Muon_Energy_L1[MuonEnergyThreshMask])
np.save(config['Muon_Energy_L2'], Muon_Energy_L2[MuonEnergyThreshMask])
np.save(config['Muon_Energy_L3'], Muon_Energy_L3[MuonEnergyThreshMask])
np.save(config['Muon_Energy_L4'], Muon_Energy_L4[MuonEnergyThreshMask])
np.save(config['Muon_Zenith_L1'], Muon_Zenith_L1[MuonEnergyThreshMask])
np.save(config['Muon_Zenith_L2'], Muon_Zenith_L2[MuonEnergyThreshMask])
np.save(config['Muon_Zenith_L3'], Muon_Zenith_L3[MuonEnergyThreshMask])
np.save(config['Muon_Zenith_L4'], Muon_Zenith_L4[MuonEnergyThreshMask])
np.save(config['Muon_Depth_L1'], Muon_Depth_L1[MuonEnergyThreshMask])
np.save(config['Muon_Depth_L2'], Muon_Depth_L2[MuonEnergyThreshMask])
np.save(config['Muon_Depth_L3'], Muon_Depth_L3[MuonEnergyThreshMask])
np.save(config['Muon_Depth_L4'], Muon_Depth_L4[MuonEnergyThreshMask])  
np.save(config['Total_Muon_energy'], Total_Muon_energy[MuonEnergyThreshMask])

Weights = calc_weight()
np.save(config['Weights'], Weights[MuonEnergyThreshMask])
