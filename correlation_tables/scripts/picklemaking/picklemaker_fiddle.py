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

filename = config['corsika_sample']

if (path.exists(filename)):
    try:
        print('I EXIST')
        with h5py.File(filename, 'r') as hdf:
    
            #test_primary_energy = hdf_file['PolyplopiaPrimary']['energy'][:]  # Load data
            #test_weights = hdf_file['weights'][:]  # Load weights data
        
        #hdf=h5py.File(filename, 'r')
########ENERGY, ZENITH, AND WEIGHTING INFORMATION#######
            GaisserH4a_weight = np.append(GaisserH4a_weight,
                                          np.asarray(hdf['weights'][:]))

            Zenith_Shower_Neutrino = np.append(Zenith_Shower_Neutrino,
                                               np.asarray((hdf['shower_neutrino_zenith'][:])))

            Flavour_Shower_Neutrino = np.append(Flavour_Shower_Neutrino, 
                                                np.asarray(hdf['shower_neutrino_type'][:]))
            Energy_Shower_Neutrino = np.append(Energy_Shower_Neutrino, 
                                               np.asarray(hdf['shower_neutrino_energy'][:]))
            Depth_Shower_Neutrino = np.append(Depth_Shower_Neutrino,
                                              np.asarray(hdf['shower_neutrino_depth'][:]))
            print('oh ima goin')
        

#########MUON INFORMATION############
        MuonMultiplicity = np.append(MuonMultiplicity,np.asarray(hdf.get('MuonMultiplicity')['value']))

        Muon_Energy_L1 = np.append(Muon_Energy_L1, np.asarray(hdf.get('Muon_Energy_L1')['value']))
        Muon_Energy_L2 = np.append(Muon_Energy_L2, np.asarray(hdf.get('Muon_Energy_L2')['value']))
        Muon_Energy_L3 = np.append(Muon_Energy_L3, np.asarray(hdf.get('Muon_Energy_L3')['value']))
        Muon_Energy_L4 = np.append(Muon_Energy_L4, np.asarray(hdf.get('Muon_Energy_L4')['value']))
        Muon_Energy_L5 = np.append(Muon_Energy_L5, np.asarray(hdf.get('Muon_Energy_L5')['value']))
        print('making progress')
        Muon_Depth_L1 = np.append(Muon_Depth_L1, np.asarray(hdf.get('Muon_L1_Depth')['value']))
        Muon_Depth_L2 = np.append(Muon_Depth_L2, np.asarray(hdf.get('Muon_L2_Depth')['value']))
        Muon_Depth_L3 = np.append(Muon_Depth_L3, np.asarray(hdf.get('Muon_L3_Depth')['value']))
        Muon_Depth_L4 = np.append(Muon_Depth_L4, np.asarray(hdf.get('Muon_L4_Depth')['value']))
        Muon_Depth_L5 = np.append(Muon_Depth_L5, np.asarray(hdf.get('Muon_L5_Depth')['value']))

        Muon_Zenith_L1=np.append(Muon_Zenith_L1,np.asarray(hdf.get('Muon_L1')['zenith']))
        Muon_Zenith_L2=np.append(Muon_Zenith_L2,np.asarray(hdf.get('Muon_L2')['zenith']))
        Muon_Zenith_L3=np.append(Muon_Zenith_L3,np.asarray(hdf.get('Muon_L3')['zenith']))
        Muon_Zenith_L4=np.append(Muon_Zenith_L4,np.asarray(hdf.get('Muon_L4')['zenith']))
        Muon_Zenith_L5=np.append(Muon_Zenith_L5,np.asarray(hdf.get('Muon_L5')['zenith']))
        print('UAU')
        Muon_Neutrino_separation_L1 = np.append(Muon_Neutrino_separation_L1, np.asarray(hdf.get('Muon_L1_Neutrino_separation')['value']))
        Muon_Neutrino_separation_L2 = np.append(Muon_Neutrino_separation_L2, np.asarray(hdf.get('Muon_L2_Neutrino_separation')['value']))
        Muon_Neutrino_separation_L3 = np.append(Muon_Neutrino_separation_L3, np.asarray(hdf.get('Muon_L3_Neutrino_separation')['value']))
        Muon_Neutrino_separation_L4 = np.append(Muon_Neutrino_separation_L4, np.asarray(hdf.get('Muon_L4_Neutrino_separation')['value']))
        Muon_Neutrino_separation_L5 = np.append(Muon_Neutrino_separation_L5, np.asarray(hdf.get('Muon_L5_Neutrino_separation')['value']))
        print('almost therrrrr')
        Muon_Shower_separation_L1 = np.append(Muon_Shower_separation_L1, np.asarray(hdf.get('Muon_L1_Shower_separation')['value']))
        Muon_Shower_separation_L2 = np.append(Muon_Shower_separation_L2, np.asarray(hdf.get('Muon_L2_Shower_separation')['value']))
        Muon_Shower_separation_L3 = np.append(Muon_Shower_separation_L3, np.asarray(hdf.get('Muon_L3_Shower_separation')['value']))
        Muon_Shower_separation_L4 = np.append(Muon_Shower_separation_L4, np.asarray(hdf.get('Muon_L4_Shower_separation')['value']))
        Muon_Shower_separation_L5 = np.append(Muon_Shower_separation_L5, np.asarray(hdf.get('Muon_L5_Shower_separation')['value']))

        Total_Muon_energy = np.append(Total_Muon_energy,np.asarray(hdf.get('Total_Muon_energy')['value']))
        print('I FINISHED LOADING MY ARRAYS WITH HDF INFO')
        hdf.close()
    except Exception as e:
        print(filename+" Is Faulty")
        print(e)
        print(f"Occurred on line: {traceback.extract_tb(e.__traceback__)[-1][1]}")

        hdf.close()

##MASKING AND SAVING INFORMATION TO SAVE TIME LATER ON                
            
MuonEnergyThreshMask=[energy>10. for energy in Muon_Energy_L1]


# print('Flavour masks')
Flav_mask_NuE=np.where(np.abs(Flavour_Shower_Neutrino[MuonEnergyThreshMask])==12., True, False)
np.save(config['flavour_mask_NuE_path'],Flav_mask_NuE)
Flav_mask_NuMu=np.where(np.abs(Flavour_Shower_Neutrino[MuonEnergyThreshMask])==14., True, False)
np.save(config['flavour_mask_NuMu_path'],Flav_mask_NuMu)
Flav_mask_NuTau=np.where(np.abs(Flavour_Shower_Neutrino[MuonEnergyThreshMask])==16., True, False)
np.save(config['flavour_mask_NuTau_path'],Flav_mask_NuTau)
##my plots do not express correlated muons by energy when I use these values -- try expressly using uncorrelated muons
Muon_Energy_L2 = Muon_Energy_L2[MuonEnergyThreshMask]
Muon_Energy_L3 = Muon_Energy_L3[MuonEnergyThreshMask]
Muon_Energy_L4 = Muon_Energy_L4[MuonEnergyThreshMask]
Muon_Energy_L5 = Muon_Energy_L5[MuonEnergyThreshMask]

L2_mask=np.where(Muon_Energy_L2==0, True, False).tolist()
L3_mask=np.where(Muon_Energy_L3==0, True, False).tolist()
L4_mask=np.where(Muon_Energy_L4==0, True, False).tolist()
L5_mask=np.where(Muon_Energy_L5==0, True, False).tolist()

SingleMuonMask=L2_mask
np.save(config['single_muon_mask_path'],SingleMuonMask)
DoubleMuonMask=np.logical_and(L3_mask,np.logical_not(L2_mask))
np.save(config['double_muon_mask_path'],DoubleMuonMask)
TripleMuonMask=np.logical_and.reduce((L4_mask,np.logical_not(L3_mask),np.logical_not(L2_mask)))
np.save(config['triple_muon_mask_path'],TripleMuonMask)
QuadMuonMask=np.logical_and.reduce((L5_mask,np.logical_not(L4_mask),np.logical_not(L3_mask),np.logical_not(L2_mask)))
np.save(config['quad_muon_mask_path'],QuadMuonMask)

##DEFINE ENERGY, ZENITH, AND DEPTH BINS HERE 
E_nu_bins = eval(config['E_nu_bins'])
nu_bin_centers = eval(config['nu_bin_centers'])
E_mu_bins = eval(config['E_mu_bins'])
mu_bin_centers = eval(config['mu_bin_centers'])
angles_space= eval(config['angles'])
depth_space= eval(config['depths'])

#########REPLACE WITH PRIMARY DEPTHMC ZENITH AND DEPTH ASAP
Zenith_Shower_Neutrino = Zenith_Shower_Neutrino[MuonEnergyThreshMask]
print(len(Zenith_Shower_Neutrino))
print(angles_space)
angle_vals_prim=np.digitize(np.cos(Zenith_Shower_Neutrino),angles_space)-1

#Depth_BrightestMuon = Depth_BrightestMuon[MuonEnergyThreshMask]
Depth_Shower_Neutrino = Depth_Shower_Neutrino[MuonEnergyThreshMask]

depths_prim=np.digitize(Depth_Shower_Neutrino/1000,depth_space)

print('angle masks')
angle_masks_prim = [angle_vals_prim == i for i in range(len(angles_space))]
for i,angle in enumerate(angles_space):
    print(np.sum(angle_masks_prim[i]),',',angle) 
    np.save(config['angle_masks_base_path'] + str(np.around(angle,decimals=2))+'.npy',angle_masks_prim[i])

print('depth masks')
depth_masks_prim = [depths_prim == i for i in range(len(depth_space))]
for i,depth in enumerate(depth_space):
    print(np.sum(depth_masks_prim[i]),',',depth)
    np.save(config['depth_masks_base_path'] +str(np.around(depth,decimals=2))+'.npy',depth_masks_prim[i])
