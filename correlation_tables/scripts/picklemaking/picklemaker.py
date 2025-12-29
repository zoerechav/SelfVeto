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
        #Muon_Zenith_L5=np.append(Muon_Zenith_L5,np.asarray(hdf.get('Muon_L5')['zenith']))
        
        #Muon_Neutrino_separation_L1 = np.append(Muon_Neutrino_separation_L1, np.asarray(hdf.get('Muon_L1_Neutrino_separation')['value']))
        #Muon_Neutrino_separation_L2 = np.append(Muon_Neutrino_separation_L2, np.asarray(hdf.get('Muon_L2_Neutrino_separation')['value']))
        #Muon_Neutrino_separation_L3 = np.append(Muon_Neutrino_separation_L3, np.asarray(hdf.get('Muon_L3_Neutrino_separation')['value']))
        #Muon_Neutrino_separation_L4 = np.append(Muon_Neutrino_separation_L4, np.asarray(hdf.get('Muon_L4_Neutrino_separation')['value']))
        #Muon_Neutrino_separation_L5 = np.append(Muon_Neutrino_separation_L5, np.asarray(hdf.get('Muon_L5_Neutrino_separation')['value']))
        print('almost therrrrr')
        #Muon_Shower_separation_L1 = np.append(Muon_Shower_separation_L1, np.asarray(hdf.get('Muon_L1_Shower_separation')['value']))
        #Muon_Shower_separation_L2 = np.append(Muon_Shower_separation_L2, np.asarray(hdf.get('Muon_L2_Shower_separation')['value']))
        #Muon_Shower_separation_L3 = np.append(Muon_Shower_separation_L3, np.asarray(hdf.get('Muon_L3_Shower_separation')['value']))
        #Muon_Shower_separation_L4 = np.append(Muon_Shower_separation_L4, np.asarray(hdf.get('Muon_L4_Shower_separation')['value']))
        #Muon_Shower_separation_L5 = np.append(Muon_Shower_separation_L5, np.asarray(hdf.get('Muon_L5_Shower_separation')['value']))

        Total_Muon_energy = np.append(Total_Muon_energy,np.asarray(hdf.get('Total_Muon_Energy')['value']))
        print('I FINISHED LOADING MY ARRAYS WITH HDF INFO')
        hdf.close()
    except Exception as e:
        print(filename+" Is Faulty")
        print(e)
        print(f"Occurred on line: {traceback.extract_tb(e.__traceback__)[-1][1]}")

        hdf.close()

##MASKING AND SAVING INFORMATION TO SAVE TIME LATER ON                
            
MuonEnergyThreshMask=[energy>10. for energy in Total_Muon_energy]
print('MuonEnergyThreshMask : ', len(MuonEnergyThreshMask))


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
#Muon_Energy_L5 = Muon_Energy_L5[MuonEnergyThreshMask]

L2_mask=np.where(Muon_Energy_L2==0, True, False).tolist()
L3_mask=np.where(Muon_Energy_L3==0, True, False).tolist()
L4_mask=np.where(Muon_Energy_L4==0, True, False).tolist()
#L5_mask=np.where(Muon_Energy_L5==0, True, False).tolist()

SingleMuonMask=L2_mask
np.save(config['single_muon_mask_path'],SingleMuonMask)
DoubleMuonMask=np.logical_and(L3_mask,np.logical_not(L2_mask))
np.save(config['double_muon_mask_path'],DoubleMuonMask)
TripleMuonMask=np.logical_and.reduce((L4_mask,np.logical_not(L3_mask),np.logical_not(L2_mask)))
np.save(config['triple_muon_mask_path'],TripleMuonMask)
QuadMuonMask=np.logical_and.reduce((np.logical_not(L4_mask),np.logical_not(L3_mask),np.logical_not(L2_mask)))
np.save(config['quad_muon_mask_path'],QuadMuonMask)

##DEFINE ENERGY, ZENITH, AND DEPTH BINS HERE 
E_nu_bins = eval(config['E_nu_bins'])
nu_bin_centers = eval(config['nu_bin_centers'])
E_mu_bins = eval(config['E_mu_bins'])
mu_bin_centers = eval(config['mu_bin_centers'])

nue_angles_space= eval(config['nue_angles'])
numu_angles_space= eval(config['numu_angles'])
depth_space= eval(config['depths'])


#########REPLACE WITH PRIMARY DEPTHMC ZENITH AND DEPTH ASAP
Zenith_Shower_Neutrino = Zenith_Shower_Neutrino[MuonEnergyThreshMask]
NuE_Zenith_Shower_Neutrino = Zenith_Shower_Neutrino[Flav_mask_NuE]
NuMu_Zenith_Shower_Neutrino = Zenith_Shower_Neutrino[Flav_mask_NuMu]
print(len(Zenith_Shower_Neutrino))
print(len(NuE_Zenith_Shower_Neutrino))
print(len(NuMu_Zenith_Shower_Neutrino))

#print(angles_space)


#Depth_BrightestMuon = Depth_BrightestMuon[MuonEnergyThreshMask]
Depth_Shower_Neutrino = Depth_Shower_Neutrino[MuonEnergyThreshMask]

depths_prim=np.digitize(Depth_Shower_Neutrino/1000,depth_space)

print('NuE angle masks')
nue_angle_vals_prim=np.digitize(np.cos(Zenith_Shower_Neutrino),nue_angles_space)-1
nue_angle_masks_prim = [nue_angle_vals_prim == i for i in range(len(nue_angles_space))]
for i,angle in enumerate(nue_angles_space):
    print(np.sum(nue_angle_masks_prim[i]),',',angle) 
    np.save(config['angle_masks_base_path'] + 'nue_'+
    str(np.around(angle,decimals=2))+'.npy',nue_angle_masks_prim[i])
    

print('NuMu angle masks')
numu_angle_vals_prim=np.digitize(np.cos(Zenith_Shower_Neutrino),numu_angles_space)-1
numu_angle_masks_prim = [numu_angle_vals_prim == i for i in range(len(numu_angles_space))]
for i,angle in enumerate(numu_angles_space):
    print(np.sum(numu_angle_masks_prim[i]),',',angle) 
    np.save(config['angle_masks_base_path'] + 'numu_'+ 
    str(np.around(angle,decimals=2))+'.npy',numu_angle_masks_prim[i])
# angle_vals_prim=np.digitize(np.cos(Zenith_Shower_Neutrino),angles_space)-1
# angle_masks_prim = [angle_vals_prim == i for i in range(len(angles_space))]
# for i,angle in enumerate(angles_space):
#     print(np.sum(angle_masks_prim[i]),',',angle) 
#     np.save(config['angle_masks_base_path'] + str(np.around(angle,decimals=2))+'.npy',angle_masks_prim[i])

print('depth masks')
depth_masks_prim = [depths_prim == i for i in range(len(depth_space))]
for i,depth in enumerate(depth_space):
    print(np.sum(depth_masks_prim[i]),',',depth)
    np.save(config['depth_masks_base_path'] +str(np.around(depth,decimals=2))+'.npy',depth_masks_prim[i])
