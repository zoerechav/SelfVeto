#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/icetray-start
#METAPROJECT /cvmfs/icecube.opensciencegrid.org/users/vbasu/meta-projects/combo2/build
import numpy as np
import math
import os.path
from os import path

import tables
import h5py
import numpy as np
import math
import os.path
from os import path

import tables
import glob
import h5py
import pandas as pd

import yaml
import argparse
import os

import simweights
 
with open('/home/zrechav/SelfVeto_Correlation_Tables/scripts/config.yaml', 'r') as yaml_file:
    config = yaml.safe_load(yaml_file)


E_nu_bins=eval(config['E_nu_bins'])
nu_bin_centers = eval(config['nu_bin_centers'])
E_mu_bins=eval(config['E_mu_bins'])
mu_bin_centers = eval(config['mu_bin_centers'])


SingleMuonMask=np.load(config['single_muon_mask_path'],mmap_mode='r',allow_pickle=True)


Muon_Energy_L1 = np.array([])
Energy_Shower_Neutrino = np.array([])   
Zenith_Shower_Neutrino = np.array([])   
Depth_Shower_Neutrino = np.array([])   
Flavour_Shower_Neutrino = np.array([])
GaisserH4a_weight = np.array([])

filename = config['corsika_sample']

nfiles = eval(config['nfiles'])

def calc_weight(inp = filename,numfiles = nfiles):
    with pd.HDFStore(inp, "r") as hdffile:
        weighter = simweights.CorsikaWeighter(hdffile, nfiles=nfiles)
        flux = simweights.GaisserH4a()
        weights = weighter.get_weights(flux)  
    return weights

if (path.exists(filename)):
    try:
        print('I EXIST')
        GaisserH4a_weight = calc_weight()
        hdf=h5py.File(filename, 'r')

        #GaisserH4a_weight = np.append(GaisserH4a_weight,np.asarray(hdf.get('GaisserH4a_weight')['item']))
        Zenith_Shower_Neutrino = np.append(Zenith_Shower_Neutrino,np.asarray(hdf.get('shower_neutrino_zenith')['value']))
        Flavour_Shower_Neutrino = np.append(Flavour_Shower_Neutrino,np.asarray(hdf.get('shower_neutrino_type')['value']))
        Energy_Shower_Neutrino = np.append(Energy_Shower_Neutrino,np.asarray(hdf.get('shower_neutrino_energy')['value']))
        Depth_Shower_Neutrino = np.append(Depth_Shower_Neutrino,np.asarray(hdf.get('shower_neutrino_depth')['value']))
                           
        Muon_Energy_L1=np.append(Muon_Energy_L1,np.asarray(hdf.get('Muon_Energy_L1')['value']))
       
        hdf.close()
    except Exception as e:
        print(filename+" Is Faulty")
        print(e)

        hdf.close()
        
MuonEnergyThreshMask=[energy>10. for energy in Muon_Energy_L1]
Energy_Shower_Neutrino=np.float32(Energy_Shower_Neutrino[MuonEnergyThreshMask])
Flavour_Shower_Neutrino=np.float32(Flavour_Shower_Neutrino[MuonEnergyThreshMask])
Muon_Energy_L1 = Muon_Energy_L1[MuonEnergyThreshMask]
GaisserH4a_weight = GaisserH4a_weight[MuonEnergyThreshMask]


flavours = (config['flavours'])
angles= eval(config['angles'])
depths = eval(config['depths'])

print(angles)
for flavour in flavours:
    Flav_mask = np.load(config['flavour_masks_base_path'] + str(flavour) +'.npy', mmap_mode = 'r', allow_pickle=True)
    for angle in angles:
        angle_mask = np.load(config['angle_masks_base_path'] + str(np.round(angle,2)) + '.npy', mmap_mode = 'r', allow_pickle=True)
        for depth in depths:
            depth_mask = np.load(config['depth_masks_base_path'] + str(np.round(depth,2)) + '.npy', mmap_mode='r',allow_pickle=True)

            print(flavour,angle,depth)
            MultiplicityMasks=[np.logical_and(SingleMuonMask,Flav_mask)]
            for multmask in MultiplicityMasks:
                ADmask=np.logical_and.reduce((angle_mask,depth_mask,multmask))
                print('AD mask sum',np.sum(ADmask))

                if np.sum(ADmask)>0:
                    print('I AM IN HERE')
                    Hist2D,xedges,yedges = np.histogram2d(Energy_Shower_Neutrino[ADmask], Muon_Energy_L1[ADmask], bins = [E_nu_bins,E_mu_bins], weights = GaisserH4a_weight[ADmask])
                    Hist2D = np.nan_to_num(Hist2D)
                    ##LEAVE HISTOGRAM UNNORMALIZED, WILL BE NORMALIZED LATER
                    #Hist2D= Hist2D/np.sum(Hist2D)
                    filename = config['single_hists_base'] + flavour + '_Single_Zen_' + str(np.round(angle,2)) + '_Depth_' + str(np.round(depth,2)) + '.npy'
                    np.save(filename, Hist2D)
                else:
                    print('Empty array')
    