#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/icetray-start
#METAPROJECT /home/zrechav/i3/icetray/build

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
 
parser = argparse.ArgumentParser()


with open('/home/zrechav/SelfVeto_Correlation_Tables/scripts/config.yaml', 'r') as yaml_file:
    config = yaml.safe_load(yaml_file)


E_nu_bins=eval(config['E_nu_bins'])
nu_bin_centers = eval(config['nu_bin_centers'])
E_mu_bins=eval(config['E_mu_bins'])
mu_bin_centers = eval(config['mu_bin_centers'])
#mult_bins = eval(config['mult_bins'])

Muon_Energy_L1 = np.array([])

Energy_Shower_Neutrino = np.array([])   
Zenith_Shower_Neutrino = np.array([])   
Depth_Shower_Neutrino = np.array([])   
Flavour_Shower_Neutrino = np.array([])
GaisserH4a_weight = np.array([])
Total_Muon_energy = np.array([])
MuonMultiplicity = np.array([])

import argparse

parser = argparse.ArgumentParser()

default_filename = config['corsika_sample']
parser.add_argument("-i", "--input", action="store", type=str, default = default_filename, dest = "INPUT", help = "input i3 file directory")

args = parser.parse_args()

print(args)

#filename = path(args.INPUT)
filename = config['corsika_sample']

print(filename)

nfiles = eval(config['nfiles'])
print(nfiles)

def calc_weight(inp = filename,numfiles = nfiles):
    with pd.HDFStore(inp, "r") as hdffile:
        weighter = simweights.CorsikaWeighter(hdffile, nfiles=nfiles)
        flux = simweights.GaisserH4a()
        weights = weighter.get_weights(flux)  
    return weights

if (path.exists(filename)):
    try:
        GaisserH4a_weight = calc_weight()
        print('weight calced')
        hdf=h5py.File(filename, 'r')
        
        #GaisserH4a_weight = np.append(GaisserH4a_weight,np.asarray(hdf.get('GaisserH4a_weight')['item']))
        Zenith_Shower_Neutrino = np.append(Zenith_Shower_Neutrino,np.asarray(hdf.get('shower_neutrino_zenith')['value']))
        Flavour_Shower_Neutrino = np.append(Flavour_Shower_Neutrino,np.asarray(hdf.get('shower_neutrino_type')['value']))
        Energy_Shower_Neutrino = np.append(Energy_Shower_Neutrino,np.asarray(hdf.get('shower_neutrino_energy')['value']))
        Depth_Shower_Neutrino = np.append(Depth_Shower_Neutrino,np.asarray(hdf.get('shower_neutrino_depth')['value']))

        Muon_Energy_L1=np.append(Muon_Energy_L1,np.asarray(hdf.get('Muon_Energy_L1')['value']))
        
        MuonMultiplicity=np.append(MuonMultiplicity,np.asarray(hdf.get('MuonMultiplicity')['value']))
        print('Im loaded')
        hdf.close()
    except Exception as e:
        print(filename+" Is Faulty")
        print(e)

        hdf.close()

    
mult_bins=eval(config['mult_bins'])
print(len(MuonMultiplicity))

MuonEnergyThreshMask=[energy>10. for energy in Muon_Energy_L1]

Energy_Shower_Neutrino=np.float32(Energy_Shower_Neutrino[MuonEnergyThreshMask])
Flavour_Shower_Neutrino=np.float32(Flavour_Shower_Neutrino[MuonEnergyThreshMask])
Muon_Energy_L1 = Muon_Energy_L1[MuonEnergyThreshMask]
GaisserH4a_weight = GaisserH4a_weight[MuonEnergyThreshMask]
MuonMultiplicity = MuonMultiplicity[MuonEnergyThreshMask]


flavours = config['flavours']
angles= eval(config['angles'])
depths = eval(config['depths'])


#for flavour in flavours:
#    flav_mask = np.load(config['flavour_masks_base_path']+str(flavour)+'.npy',mmap_mode='r',allow_pickle=True)
for angle in angles:
    angle_mask = np.load(config['angle_masks_base_path']+str(np.around(angle,2))+'.npy',mmap_mode='r',allow_pickle=True)
    for depth in depths:
        depth_mask = np.load(config['depth_masks_base_path']+str(np.around(depth,2))+'.npy',mmap_mode='r',allow_pickle=True)

        ADmask=np.logical_and.reduce((angle_mask,depth_mask))

        if np.sum(ADmask)>0:
            stackarray=np.stack((MuonMultiplicity[ADmask],Energy_Shower_Neutrino[ADmask]),axis=1)
            Hist2D,edges=np.histogramdd(stackarray,bins=(mult_bins,E_nu_bins),weights=GaisserH4a_weight[ADmask])
            #Hist2D = Hist2D /np.sum(Hist2D)
            Hist2D = np.nan_to_num(Hist2D)
            Hist2D = Hist2D.T
            filename = config['multi_base'] + 'Mult_Zen_' + str(np.around(angle,2)) + '_Depth_' + str(np.around(depth,2))+'.npy'

            print('Zen_'+str(angle)+'_Depth_'+str(depth)+' '+str(np.sum(Hist2D)))
            np.save(filename, Hist2D)
        else:
            print('Empty array')
