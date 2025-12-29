#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/icetray-start
#METAPROJECT /home/zrechav/i3/icetray/build
import pyarrow
from icecube.icetray import I3Tray
from icecube import icetray, phys_services, dataio, dataclasses,MuonGun
import simweights
import pandas as pd
from pathlib import Path
import simweights
import h5py
import numpy as np
import traceback


import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", type=str, default = "/data/user/zrechav/compiled_hdf5s/airshowered_corsika/22803/22803_0000000-0000999_unweighted.hdf5", dest = "INPUT", help = "input i3 file directory")
parser.add_argument("-o", "--output", action="store", type=str, default ='/data/user/zrechav/compiled_hdf5s/airshowered_corsika/22803/22803_0000000-0000999_weighted.hdf5', dest = "OUTPUT", help = "Output i3 file")

args = parser.parse_args()

print(args)

input_file = args.INPUT
output_file = args.OUTPUT

nfiles = 1000 ##number of Corsika files used, trigger level!
print(input_file, '\n', nfiles)


def calc_weight(inp = input_file,numfiles = nfiles):
    with pd.HDFStore(inp, "r") as hdffile:
        weighter = simweights.CorsikaWeighter(hdffile, nfiles=nfiles)
        flux = simweights.GaisserH4a()
        weights = weighter.get_weights(flux)  
    return weights
    
weights = calc_weight()

with h5py.File(input_file, 'r') as existing_hdf:
    with h5py.File(output_file, 'w') as new_hdf:
        for key in existing_hdf.keys():
            existing_hdf.copy(key, new_hdf)

        new_hdf.create_dataset('weights', data=weights, maxshape=(None,), chunks=True)

print("done")
