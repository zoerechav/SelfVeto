from I3Tray import *
from icecube import icetray, dataio, dataclasses
import sys, math
import copy
import glob
from icecube.simprod import segments

import os
from icecube.dataclasses import I3Double, I3Particle, I3Direction, I3Position, I3VectorI3Particle, I3Constants, I3VectorOMKey
from icecube.simclasses import I3MMCTrack
from icecube import MuonGun, simclasses,millipede
import numpy as np
from icecube.photonics_service import I3PhotoSplineService

from icecube.dataclasses import I3MapStringDouble
from icecube import snowstorm
from icecube.hdfwriter import I3HDFTableService, I3HDFWriter

import time
#sys.path.append('/cvmfs/icecube.opensciencegrid.org/users/vbasu/scripts/NewIce/')

import collections
import matplotlib.path as mpltPath

start_time = time.asctime()
print ('Started:', start_time)
# handling of command line arguments  
from optparse import OptionParser
parser = OptionParser()
usage = """%prog [options]"""
parser.set_usage(usage)
parser.add_option("-i", "--input", action="store", type="string", dest="INPUT", help="Input i3 file")
parser.add_option("-o", "--output", action="store", type="string", dest="OUTPUT", help="Output i3 file", default = '/home/zrechav/test/test')


# parse cmd line args, bail out if anything is not understood
(options,args) = parser.parse_args()

FileInput=options.INPUT
print(FileInput)
FileOutput=options.OUTPUT
print(FileOutput)
##omit first two characters, due to submit and shell script formatting
FileInput=options.INPUT[2:]
FileOutput=options.OUTPUT[2:]

infiles = (sorted(glob.glob(FileInput)))
print(FileInput)
print(len(infiles))
#print(infiles.insert(0,options.GCD))
outfile =FileOutput



tray = I3Tray()

#tray.Add("Rename","rename", Keys = ['MuonWeight', 'MuonWeight_Zoe'])

tray.AddSegment(dataio.I3Reader, 'reader', FilenameList=infiles)
    
print ("Reading input file...", len(infiles))


from icecube.dataclasses import I3MapStringDouble
from icecube import snowstorm, phys_services
    
table_keys = [
    "CorsikaWeightMap",
    "I3EventHeader",
    "MuonMultiplicity",
    "Muon_Energy_L1",
    "Muon_Energy_L2",
    "Muon_Energy_L3",
    "Muon_Energy_L4",
    "Muon_L1",
    "Muon_L1_Depth",
    "Muon_L1_Neutrino_separation",
    "Muon_L1_Shower_separation",
    "Muon_L1_rho",
    "Muon_L1_rho_inter",
    "Muon_L1_x_inter",
    "Muon_L1_y_inter",
    "Muon_L1_z_inter",
    "Muon_L2",
    "Muon_L2_Depth",
    "Muon_L2_Neutrino_separation",
    "Muon_L2_Shower_separation",
    "Muon_L2_rho",
    "Muon_L2_rho_inter",
    "Muon_L2_x_inter",
    "Muon_L2_y_inter",
    "Muon_L2_z_inter",
    "Muon_L3",
    "Muon_L3_Depth",
    "Muon_L3_Neutrino_separation",
    "Muon_L3_Shower_separation",
    "Muon_L3_rho",
    "Muon_L3_rho_inter",
    "Muon_L3_x_inter",
    "Muon_L3_y_inter",
    "Muon_L3_z_inter",
    "Muon_L4",
    "Muon_L4_Depth",
    "Muon_L4_Neutrino_separation",
    "Muon_L4_Shower_separation",
    "Muon_L4_rho",
    "Muon_L4_rho_inter",
    "Muon_L4_x_inter",
    "Muon_L4_y_inter",
    "Muon_L4_z_inter",
    "Muon_L5_rho",
    "Muon_L5_rho_inter",
    "Muon_L5_x_inter",
    "Muon_L5_y_inter",
    "Muon_L5_z_inter",
    "NeutrinoParent",
    "PolyplopiaPrimary",
    "PrimaryMass",
    "ShowerNeutrino",
    "Total_Muon_Energy",
    "__I3Index__",
    "num_neutrinos",
    "rho_inter",
    "rho_neutrino",
    "shower_neutrino_depth",
    "shower_neutrino_energy",
    "shower_neutrino_type",
    "shower_neutrino_zenith",
    "total_neutrino_energy",
    "weights",
    "x_inter",
    "y_inter",
    "z_inter"
]

    
print(table_keys)
tray.Add(I3HDFWriter,'HDFwriter',
               Output=FileOutput+".hdf5",
               keys         = table_keys,
               SubEventStreams = ['topological_split','InIceSplit']
)
print('I wrote to an hdf file')

print(outfile)

tray.AddModule('TrashCan', 'thecan')

tray.Execute()
tray.Finish()


del tray

stop_time = time.asctime()

print ('Started:', start_time)
print ('Ended:', stop_time)
