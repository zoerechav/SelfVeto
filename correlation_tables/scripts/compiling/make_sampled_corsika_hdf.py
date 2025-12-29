#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/icetray-start
#METAPROJECT /home/zrechav/i3/icetray/build

from pathlib import Path
from icecube.icetray import I3Tray
from icecube import snowstorm, icetray, dataio, dataclasses, hdfwriter, simclasses,phys_services
from icecube.dataclasses import I3MapStringDouble
from icecube.icetray import traysegment, I3Module
from icecube.simprod import segments
from icecube.photonics_service import I3PhotoSplineService
from I3Tray import *
from icecube import MuonGun, simclasses
import numpy as np
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", type=str, default = "/data/user/zrechav/airshowered_corsika/22803/0000000-0000999", dest = "INPUT", help = "input i3 file directory")
parser.add_argument("-o", "--output", action="store", type=str, default ='/data/user/zrechav/compiled_hdf5s/airshowered_corsika/22803/22803_0000000-0000999_unweighted.hdf5', dest = "OUTPUT", help = "Output i3 file")

args = parser.parse_args()

print(args)

infile = args.INPUT
outfile = args.OUTPUT
FILE_DIR = Path(infile)

files = sorted(str(f) for f in FILE_DIR.glob("*.i3.zst"))
print(files)
tray = I3Tray()
tray.Add("I3Reader", FileNameList=files)



#def get_total_neutrino_energy(frame):
#    if 'I3MCTree' in frame:
#        mctree = frame['I3MCTree']
#        energy = 0
#        for p in mctree:
#            if p.type in [p.NuE, p.NuEBar, p.NuMu, p.NuMuBar,p.NuTau,p.NuTauBar] and p.energy > 500.:
#                 energy += p.energy#
#
#        frame['total_neutrino_energy'] = dataclasses.I3Double(energy)
#tray.Add(get_total_neutrino_energy,'get_total_neutrino_energy')


def Frame_Yeeter(frame):
    if frame.Has("Total_Muon_Energy") and frame.Has('I3EventHeader') and frame.Has('Muon_L1'):
        return True  # Keep the frame
    else:
        return False  # Kick the frame out

tray.Add(Frame_Yeeter,'Frame_Yeeter')


try:
    tray.Add(
        hdfwriter.I3HDFWriter,
        SubEventStreams=["InIceSplit"],
        keys=[
            "Homogenized_QTot",
            "OneWeight",
            "NEvents",
            "PrimaryMass",
            'I3PrimaryInjectorInfo',
            "PolyplopiaPrimary",
            "I3MCWeightDict",
            "CorsikaWeightMap",
            'OneWeight',        
            "MCPrimaryType",
            "MCPrimary",
            "I3EventHeader",
            "intersection_rho_small",
            "intersection_rho_large",
            "intersection_rho_bound",
            'PolyplopiaInfo',
            'MuonMultiplicity',
            'ShowerNeutrino',
            'Total_Muon_Energy',
            'shower_neutrino_depth',
            'shower_neutrino_energy',
            'shower_neutrino_type',
            'shower_neutrino_zenith',
            'z_inter ',
            'Muon_Energy_L1',
            'Muon_L1_Depth',
            'Muon_L1_Neutrino_separation',
            'Muon_L1_Shower_separation',
            'Muon_L1',
            'Muon_Energy_L2',
            'Muon_L2_Depth',
            'Muon_L2_Neutrino_separation',
            'Muon_L2_Shower_separation',
            'Muon_L2',
            'Muon_Energy_L3',
            'Muon_L3_Depth',
            'Muon_L3_Neutrino_separation',
            'Muon_L3_Shower_separation',
            'Muon_L3',
            'Muon_Energy_L4',
            'Muon_L4_Depth',
            'Muon_L4_Neutrino_separation',
            'Muon_L4_Shower_separation',
            'Muon_L4',
            #'I3MCTree',
            'total_neutrino_energy',
            'z_inter',
           
           
            'NeutrinoParent',
            'num_neutrinos',
    
            ###new keys to add
            'rho_neutrino',   
            'rho_inter',
            'x_inter',
            'y_inter',
            'Muon_L1_rho',
            'Muon_L1_rho_inter',
            'Muon_L1_z_inter',
            'Muon_L1_y_inter',
            'Muon_L1_x_inter',
            'Muon_L2_rho',
            'Muon_L2_rho_inter',
            'Muon_L2_z_inter',
            'Muon_L2_y_inter',
            'Muon_L2_x_inter',
            'Muon_L3_rho',
            'Muon_L3_rho_inter',
            'Muon_L3_z_inter',
            'Muon_L3_y_inter',
            'Muon_L3_x_inter',
            'Muon_L4_rho',
            'Muon_L4_rho_inter',
            'Muon_L4_z_inter',
            'Muon_L4_y_inter',
            'Muon_L4_x_inter',
            #'Muon_L5_rho',
            #'Muon_L5_rho_inter',
            #'Muon_L5_z_inter',
            #'Muon_L5_y_inter',
            #'Muon_L5_x_inter',
             
            
        ],
        output=outfile,
    )

    tray.Execute()

except Exception as e:
    #print(f"Error occurred while processing file: {files[tray.curframe]}")
    print(f"Error message: {str(e)}")
#continue