#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/icetray-start
#METAPROJECT /home/zrechav/i3/icetray/build


from icecube.icetray import I3Tray
from icecube import icetray, phys_services, dataio, dataclasses,MuonGun
from icecube.simprod import segments
import sys, math
import copy
import glob
import numpy as np
import collections
import os
import random
import ROOT
import time
from icecube.hdfwriter import I3HDFTableService, I3HDFWriter
import argparse
sys.path.append('/cvmfs/icecube.opensciencegrid.org/users/vbasu/scripts/NewIce/')
sys.path.append('/home/vbasu/scripts/')
import meseUtils
import MCWeightUtils_Airshowers
from icecube import dataclasses

import yaml


with open('/home/zrechav/SelfVeto_Correlation_Tables/scripts/config.yaml', 'r') as yaml_file:
    config = yaml.safe_load(yaml_file)



mass_dict=    {'2212': 1, #PPlus
               '1000020040': 4,#Alpha
               '1000070140': 14,#N14
               '1000130270': 27,#Al27
               '1000260560': 56#Fe56
              }

start_time = time.asctime()
print ('Started:', start_time)

# handling of command line arguments  

parser = argparse.ArgumentParser()

corsika_filelist_length = 1

gcd = config['gcd']

parser.add_argument("-i", "--input", action="store", type=str, dest = "INPUT", help = "Output i3 file")
parser.add_argument("-o", "--output", action="store", type=str, dest = "OUTPUT", help = "Output i3 file")
parser.add_argument("-g", "--gcd", action="store", type=str, default = gcd, dest="GCD", help="GCD file for input i3 file")


args = parser.parse_args()
print(args)
#infile = args.INPUT
infile = (args.INPUT)
#infile = args.INPUT
infiles = [gcd, infile]
print(infiles)
#print(infile)
outfile = args.OUTPUT
gcdFile= args.GCD 
tray = I3Tray()

print ("Reading input file...",  args.INPUT)
tray.AddSegment(dataio.I3Reader, 'reader', FilenameList=infiles)

#########################################
## GATHER SOME MC TRUTH
#########################################   


NuTypes=[                    
                          dataclasses.I3Particle.ParticleType.NuE,
                          dataclasses.I3Particle.ParticleType.NuEBar,
                          dataclasses.I3Particle.ParticleType.NuMu,
                          dataclasses.I3Particle.ParticleType.NuMuBar,
                          dataclasses.I3Particle.ParticleType.NuTau,
                          dataclasses.I3Particle.ParticleType.NuTauBar,
                          ]
MuTypes=[         
                          dataclasses.I3Particle.ParticleType.MuPlus,
                          dataclasses.I3Particle.ParticleType.MuMinus,
                        
                          ]
gcdFile= config['gcd']

import matplotlib.path as mpltPath
##need this line of code for i3files with no pframes
def add_eventheader(frame):
    eh = dataclasses.I3EventHeader()
    eh.run_id = 1
    eh.event_id = add_eventheader.event_id
    add_eventheader.event_id += 1
    frame['I3EventHeader'] = eh
    
add_eventheader.event_id = 1
tray.Add(add_eventheader, Streams=[icetray.I3Frame.DAQ])

tray.AddModule("I3NullSplitter", "fullevent", SubEventStreamName='InIceSplit')

###surface for cylinder of partially contained boundary, with 50 meter padding on cylinder top as an extra bin --> this surface_det configuration is DIFFERENT from contained vs. partially contained id, independent of string geometry    
def get_surface_det(gcdFile=None):
    
    from icecube import MuonGun
    bound_2D=[]
    surface_det = MuonGun.Cylinder(1100,700)   
    t = np.linspace(0, 2 * np.pi, 100)
    surface_det_x = 700 * np.cos(t)
    surface_det_y = 700 * np.sin(t)
    x=[(surface_det_x[i],surface_det_y[i])for i in range(len(surface_det_x))]
    bound_2D= mpltPath.Path(x)#Projection of detector on x,y plane
    return bound_2D, surface_det
    
def boundary_check(particle1,gcdFile=None):
    ####checks if particle is inside the detector###
    gcdFile=gcdFile
    bound_2D,surface_det = get_surface_det(gcdFile=gcdFile)
    cyl_top = 550.
    cyl_bot = -500.
    inlimit = False  
    print(particle1.pos.z)
    if (((particle1.pos.z <=cyl_top) and (particle1.pos.z>=cyl_bot))) and bound_2D.contains_points([(particle1.pos.x, particle1.pos.y)]):
            inlimit=True            
    return inlimit

bound_2D,surface_det = get_surface_det()


def NeutrinoSelector(frame):
    print('IM in NEUTRINO SELECTOR')
    if 'I3MCTree' in frame:# and 'Homogenized_QTot' in frame: ##to make sure i am grabbing the in ice split pframe and no others
        print('im in da loop')
        #print(frame['FilterMask'])                                                        
        frame['I3MCTree_copy']=frame['I3MCTree']
        mctree = frame['I3MCTree_copy']
        #mctree = frame['I3MCTree_preMuonProp']
        e_dep_total=0
        neutrino = None
        neutrinos=[]
        primary=None
        for p in mctree:#search for neutrino
            ##conditional if loop --> only if there is no polyplopiaprimary information for some reason
            if mctree.depth(p)==0 and primary==None:
                primary=p
                if 'PolyplopiaPrimary' in frame:
                    frame['PolyplopiaPrimaryCopy'] = frame['PolyplopiaPrimary']
                    del frame['PolyplopiaPrimary']
                frame["PolyplopiaPrimary"]=dataclasses.I3Particle(primary)
                frame["PrimaryMass"]=dataclasses.I3Double(mass_dict[str(primary.pdg_encoding)])

            if (p.type in NuTypes) and (p.energy > 500.): ##273 propagation cutoff, 500 DNNCascades cutoff
                neutrinos.append(p)                       
        if neutrinos==[]:
            print('IM EMPTY OF NEUTRINOS, nEXT')
            return False 
        num_neutrinos = len(neutrinos)
        frame['num_neutrinos'] = dataclasses.I3Double(num_neutrinos)
        neutrino=random.choice(neutrinos)
        neutrino_parent=mctree.parent(neutrino)
        frame['ShowerNeutrino'] =dataclasses.I3Particle(neutrino)
        
        
        intersection=surface_det.intersection(neutrino.pos, neutrino.dir)#points of intersection
        z_inter=neutrino.pos.z-intersection.first*np.cos(neutrino.dir.zenith)
        n_rho = np.sqrt(neutrino.pos.x**2 + neutrino.pos.y**2)
        rho_inter = n_rho - intersection.first * np.sin(neutrino.dir.zenith)
        # Calculate x and y from rho
        x_inter = rho_inter * np.cos(neutrino.dir.azimuth)
        y_inter = rho_inter * np.sin(neutrino.dir.azimuth)
        #print(z_inter)
        depth = 1948.07 - z_inter
        #print(depth)
        frame['rho_neutrino'] = dataclasses.I3Double(n_rho)
        frame['z_inter'] = dataclasses.I3Double(z_inter)
        frame['rho_inter'] = dataclasses.I3Double(rho_inter)
        frame['x_inter'] = dataclasses.I3Double(x_inter)
        frame['y_inter'] = dataclasses.I3Double(y_inter)
        frame['shower_neutrino_depth'] = dataclasses.I3Double(1948.07 - z_inter)
        frame['shower_neutrino_energy']=dataclasses.I3Double(neutrino.energy)
        frame['shower_neutrino_zenith']=dataclasses.I3Double(neutrino.dir.zenith)
        frame['shower_neutrino_type']=dataclasses.I3Double(neutrino.type)
        frame['NeutrinoParent'] =dataclasses.I3Particle(neutrino_parent)
        print('I FINISHED SELECTING?')
        mctree.erase_children(neutrino)
                
        

        
tray.Add(NeutrinoSelector)
#required because not interested in NuMu cc muons :/
def issameparticle(particle1, particle2):
        isequal = True
        if (particle1.pos-particle2.pos).magnitude > 1e-3:
            isequal = False
        if np.abs(particle1.dir.zenith-particle2.dir.zenith) > 1e-3:
            isequal = False
        if np.abs(particle1.dir.azimuth-particle2.dir.azimuth) > 1e-3:
            isequal = False
        if np.abs(particle1.energy-particle2.energy) > 1e-3:
            isequal = False
        if np.abs(particle1.time-particle2.time) > 1e-3:
            isequal = False
        if particle1.type != particle2.type:
            isequal = False        
        return isequal
def get_lateral_separation(muon, particle2):

    # c = 2.99792458e8 #* I3Units::m / (I3Units::second)
    gamma_muon=muon.energy/0.10566 #MuonMass=105.66MeV/c^2
    p_muon=np.sqrt((muon.energy)**2-(0.10566)**2)/dataclasses.I3Constants.c
    theta=muon.dir.angle(particle2.dir)
    pt_muon=p_muon*np.sin(theta)
    muon_path_length=muon.pos.z/np.cos(muon.dir.zenith)
    lateral_separation=pt_muon*dataclasses.I3Constants.c*muon_path_length/muon.energy
    
    return lateral_separation


def get_deposit_energy(frame):
    print('IM IN DEPOSIT ENERGY')
    if 'I3MCTree_copy' in frame:
        mctree=frame['I3MCTree_copy']
        muon_multiplicity=0

        e_muon_total=0
        muon_list=[]
        muon_energy_list=[]
        for track in MuonGun.Track.harvest(mctree, frame['MMCTrackList']):
            if (boundary_check(track,gcdFile)): continue
       
            track_energy_at_det=track.get_energy(surface_det.intersection(track.pos, track.dir).first)
            if track_energy_at_det>10.: ##5 GeV minimum energy
                    e_muon_total+=track_energy_at_det
                    muon_list.append(track)
                    muon_energy_list.append(track_energy_at_det)
                    muon_multiplicity+=1

        if muon_list == []: ##adding new line here 2 Sept 2024
            print('IM EMPTY OF MUONS, nEXT')
            return False
        frame['Total_Muon_Energy'] =dataclasses.I3Double(e_muon_total)
        frame['MuonMultiplicity'] =dataclasses.I3Double(muon_multiplicity)

        ##try random sampling instead of highest energy sampling test!
        combined_list=sorted(zip(muon_energy_list, muon_list))
        #combined_list = list(zip(muon_energy_list, muon_list))
        sample_size = min(len(combined_list), 5)
        sorted_energy_particles = random.sample(combined_list, sample_size)
        #sorted_energy_particles = random.sample(sorted_energy_particles, 5)
        ctr=0
        for energy,muon in sorted_energy_particles:
            ctr+=1
            energy_name="Muon_Energy_L{}".format(ctr)
            muon_name="Muon_L{}".format(ctr)
            muon_depth_name="Muon_L{}_Depth".format(ctr)
            rho_name="Muon_L{}_rho".format(ctr)
            z_inter_name ="Muon_L{}_z_inter".format(ctr)
            y_inter_name ="Muon_L{}_y_inter".format(ctr)
            x_inter_name ="Muon_L{}_x_inter".format(ctr)
            rho_inter_name ="Muon_L{}_rho_inter".format(ctr)
            
            frame[energy_name]=dataclasses.I3Double(energy)
            frame[muon_name]=dataclasses.I3Particle(muon)
            
            intersection=surface_det.intersection(muon.pos, muon.dir)#points of intersection
            z_inter=muon.pos.z-intersection.first*np.cos(muon.dir.zenith)
            depth=1948.07-z_inter
            muon_rho = np.sqrt(muon.pos.x**2 + muon.pos.y**2)
            rho_inter = muon_rho - intersection.first * np.sin(muon.dir.zenith)
        # Calculate x and y from rho
            x_inter = rho_inter * np.cos(muon.dir.azimuth)
            y_inter = rho_inter * np.sin(muon.dir.azimuth)
        
            frame[rho_name] = dataclasses.I3Double(muon_rho)
            frame[z_inter_name] = dataclasses.I3Double(z_inter)
            frame[rho_inter_name] = dataclasses.I3Double(rho_inter)
            frame[x_inter_name] = dataclasses.I3Double(x_inter)
            frame[y_inter_name] = dataclasses.I3Double(y_inter)
            
            frame[muon_depth_name]=dataclasses.I3Double(depth)
            muon_sep_neut_name="Muon_L{}_Neutrino_separation".format(ctr)
            muon_sep_shower_name="Muon_L{}_Shower_separation".format(ctr)
   
            sep1=get_lateral_separation(muon,frame["ShowerNeutrino"])
            frame[muon_sep_neut_name]=dataclasses.I3Double(sep1)
            sep2=get_lateral_separation(muon,frame["PolyplopiaPrimary"])
            frame[muon_sep_shower_name]=dataclasses.I3Double(sep2)
            print('I FINISHED MUON SELECTING')
tray.Add(get_deposit_energy)


###########book keys for hdf file###########
table_keys=[
            'GaisserH4a_weight',
            #"HomogenizedQTot",
            "OneWeight",
            "NEvents",
            "PrimaryMass",
            'I3PrimaryInjectorInfo',
            "PolyplopiaPrimary",
            'PolyplopiaInfo',
            "I3MCWeightDict",
            "CorsikaWeightMap",
            'PrimaryFlavour',
            #'I3CorsikaInfo',
            'Interaction_Type',
            'OneWeight',
            'MuonWeight',
            'MuonWeight2',
            "MCPrimaryType",
            "MCPrimary",
            'Muon_Energy_L1',
            'Muon_Energy_L2',
            'Muon_Energy_L3',
            'Muon_Energy_L4',
            'Muon_Energy_L5',
            'Muon_L1',
            'Muon_L2',
            'Muon_L3',
            'Muon_L4',
            'Muon_L5',
            'Muon_L1_Depth',
            'Muon_L2_Depth',
            'Muon_L3_Depth',
            'Muon_L4_Depth',
            'Muon_L5_Depth',
            'Muon_L1_Neutrino_separation',
            'Muon_L2_Neutrino_separation',
            'Muon_L3_Neutrino_separation',
            'Muon_L4_Neutrino_separation',
            'Muon_L5_Neutrino_separation',
            'Muon_L1_Shower_separation',
            'Muon_L2_Shower_separation',
            'Muon_L3_Shower_separation',
            'Muon_L4_Shower_separation',
            'Muon_L5_Shower_separation',
            'MuonMultiplicity',
            
            'ShowerNeutrino',
            'Total_Muon_Energy',
            'BrightestMuonAtDepth',
            'Energy_BrightestMuonAtDepth',
            "I3EventHeader",
            "NEvents",
            'PolyplopiaPrimaryCopy',
            'z_inter',
            'shower_neutrino_depth',
            'shower_neutrino_energy',
            'shower_neutrino_zenith',
            'shower_neutrino_type',
            #'NeutrinoParent',
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
            'Muon_L5_rho',
            'Muon_L5_rho_inter',
            'Muon_L5_z_inter',
            'Muon_L5_y_inter',
            'Muon_L5_x_inter',
            #'I3MCTree'
           ]

tray.AddModule('I3Writer', 'writer',
    DropOrphanStreams=[icetray.I3Frame.DAQ],
    Streams=[  icetray.I3Frame.DAQ, icetray.I3Frame.Physics],
    filename=outfile+'.i3.zst')

# tray.Add(I3HDFWriter,'HDFwriter',
#               Output=outfile+".hdf5",
#               keys         = table_keys,
#               SubEventStreams = ['InIceSplit']
# )
tray.AddModule('TrashCan', 'thecan')

tray.Execute()
tray.Finish()

del tray

stop_time = time.asctime()

print ('Started:', start_time)
print ('Ended:', stop_time)

