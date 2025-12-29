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
#print ("Reading input file...",  options.INPUT)  
#print(options.GCD)


##need this line of code for i3files with no pframes
# def add_eventheader(frame):
#     eh = dataclasses.I3EventHeader()
#     eh.run_id = 1
#     eh.event_id = add_eventheader.event_id
#     add_eventheader.event_id += 1
#     frame['I3EventHeader'] = eh
    
# add_eventheader.event_id = 1
#tray.Add(add_eventheader, Streams=[icetray.I3Frame.DAQ])

#tray.AddModule("I3NullSplitter", "fullevent", SubEventStreamName='InIceSplit')
#tray.AddModule("I3NullSplitter",'InIceSplit')
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
def select(geometry):
        r"""Select IceCube DOMs.
        Select all DOMs with an OM type `IceCube` from the given
        detector geometry and sort the selected DOMs per string in
        ascending depth.
        Parameters
        ----------
        geometry : I3OMGeoMap or dict(OMKey, tuple(I3OMGeo, ...))
            Detector geometry
        Returns
        -------
        dict(int, list(tuple(OMKey, I3OMGeo)))
            Mapping of string numbers to sequences of IceCube DOMs
            arranged in ascending depth.
        """
        strings = collections.defaultdict(list)
        # print (type(geometry))
        for omkey, omgeo in geometry.items():
            if np.iterable(omgeo):
                omgeo = omgeo[0]

            if omgeo.omtype == dataclasses.I3OMGeo.IceCube:
                strings[omkey.string].append((omkey, omgeo))

        for doms in strings.values():
            doms.sort(
                key=lambda omgeo: omgeo[1].position.z, reverse=True)

        #print(strings)
        return strings
    
    
def boundaries(geometry):
#         Side and top boundaries
#         Find the veto's side and top boundaries.
#         Parameters
#         ----------
#         geometry : I3OMGeoMap or dict(OMKey, tuple(I3OMGeo, ...))
#             IC79 or IC86 detector geometry
#         Returns

#         -------
#         sides : set(int)
#             Sequence of string numbers of the outermost strings
#         top : float
#             Depth in detector coordinates of the first DOM on the
#             deepest non-DeepCore string minus the thickness given
#             by `top_layer`
        
        top_layer=90.*icetray.I3Units.m,
        dust_layer=(-220.*icetray.I3Units.m,-100.*icetray.I3Units.m)
        strings = select(geometry)
        top = min(strings[s][0][1].position.z for s in strings if s <= 78)
        dmax = 160.*icetray.I3Units.m

        string_pos=[]

        for string in strings:
            pos = strings[string][0][1].position
            string_pos.append([pos.x,pos.y,pos.z])
            
        #I chose these strings, 2nd outermost layer
        manual_sides = [9,10,11,12,20,29,39,49,58,66,65,64,71,70,
                       69,61,52,42,32,23,15,8] 
        boundary_x=[]
        boundary_y=[]

        
        for side_string in manual_sides:
            pos=strings[side_string][0][1].position
            boundary_x.append(pos.x)
            boundary_y.append(pos.y)
        boundary_x.append(boundary_x[0])
        boundary_y.append(boundary_y[0])
        return boundary_x,boundary_y
    
def get_surface_det_og(gcdFile=None):
    
    from icecube import MuonGun
    gcdFile=config['gcd']
    bound_2D=[]
    MuonGunGCD= gcdFile
    surface_det = MuonGun.ExtrudedPolygon.from_file(MuonGunGCD, padding=0)##Build Polygon from I3Geometr
    f = dataio.I3File(MuonGunGCD)
    omgeo = f.pop_frame(icetray.I3Frame.Geometry)['I3Geometry'].omgeo
    surface_det_x,surface_det_y=boundaries(omgeo)
    x=[(surface_det_x[i],surface_det_y[i])for i in range(len(surface_det_x))]
    bound_2D= mpltPath.Path(x)#Projection of detector on x,y plane
    #print(bound_2D)
    return bound_2D, surface_det


def extend_diagonally(point, radial_limit=700):
    x, y = point
    current_radius = np.sqrt(x**2 + y**2)
    if current_radius == 0:
        return [radial_limit, 0]  # Special case for origin
    scale_factor = radial_limit / current_radius
    new_x = x * scale_factor
    new_y = y * scale_factor
    return [new_x, new_y]

def get_surface_det(gcdFile=None):
    
    from icecube import MuonGun
    gcdFile=config['gcd']
    bound_2D=[]
    MuonGunGCD= gcdFile
    surface_det = MuonGun.ExtrudedPolygon.from_file(MuonGunGCD, padding=0)##Build Polygon from I3Geometr
    f = dataio.I3File(MuonGunGCD)
    omgeo = f.pop_frame(icetray.I3Frame.Geometry)['I3Geometry'].omgeo
    surface_det_x,surface_det_y=boundaries(omgeo)
    original_points = np.column_stack((surface_det_x, surface_det_y))
    extended_points = np.apply_along_axis(extend_diagonally, axis=1, arr=original_points)
    surface_det_x = extended_points[:, 0]
    surface_det_y = extended_points[:, 1]
    x=[(surface_det_x[i],surface_det_y[i])for i in range(len(surface_det_x))]
    bound_2D= mpltPath.Path(x)#Projection of detector on x,y plane
    #print(bound_2D)
    return bound_2D, surface_det

def boundary_check(particle1,gcdFile=None):
    ####checks if particle is inside the detector###
    gcdFile=gcdFile
    bound_2D,surface_det = get_surface_det(gcdFile=gcdFile)
    inlimit = False  
    if (((particle1.pos.z <=max(surface_det.z)) and (particle1.pos.z>=min(surface_det.z)))) and bound_2D.contains_points([(particle1.pos.x, particle1.pos.y)]):
            inlimit=True       
    #frame["In_Boundary"] = icetray.I3Bool(inlimit)
    return inlimit

bound_2D,surface_det = get_surface_det()


def NeutrinoSelector(frame):
    print('IM in NEUTRINO SELECTOR')
    if 'I3MCTree' in frame and 'Homogenized_QTot' in frame: ##to make sure i am grabbing the in ice split pframe and no others
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
            return False 
        num_neutrinos = len(neutrinos)
        frame['num_neutrinos'] = dataclasses.I3Double(num_neutrinos)
        neutrino=random.choice(neutrinos)
        neutrino_parent=mctree.parent(neutrino)
        frame['ShowerNeutrino'] =dataclasses.I3Particle(neutrino)
        
        #surface = MuonGun.Cylinder(1000,700) ###change!!!!!
        intersection=surface_det.intersection(neutrino.pos, neutrino.dir)#points of intersection
        z_inter=neutrino.pos.z-intersection.first*np.cos(neutrino.dir.zenith)
        #print(z_inter)
        depth = 1948.07 - z_inter
        #print(depth)
        frame['z_inter'] = dataclasses.I3Double(z_inter)
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
            return False
        frame['Total_Muon_Energy'] =dataclasses.I3Double(e_muon_total)
        frame['MuonMultiplicity'] =dataclasses.I3Double(muon_multiplicity)

        
        sorted_energy_particles=sorted(zip(muon_energy_list, muon_list),reverse=True)[:5]
        ctr=0
        for energy,muon in sorted_energy_particles:
            ctr+=1
            energy_name="Muon_Energy_L{}".format(ctr)
            muon_name="Muon_L{}".format(ctr)
            muon_depth_name="Muon_L{}_Depth".format(ctr)
            frame[energy_name]=dataclasses.I3Double(energy)
            frame[muon_name]=dataclasses.I3Particle(muon)
            
            intersection=surface_det.intersection(muon.pos, muon.dir)#points of intersection
            z_inter=muon.pos.z-intersection.first*np.cos(muon.dir.zenith)
            depth=1948.07-z_inter
            frame[muon_depth_name]=dataclasses.I3Double(depth)
            muon_sep_neut_name="Muon_L{}_Neutrino_separation".format(ctr)
            muon_sep_shower_name="Muon_L{}_Shower_separation".format(ctr)
   
            sep1=get_lateral_separation(muon,frame["ShowerNeutrino"])
            frame[muon_sep_neut_name]=dataclasses.I3Double(sep1)
            sep2=get_lateral_separation(muon,frame["PolyplopiaPrimary"])
            frame[muon_sep_shower_name]=dataclasses.I3Double(sep2)
            
tray.Add(get_deposit_energy)


###########book keys for hdf file###########
table_keys=[
            'GaisserH4a_weight',
            "HomogenizedQTot",
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
            'NeutrinoParent',
            'num_neutrinos'
            #'I3MCTree'
           ]

tray.AddModule('I3Writer', 'writer',
    DropOrphanStreams=[icetray.I3Frame.DAQ],
    Streams=[  icetray.I3Frame.DAQ, icetray.I3Frame.Physics],
    filename=outfile+'.i3.zst')

tray.Add(I3HDFWriter,'HDFwriter',
              Output=outfile+".hdf5",
              keys         = table_keys,
              SubEventStreams = ['InIceSplit']
)
tray.AddModule('TrashCan', 'thecan')

tray.Execute()
tray.Finish()

del tray

stop_time = time.asctime()

print ('Started:', start_time)
print ('Ended:', stop_time)

