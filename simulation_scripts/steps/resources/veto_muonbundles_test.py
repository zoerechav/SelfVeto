from __future__ import division
import numpy as np
np.set_printoptions(formatter={'all':lambda x: str(x)})
import click
import os
from I3Tray import I3Tray
from icecube import icetray, dataclasses,MuonGun, dataio, simclasses, phys_services
from icecube.phys_services import I3Calculator, ExtrudedPolygon
import sys, math
from icecube.dataclasses import I3Double, I3Particle, I3Direction, I3Position, I3VectorI3Particle, I3Constants, I3VectorOMKey
from icecube.simclasses import I3MMCTrack
import collections
from icecube.hdfwriter import I3HDFTableService, I3HDFWriter
import time
import matplotlib.pyplot as plt
from ic3_labels.labels.utils import muon as mu_utils
from ic3_labels.labels.utils import detector, geometry
import matplotlib.path as mpltPath
from scipy.spatial import ConvexHull
import yaml
import traceback
import pickle


with open('/home/zrechav/SelfVeto_Correlation_Tables/SelfVeto_Correlation_Tables/scripts/config.yaml', 'r') as yaml_file:
    config = yaml.safe_load(yaml_file)
    
gcd = config['gcd']

angles_space = eval(config['angles'])
nue_angles = eval(config['nue_angles'])
numu_angles = eval(config['numu_angles'])
depths = eval(config['depths'])
flav_dict={}
flav_dict=    {12:'NuE',-12:'NuE',14:'NuMu',-14:'NuMu',16:'NuTau',-16:'NuTau'}


#deep=2.06
#######FITTING#######
#Parameters:
#    - x: Input data.
#    - A: Amplitude.
#    - mu: Mean in log-space.
#    - sigma1: Spread below the peak.
#    - sigma2: Spread above the peak.
#    - p: Exponent for Gaussian terms.
#    - q: Exponent for exponential modifier on the right side.

def two_sided_log_gaussian(x, A, mu, sigma1, sigma2):
    log_x = np.log(x)
    return np.where(
        log_x < mu,
        A * np.exp(-0.5 * ((log_x - mu) / sigma1) ** 2),
        A * np.exp(-0.5 * ((log_x - mu) / sigma2) ** 2)
    )

def modified_two_sided_log_gaussian(x, A, mu, sigma1, sigma2, p):
    log_x = np.log(x)
    return np.where(
        log_x < mu,
        A * np.exp(-0.5 * (np.abs((log_x - mu) / sigma1) ** p)),
        A * np.exp(-0.5 * (np.abs((log_x - mu) / sigma2) ** p))
    )

def modified_two_sided_log_gaussian_exp(x, A, mu, sigma1, sigma2, p, q):
    """
    Modified two-sided log Gaussian with an exponential modifier on the right side.

    """
    log_x = np.log(x)
    left_side = A * np.exp(-0.5 * (np.abs((log_x - mu) / sigma1) ** p))  # Left side Gaussian
    right_side = A * np.exp(-0.5 * (np.abs((log_x - mu) / sigma2) ** p)) * np.exp(-q * (log_x - mu))  # Right side Gaussian + exponential
    return np.where(log_x < mu, left_side, right_side)

def get_info_backup(pickle):
    bin_centers = None
    params = None
    fit = None
    try:
        best_function = pickle.get('best_function')  # Correct access

        if best_function == 'mge':
            bin_centers = pickle['mge']['bin_centers']
            params = pickle['mge']['params']
            fit = modified_two_sided_log_gaussian_exp(bin_centers, *params)
        elif best_function == 'mg':
            bin_centers = pickle['mg']['bin_centers']
            params = pickle['mg']['params']
            fit = modified_two_sided_log_gaussian(bin_centers, *params)
        elif best_function == 'g':
            bin_centers = pickle['g']['bin_centers']
            params = pickle['g']['params']
            fit = two_sided_log_gaussian(bin_centers, *params)
    except Exception as e:  # Change 'error' to 'Exception'
        print(e)
    return bin_centers, fit, best_function

def get_info(pickle):
    bin_centers = None
    params = None
    fit = None
    try:
        best_function = pickle['best_function']['function']  # Use consistent access

        if best_function == 'mge':
            bin_centers = pickle['mge']['bin_centers']
            params = pickle['mge']['params']
            fit = modified_two_sided_log_gaussian_exp(bin_centers, *params)
        elif best_function == 'mg':
            bin_centers = pickle['mg']['bin_centers']
            params = pickle['mg']['params']
            fit = modified_two_sided_log_gaussian(bin_centers, *params)
        elif best_function == 'g':
            bin_centers = pickle['g']['bin_centers']
            params = pickle['g']['params']
            fit = two_sided_log_gaussian(bin_centers, *params)
    except Exception as e:  # Change 'error' to 'Exception'
        print(e)
        return get_info_backup(pickle)
    return bin_centers, fit, best_function

def convert_pos_list(points):
    """Convert an array like list of points to list of I3Position

    Parameters
    ----------
    points : array_like
        Points (x, y, z) given as array .
        Shape: [n_points, 3]

    Returns
    -------
    list[I3Position]
        A list of I3Position.
    """
    return [dataclasses.I3Position(*point) for point in points]

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
    cyl_bot = -700.
    inlimit = False  
    print(particle1.pos.z)
    if (((particle1.pos.z <=cyl_top) and (particle1.pos.z>=cyl_bot))) and bound_2D.contains_points([(particle1.pos.x, particle1.pos.y)]):
            inlimit=True            
    return inlimit

bound_2D,surface_det = get_surface_det()

def default_sampler( e_min=1e1, e_max=1e7, gamma=3,multiplicity=1):
        """Sample from Powerlaw Distribution

        Sample `num` events from a power law with index gamma between e_min
        and e_max by using the analytic inversion method.
        The power law pdf is given by
        .. math::
           \mathrm{pdf}(\gamma) = x^{-\gamma} / \mathrm{norm}
        where norm ensures an area under curve of one. Positive spectral index
        gamma means a falling spectrum.
        Note: When :math:`\gamma=1` the integral is
        .. math::
           \int 1/x \mathrm{d}x = ln(x) + c
        This case is also handled.

        Parameters
        ----------
        e_min : float
            The lower bound of the PDF, needed for proper normalization.
        e_max : float
            The upper bound of the PDF, needed for proper normalization.
        gamma : float, optional
            Power law index.
        num : int
            Number of random numbers to generate.

        Returns
        -------
        float
            The sampled energy from the specified powerlaw distribution.
        """
        u = np.random.uniform(0., 1.,multiplicity)

        if gamma == 1:
            return np.exp(u * np.log(e_max / e_min)) * e_min
        else:
            radicant = (u * (e_max**(1. - gamma) - e_min**(1. - gamma))
                        + e_min**(1. - gamma))
            return radicant**(1. / (1. - gamma))
        
class InjectVetoMuons(icetray.I3ConditionalModule):

    """Class to inject an accompanying muon for a provided neutrino event

    This is intended to run after neutrino injection. It is expected that the
    first primary in the provided I3MCTree is the injected neutrino event.
    For this particle, an accompanying muon will be injected at the convex hull
    of the detector.

    Attributes
    ----------
    mctree_name : str
        The name of the I3MCTree key.
    """

    def __init__(self, context):
        """Class to inject accompanying muons

        Parameters
        ----------
        context : TYPE
            Description
        """
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddOutBox('OutBox')
        self.AddParameter(
            'mctree_name', 'The name of the I3MCTree to use.', 'I3MCTree')
        self.AddParameter(
            'n_frames_per_neutrino',
            'The number of coincidence events to inject per imported event.',
            1)
        self.AddParameter(
            'sampling_settings',
            'Settings specifying how the muon energy is sampled.',
            {'method': 'power_law', 'range': [10, 1e7], 'gamma': 2})
        self.AddParameter(
            'uncorrelated_muon_settings',
            'If provided, the injected muon will be an uncorrelated coincident'
            ' muon sampled via the provided settings. The provided keys must '
            'include: anchor_(x/y/z)_range, time_range, '
            'azimuth_range [deg], zenith_range [deg]',
            None)
        self.AddParameter(
            'random_service',
            'The random service or seed to use. If this is an '
            'integer, a numpy random state will be created with '
            'the seed set to `random_service`',
            42)
        self.AddParameter(
            'output_key',
            'The key to which the sampling information is written to.',
            'MCVetoMuonInjectionInfo')

    def Configure(self):
        """Configures MuonLossProfileFilter.
        """
        self.mctree_name = self.GetParameter('mctree_name')
        self.n_frames_per_neutrino = self.GetParameter('n_frames_per_neutrino')
        self.sampling_settings = self.GetParameter('sampling_settings')
        self.uncorrelated_muon_settings = self.GetParameter(
            'uncorrelated_muon_settings')
        self.random_service = self.GetParameter('random_service')
        self.sampling_method = self.sampling_settings['method']
        self.output_key = self.GetParameter('output_key')

        if isinstance(self.random_service, int):
            self.random_service = np.random.RandomState(self.random_service)

    def DAQ(self, frame):
        """Inject accompanying muons
        """
        print(self.n_frames_per_neutrino)
        for i in range(self.n_frames_per_neutrino):
            print()
            print(frame['OldRunNumber'].value,';',frame['OldEventID'].value)
            # get I3MCTree
            mc_tree = dataclasses.I3MCTree(frame[self.mctree_name])

            # primary particle
            primary = mc_tree.get_primaries()[0] ###BUG for Zelong

            # inject uncorrelated muon with arbitrary anchor point
            if self.uncorrelated_muon_settings is not None:
                primary = dataclasses.I3Particle(primary)
                x = self.random_service.uniform(
                    *self.uncorrelated_muon_settings['anchor_x_range'])
                y = self.random_service.uniform(
                    *self.uncorrelated_muon_settings['anchor_y_range'])
                z = self.random_service.uniform(
                    *self.uncorrelated_muon_settings['anchor_z_range'])
                t = self.random_service.uniform(
                    *self.uncorrelated_muon_settings['time_range'])
                azimuth = self.random_service.uniform(*np.deg2rad(
                    self.uncorrelated_muon_settings['azimuth_range']))
                zenith = self.random_service.uniform(*np.deg2rad(
                    self.uncorrelated_muon_settings['zenith_range']))

                primary.pos = dataclasses.I3Position(x, y, z)
                primary.dir = dataclasses.I3Direction(zenith, azimuth)
                primary.time = t
            
            # compute entry point
            bound_2D,surface_det = get_surface_det(gcdFile=gcd)
            
            ##define dummy_hull to determine if an intersection exists or not. my_hull are xyz coordinates of cylinder with radius 700m 
            min_z = -700.
            max_z = 550.
            my_hull = [
                [-156.33, -682.32, min_z],
                [67.87, -696.70,   min_z],
                [292.74, -635.85,  min_z],
                [458.33, -529.09,  min_z],
                [589.68, -377.19,  min_z],
                [679.40, -168.55,  min_z],
                [699.65, 22.09,    min_z],
                [675.64, 183.06,   min_z],
                [593.65, 370.92,   min_z],
                [469.41, 519.28,   min_z],
                [338.72, 612.59,   min_z],
                [127.48, 688.29,   min_z],
                [-39.05, 698.91,   min_z],
                [-253.58, 652.45,  min_z],
                [-423.24, 557.56,  min_z],
                [-558.08, 422.55,  min_z],
                [-663.48, 223.16,  min_z],
                [-699.97, 6.06,    min_z],
                [-678.67, -171.48, min_z],
                [-608.19, -346.57, min_z],
                [-480.95, -508.61, min_z],
                [-324.18, -620.41, min_z],
    
                [-156.33, -682.32, max_z],
                [67.87, -696.70,   max_z],
                [292.74, -635.85,  max_z],
                [458.33, -529.09,  max_z],
                [589.68, -377.19,  max_z],
                [679.40, -168.55,  max_z],
                [699.65, 22.09,    max_z],
                [675.64, 183.06,   max_z],
                [593.65, 370.92,   max_z],
                [469.41, 519.28,   max_z],
                [338.72, 612.59,   max_z],
                [127.48, 688.29,   max_z],
                [-39.05, 698.91,   max_z],
                [-253.58, 652.45,  max_z],
                [-423.24, 557.56,  max_z],
                [-558.08, 422.55,  max_z],
                [-663.48, 223.16,  max_z],
                [-699.97, 6.06,    max_z],
                [-678.67, -171.48, max_z],
                [-608.19, -346.57, max_z],
                [-480.95, -508.61, max_z],
                [-324.18, -620.41, max_z],
    
                [-156.33, -682.32, max_z]
    
                    ]
            the_hull = ExtrudedPolygon(convert_pos_list(my_hull))                               
            dummy_intersections = mu_utils.get_muon_convex_hull_intersections(primary, convex_hull=the_hull)

            
            intersections = surface_det.intersection(primary.pos,primary.dir)
            
            if len(dummy_intersections) == 0:
                # particle did not hit convex hull, use closest approach
                closest_position = I3Calculator.closest_approach_position(
                    primary, dataclasses.I3Position(0., 0., 0.))
                distance = mu_utils.get_distance_along_track_to_point(
                    primary.pos, primary.dir, closest_position)
                print('no intersection')
            else:
                print('intersection')
                distance = intersections.first
            print("my_hull calculation: ", distance)
            inj_pos = primary.pos + distance * primary.dir
            print(inj_pos)
            inj_time = primary.time + distance / dataclasses.I3Constants.c
            inj_dir = dataclasses.I3Direction(primary.dir)
            global flavour
            flavour=flav_dict[primary.pdg_encoding]
            print('measuring multiplicity')
            try:
                bundle_multiplicity,depth,binnedsampling=self._sample_multiplicity(primary)
                print('multiplicity:{}, depth:{}, binnedsample:{}'.format(bundle_multiplicity,depth,binnedsampling))
                click.echo('multiplicity measured, going to muon injection')
                            # inject new muon
                print('NeutrinoEnergy',primary.energy)
                
                mc_tree, injection_info,energysample = self._create_bundle_mc_tree(
                    inj_energy=primary.energy,
                    inj_pos=inj_pos,
                    inj_time=inj_time,
                    inj_dir=inj_dir,
                    depth=depth,
                    multiplicity=bundle_multiplicity,
                )
            except Exception as e:
                print(e)
                print(traceback.print_exc())
                

            # copy frame
            frame_copy = icetray.I3Frame(frame)

            # replace I3MCTree
            frame_copy[self.mctree_name + 'VetoMuon_preMuonProp'] = mc_tree

            # add info to frame
            try:
                injection_info['injection_counter'] = float(i)
                if binnedsampling:
                    injection_info['MultiplicitySample'] = True
                else:
                    injection_info['MultiplicitySample'] = False
                if energysample:
                    injection_info['EnergySample'] = True
                else:
                    injection_info['EnergySample'] = False
                frame_copy[self.output_key] = injection_info

            # push frame on to subsequent modules
                self.PushFrame(frame_copy)
            except Exception as e:
                print(e)
            

    
    def _powerlaw_sampler(self, e_min, e_max, gamma):
        """Sample from Powerlaw Distribution

        Sample `num` events from a power law with index gamma between e_min
        and e_max by using the analytic inversion method.
        The power law pdf is given by
        .. math::
           \mathrm{pdf}(\gamma) = x^{-\gamma} / \mathrm{norm}
        where norm ensures an area under curve of one. Positive spectral index
        gamma means a falling spectrum.
        Note: When :math:`\gamma=1` the integral is
        .. math::
           \int 1/x \mathrm{d}x = ln(x) + c
        This case is also handled.

        Parameters
        ----------
        e_min : float
            The lower bound of the PDF, needed for proper normalization.
        e_max : float
            The upper bound of the PDF, needed for proper normalization.
        gamma : float, optional
            Power law index.
        num : int
            Number of random numbers to generate.

        Returns
        -------
        float
            The sampled energy from the specified powerlaw distribution.
        """
        u = self.random_service.uniform(0., 1.)

        if gamma == 1:
            return np.exp(u * np.log(e_max / e_min)) * e_min
        else:
            radicant = (u * (e_max**(1. - gamma) - e_min**(1. - gamma))
                        + e_min**(1. - gamma))
            return radicant**(1. / (1. - gamma))

    def _sample_multiplicity(self,primary):
        """Sample Multiplicity of the Muon bundle for given neutrino

        Returns
        -------
        int
            The bundle multiplicity
        """
        from random import seed
        import random
        mult_bins=eval(config['mult_bins'])
        mult_bin_centers = np.sqrt(mult_bins[:-1] * mult_bins[1:])
        bound_2D,surface_det = get_surface_det(gcdFile=gcd)
        
        Hist2D = None


        intersection = surface_det.intersection(primary.pos,primary.dir)
        print("surface_calculation: ", (intersection.first))
        
        z_inter=primary.pos.z-intersection.first*np.cos(primary.dir.zenith)
        depth=1948.07-z_inter
        print('depth ',depth)
        
        depth=np.around(depths[np.digitize(depth/1000,depths)-1],decimals=2)
        coszen = None
        
        
        flavour == 'NuE' #### Hardcoded for TEST ######
        
        
        
        if flavour == 'NuMu':
            coszen=np.around(angles_space[np.digitize(np.cos(primary.dir.zenith),numu_angles)-1],
                             decimals=2)
            if coszen == 1.0:
                coszen = 0.8
        if flavour == 'NuE':
            coszen=np.around(angles_space[np.digitize(np.cos(primary.dir.zenith),nue_angles)-1],
                             decimals=2)
            if coszen == 1.0:
                coszen = 0.85
        #coszen=np.around(angles_space[np.digitize(np.cos(primary.dir.zenith),angles_space)-1],decimals=2)
        print('coszen ',coszen)
        E_nu_bins=eval(config['E_nu_bins'])
        print('E_nu_bins ', E_nu_bins)
        print('depth in multiplicity',depth)
        
        
        binnedsample=False
        
        try:
            print('MultFile: ',
                      '/data/user/zrechav/output_SelfVeto_Correlation_Tables/multiplicity_mapping/'
                      +flavour+'_Mult_Zen_'+str(coszen)+'_Depth_'+str(depth)+'.npy')
            #sample multiplicity from pdf for specific flavour of neutrinos 
            if os.path.exists('/data/user/zrechav/output_SelfVeto_Correlation_Tables/multiplicity_mapping/'
                              +flavour+'_Mult_Zen_'+str(coszen)+'_Depth_'+str(depth)+'.npy'):
                Hist2D=np.load('/data/user/zrechav/output_SelfVeto_Correlation_Tables/multiplicity_mapping/'
                               +flavour+'_Mult_Zen_'+str(coszen)+'_Depth_'+str(depth)+'.npy')
                #Hist2D = np.nan_to_num(Hist2D)
                print('I the multiplicity file, exist')
            
            Hist2D = np.nan_to_num(Hist2D)
            energy_index=np.digitize(primary.energy,E_nu_bins)
            if energy_index == 6: ##no multiplicity statistics for highest energy bin, assume second high
                energy_index = 5
            print('mult energy index ', primary.energy, ' ',energy_index)
            if (energy_index) > len(E_nu_bins):
                print('length of energy index was larger than length of E_nu_bins')
                #continue
            
            print('HIST2D ', Hist2D)
            Hist1D = Hist2D[energy_index,:]
            print('HIST1D ',Hist1D)
            Hist1D=Hist1D/np.sum(Hist1D)
            print(Hist1D)
            print('length mult_bin centers ',len(mult_bin_centers))
            print('length Hist1D ', len(Hist1D))
            multiplicity= np.random.choice(np.arange(len(mult_bin_centers)), 1, p=Hist1D)[0]
            print('len multiplicity ', multiplicity)
            if multiplicity<1:
                multiplicity=1
            binnedsample=True
        except Exception as e:
            print(e)
            multiplicity=1
        print('multiplicity:{}, depth:{}, binnedsample:{}'.format(multiplicity,depth,binnedsample))
        return multiplicity,depth,binnedsample
        
    def _sample_bundle_energies(self,multiplicity,neutrino_energy,coszen,depth):
        """Sample Energy of the injected Muon

        Returns
        -------
        double
            The sampled energy
        """
        flavor = None
        
        flavour = 'NuE' #####Hardcode for TEST ######
        
        
        
        if flavour == 'NuMu':
            flavor = 'numu'
        if flavour == 'NuE':
            flavor = 'nue'
        E_nu_bins=eval(config['E_nu_bins'])
        energy_index=np.digitize(neutrino_energy,E_nu_bins)-1
        print('og energy index: ', energy_index)
        
        energy_index += 1 #off by one error
        print('new energy index: ', energy_index)
        depth_index=np.digitize(depth,depths)-1
        depth_str = None
        coszen_index=np.digitize(coszen,angles_space)
        print(coszen)
        inj_muon_energies=[]
        if coszen<=0:
            print('CosZen <0')
            return inj_muon_energies,False
        coszen=angles_space[np.digitize(coszen,angles_space)]
        depth=depths[np.digitize(depth,depths)-1]
        coszen=np.around(coszen,decimals=2)
        depth=np.around(depth,decimals=2)
        print('energy ', neutrino_energy)
        print('energy index ', energy_index)
        print('depth ', depth, depth_index)
        print('coszen ',coszen, coszen_index)
        
        if depth_index == 0:
            depth_str = 'top'
        if depth_index == 1:
            depth_str = 'middle'
        if depth_index == 2:
            depth_str = 'middle'
        if depth_index == 3:
            depth_str = 'bottom'
        if depth_index == 4:
            depth_str = 'bottom'
        print(depth_str)
        energysample=True
        print('Bundle Energy Sampler')
        try:
            if multiplicity==1:
                filename=f'/data/user/zrechav/output_SelfVeto_Correlation_Tables/binning/fits/{flavor}/m_one/muon1/{flavor}_one_{depth_str}_nubin{energy_index}_zenbin{coszen_index}.pkl'
                print(filename)
                
                with open(filename,'rb') as f:
                    histogram = pickle.load(f)
                bins,prob,function = get_info(histogram)
                print('bins:',bins)
                print('prob:',prob)
                prob_norm = prob/np.sum(prob) ##decimal cleaning
                sample_energy = (np.random.choice(bins,p=prob_norm))
                ##optional, to know the probability of the sampled energy
                index = np.argmin(np.abs(bins-sample_energy))
                sample_prob = prob_norm[index]
                inj_muon_energies=[sample_energy]
                
            elif multiplicity==2:
                filename1=f'/data/user/zrechav/output_SelfVeto_Correlation_Tables/binning/fits/{flavor}/m_two/muon1/{flavor}_two_{depth_str}_nubin{energy_index}_zenbin{coszen_index}.pkl'
                print(filename1)
                
                with open(filename1,'rb') as f:
                    histogram1 = pickle.load(f)
                bins1,prob1,function1 = get_info(histogram1)
                print('bins:',bins1)
                print('prob:',prob1)
                prob_norm1 = prob1/np.sum(prob1) ##decimal cleaning
                sample_energy1 = (np.random.choice(bins1,p=prob_norm1))
                ##optional, to know the probability of the sampled energy
                index1 = np.argmin(np.abs(bins1-sample_energy1))
                sample_prob1 = prob_norm1[index1]
                #inj_muon_energies=[sample_energy]
                filename2=f'/data/user/zrechav/output_SelfVeto_Correlation_Tables/binning/fits/{flavor}/m_two/muon2/{flavor}_two_{depth_str}_nubin{energy_index}_zenbin{coszen_index}.pkl'
                print(filename2)
                
                with open(filename2,'rb') as f:
                    histogram2 = pickle.load(f)
                bins2,prob2,function2 = get_info(histogram2)
                print('bins:',bins2)
                print('prob:',prob2)
                prob_norm2 = prob2/np.sum(prob2) ##decimal cleaning
                sample_energy2 = (np.random.choice(bins2,p=prob_norm2))
                ##optional, to know the probability of the sampled energy
                index2 = np.argmin(np.abs(bins2-sample_energy2))
                sample_prob2 = prob_norm2[index2]
                
                inj_muon_energies=[sample_energy1,sample_energy2]

            elif multiplicity==3:
                filename1=f'/data/user/zrechav/output_SelfVeto_Correlation_Tables/binning/fits/{flavor}/m_three/muon1/{flavor}_three_{depth_str}_nubin{energy_index}_zenbin{coszen_index}.pkl'
                print(filename1)

                with open(filename1,'rb') as f:
                    histogram1 = pickle.load(f)
                bins1,prob1,function1 = get_info(histogram1)
                print('bins:',bins1)
                print('prob:',prob1)
                prob_norm1 = prob1/np.sum(prob1) ##decimal cleaning
                sample_energy1 = (np.random.choice(bins1,p=prob_norm1))
                ##optional, to know the probability of the sampled energy
                index1 = np.argmin(np.abs(bins1-sample_energy1))
                sample_prob1 = prob_norm1[index1]
                #inj_muon_energies=[sample_energy]
                filename2=f'/data/user/zrechav/output_SelfVeto_Correlation_Tables/binning/fits/{flavor}/m_three/muon2/{flavor}_three_{depth_str}_nubin{energy_index}_zenbin{coszen_index}.pkl'
                print(filename2)
                
                with open(filename2,'rb') as f:
                    histogram2 = pickle.load(f)
                bins2,prob2,function2 = get_info(histogram2)
                print('bins:',bins2)
                print('prob:',prob2)
                prob_norm2 = prob2/np.sum(prob2) ##decimal cleaning
                sample_energy2 = (np.random.choice(bins2,p=prob_norm2))
                ##optional, to know the probability of the sampled energy
                index2 = np.argmin(np.abs(bins2-sample_energy2))
                sample_prob2 = prob_norm2[index2]
                #inj_muon_energies=[sample_energy]
                filename3=f'/data/user/zrechav/output_SelfVeto_Correlation_Tables/binning/fits/{flavor}/m_three/muon3/{flavor}_three_{depth_str}_nubin{energy_index}_zenbin{coszen_index}.pkl'
                print(filename3)
                
                with open(filename3,'rb') as f:
                    histogram3 = pickle.load(f)
                bins3,prob3,function3 = get_info(histogram3)
                print('bins:',bins3)
                print('prob:',prob3)
                prob_norm3 = prob3/np.sum(prob3) ##decimal cleaning
                sample_energy3 = (np.random.choice(bins3,p=prob_norm3))
                ##optional, to know the probability of the sampled energy
                index3 = np.argmin(np.abs(bins3-sample_energy3))
                sample_prob3 = prob_norm3[index3]
                
                
                inj_muon_energies=[sample_energy1,sample_energy2,sample_energy3]
               
        

            elif multiplicity>=4:
                filename1=f'/data/user/zrechav/output_SelfVeto_Correlation_Tables/binning/fits/{flavor}/m_four/muon1/{flavor}_four_{depth_str}_nubin{energy_index}_zenbin{coszen_index}.pkl'
                print(filename1)

                with open(filename1,'rb') as f:
                    histogram1 = pickle.load(f)
                bins1,prob1,function1 = get_info(histogram1)
                print('bins:',bins1)
                print('prob:',prob1)
                prob_norm1 = prob1/np.sum(prob1) ##decimal cleaning
                sample_energy1 = (np.random.choice(bins1,p=prob_norm1))
                ##optional, to know the probability of the sampled energy
                index1 = np.argmin(np.abs(bins1-sample_energy1))
                sample_prob1 = prob_norm1[index1]
                #inj_muon_energies=[sample_energy]
                filename2=f'/data/user/zrechav/output_SelfVeto_Correlation_Tables/binning/fits/{flavor}/m_four/muon2/{flavor}_four_{depth_str}_nubin{energy_index}_zenbin{coszen_index}.pkl'
                print(filename2)
                
                with open(filename2,'rb') as f:
                    histogram2 = pickle.load(f)
                bins2,prob2,function2 = get_info(histogram2)
                print('bins:',bins2)
                print('prob:',prob2)
                prob_norm2 = prob2/np.sum(prob2) ##decimal cleaning
                sample_energy2 = (np.random.choice(bins2,p=prob_norm2))
                ##optional, to know the probability of the sampled energy
                index2 = np.argmin(np.abs(bins2-sample_energy2))
                sample_prob2 = prob_norm2[index2]
                #inj_muon_energies=[sample_energy]
                filename3=f'/data/user/zrechav/output_SelfVeto_Correlation_Tables/binning/fits/{flavor}/m_four/muon3/{flavor}_four_{depth_str}_nubin{energy_index}_zenbin{coszen_index}.pkl'
                print(filename3)
                
                with open(filename3,'rb') as f:
                    histogram3 = pickle.load(f)
                bins3,prob3,function3 = get_info(histogram3)
                print('bins:',bins3)
                print('prob:',prob3)
                prob_norm3 = prob3/np.sum(prob3) ##decimal cleaning
                sample_energy3 = (np.random.choice(bins3,p=prob_norm3))
                ##optional, to know the probability of the sampled energy
                index3 = np.argmin(np.abs(bins3-sample_energy3))
                sample_prob3 = prob_norm3[index3]
                filename4=f'/data/user/zrechav/output_SelfVeto_Correlation_Tables/binning/fits/{flavor}/m_four/muon4/{flavor}_four_{depth_str}_nubin{energy_index}_zenbin{coszen_index}.pkl'
                print(filename4)
                
                with open(filename4,'rb') as f:
                    histogram4 = pickle.load(f)
                bins4,prob4,function4 = get_info(histogram4)
                print('bins:',bins4)
                print('prob:',prob4)
                prob_norm4 = prob4/np.sum(prob4) ##decimal cleaning
                sample_energy4 = (np.random.choice(bins4,p=prob_norm4))
                ##optional, to know the probability of the sampled energy
                index4 = np.argmin(np.abs(bins4-sample_energy4))
                sample_prob4 = prob_norm4[index4]                
                
                
                inj_muon_energies=[sample_energy1,sample_energy2,sample_energy3,sample_energy4]
               

        except Exception as e:
            print(e)
            print(traceback.print_exc())
            print('THE DEFAULT SAMPLER IS BEING USED!')
            inj_muon_energies=default_sampler(multiplicity=multiplicity)
            energysample=False

            
        print('EnergySample',energysample)
        return inj_muon_energies,energysample

#     def _sample_bundle_energies(self, multiplicity, primary_nu_energy, coszen, depth):
#         """Sample Energy of the injected Muon bundle with rejection sampling.
#         """
#         flavor = None
        
#         flavour == 'NuE' #### Hardcoded for TEST ######
        
        
        
#         if flavour == 'NuMu':
#             flavor = 'numu'
#         if flavour == 'NuE':
#             flavor = 'nue'

#         E_nu_bins = eval(config['E_nu_bins'])
#         energy_index = np.digitize(primary_nu_energy, E_nu_bins) - 1
#         energy_index += 1  # off-by-one correction

#         depth_index = np.digitize(depth, depths) - 1
#         coszen_index = np.digitize(coszen, angles_space)

#         inj_muon_energies = []

#         if coszen <= 0:
            
#             bundle_diagnostics = {
#                 #"bundle_n_tries": 0.0,
#                 "bundle_energy": 0.0,
#                 "bundle_primary": float(primary_nu_energy),
#                 "bundle_coszen": float(coszen),
#                 "bundle_depth": float(depth),
#             }
#             return inj_muon_energies, False, bundle_diagnostics

#         coszen = angles_space[np.digitize(coszen, angles_space)]
#         depth = depths[np.digitize(depth, depths) - 1]
#         coszen = np.around(coszen, decimals=2)
#         depth = np.around(depth, decimals=2)

#         if depth_index == 0:
#             depth_str = 'top'
#         elif depth_index in [1, 2]:
#             depth_str = 'middle'
#         else:
#             depth_str = 'bottom'

#         energysample = True

#         #MAX_TRIES = 1
#         #n_tries = 0


#         try:

#             inj_muon_energies = []


#             if multiplicity == 1:
#                 filename = (
#                     f'/data/user/zrechav/output_SelfVeto_Correlation_Tables/binning/fits/{flavor}/m_one/muon1/{flavor}_one_{depth_str}_nubin{energy_index}_zenbin{coszen_index}.pkl'
#                 )
#                 with open(filename, 'rb') as f:
#                     histogram = pickle.load(f)

#                 bins, prob, _ = get_info(histogram)
#                 prob_norm = prob / np.sum(prob)
#                 sample_energy = np.random.choice(bins, p=prob_norm)
#                 inj_muon_energies = [sample_energy]

#             elif multiplicity == 2:
#                 energies = []
#                 for i in [1, 2]:
#                     filename = (
#                         f'/data/user/zrechav/output_SelfVeto_Correlation_Tables/binning/fits/{flavor}/m_two/muon{i}/{flavor}_two_{depth_str}_nubin{energy_index}_zenbin{coszen_index}.pkl'
#                     )
#                     with open(filename, 'rb') as f:
#                         histogram = pickle.load(f)

#                     bins, prob, _ = get_info(histogram)
#                     prob_norm = prob / np.sum(prob)
#                     energies.append(np.random.choice(bins, p=prob_norm))

#                 inj_muon_energies = energies

#             elif multiplicity == 3:
#                 energies = []
#                 for i in [1, 2, 3]:
#                     filename = (
#                         f'/data/user/zrechav/output_SelfVeto_Correlation_Tables/binning/fits/{flavor}/m_three/muon{i}/{flavor}_three_{depth_str}_nubin{energy_index}_zenbin{coszen_index}.pkl'
#                     )
#                     with open(filename, 'rb') as f:
#                         histogram = pickle.load(f)

#                     bins, prob, _ = get_info(histogram)
#                     prob_norm = prob / np.sum(prob)
#                     energies.append(np.random.choice(bins, p=prob_norm))

#                 inj_muon_energies = energies

#             elif multiplicity >= 4:
#                 energies = []
#                 for i in [1, 2, 3, 4]:
#                     filename = (
#                         f'/data/user/zrechav/output_SelfVeto_Correlation_Tables/binning/fits/{flavor}/m_four/muon{i}/{flavor}_four_{depth_str}_nubin{energy_index}_zenbin{coszen_index}.pkl'
#                     )
#                     with open(filename, 'rb') as f:
#                         histogram = pickle.load(f)

#                     bins, prob, _ = get_info(histogram)
#                     prob_norm = prob / np.sum(prob)
#                     energies.append(np.random.choice(bins, p=prob_norm))

#                 inj_muon_energies = energies

#             total_mu_energy = np.sum(inj_muon_energies)
#             bundle_diagnostics = {
#                     #"bundle_n_tries": float(n_tries),
#                     "bundle_energy": float(total_mu_energy),
#                     "bundle_primary": float(primary_nu_energy),
#                     "bundle_coszen": float(coszen),
#                     "bundle_depth": float(depth),
#                 }





#             #if total_mu_energy <= primary_nu_energy:
#              #   break  # ✅ accept
#             # else: reject and resample

#         except Exception as e:
#             print(e)
#             print(traceback.print_exc())
#             print('THE DEFAULT SAMPLER IS BEING USED!')
#             inj_muon_energies = default_sampler(multiplicity=multiplicity)
#             energysample = False

#         print('EnergySample:', energysample)
#         return inj_muon_energies, energysample, bundle_diagnostics


    def _sample_energy(self):
        """Sample Energy of the injected Muon

        Returns
        -------
        double
            The sampled energy
        """
        if self.sampling_method == 'power_law':
            energy = self._powerlaw_sampler(
                e_min=self.sampling_settings['range'][0],
                e_max=self.sampling_settings['range'][1],
                gamma=self.sampling_settings['gamma'],
            )


        else:
            raise ValueError('Unknown sampling method: {}'.format(
                self.sampling_method))
        return energy
    def _create_mc_tree(self, inj_pos, inj_time, inj_dir):
        """Inject accompanying muon in provided I3MCTree

        Parameters
        ----------
        inj_pos : I3Position
            The position at which to inject the muon.
        inj_time : double
            The time at which to inject the muon.
        inj_dir : I3Direction
            The direction of the injected muon.

        Returns
        -------
        I3MCTree
            The modified I3MCTree with the injected Muon.
        I3MapStringDouble
            Information on the injected muon.
        """
        mc_tree = dataclasses.I3MCTree()
        muon_primary = dataclasses.I3Particle()
        muon_primary.shape = dataclasses.I3Particle.ParticleShape.Primary
        muon_primary.dir = dataclasses.I3Direction(inj_dir)
        muon_primary.pos = dataclasses.I3Position(inj_pos)
        muon_primary.time = inj_time

        muon = dataclasses.I3Particle()
        muon.dir = dataclasses.I3Direction(inj_dir)
        muon.pos = dataclasses.I3Position(inj_pos)
        muon.time = inj_time
        muon.location_type = dataclasses.I3Particle.LocationType.InIce

        # sample type: MuPlus or MuMinus
        u = self.random_service.uniform(0., 1.)
        if u > 0.5:
            pdg_encoding = 13
        else:
            pdg_encoding = -13
        muon.pdg_encoding = pdg_encoding

        # sample energy
        muon.energy = self._sample_energy()

        # add muon primary to I3MCTree
        mc_tree.add_primary(muon_primary)

        # add muon as new child
        mc_tree.append_child(muon_primary, muon)

        # add info
        injection_info = dataclasses.I3MapStringDouble({
            'muon_energy': muon.energy,
            'muon_pdg_encoding': muon.pdg_encoding,
        })

        return mc_tree, injection_info
        
    def _create_bundle_mc_tree(self, inj_energy,inj_pos, inj_time, inj_dir,depth,multiplicity):
        """Inject accompanying muon in provided I3MCTree

        Parameters
        ----------
        inj_energy : Primary energy 
            Neutrino Energy.
        inj_pos : I3Position
            The position at which to inject the muons.
        inj_time : double
            The time at which to inject the muons.
        inj_dir : I3Direction
            The direction of the injected muons.
        multiplicity : int
            The number of injected muons.

        Returns
        -------
        I3MCTree
            The modified I3MCTree with the injected Muon.
        I3MapStringDouble
            Information on the injected muon.
        """
        mc_tree = dataclasses.I3MCTree()
        print('Empty MC Tree created')
        injectdict={}
        keys=['muon_1_energy','muon_1_pdg_encoding','muon_2_energy','muon_2_pdg_encoding','muon_3_energy','muon_3_pdg_encoding','muon_4_energy','muon_4_pdg_encoding']
        for key in keys:#initializing injection dictionary to zero
            injectdict[key]=0
        print('Starting Sampler')
        inj_muon_energies,energysample=self._sample_bundle_energies(multiplicity,inj_energy,np.cos(inj_dir.zenith),depth)

        #for key, val in bundle_diagnostics.items():
        #    injectdict[key] = val

        #len(inj_muon_energies) is the injected multiplicity
        print('Energies',inj_muon_energies)
        ctr=1
        if not inj_muon_energies:
            injection_info=dataclasses.I3MapStringDouble(injectdict)
            return mc_tree, injection_info,energysample
        for muon_energy in inj_muon_energies:

            muon_primary = dataclasses.I3Particle()
            muon_primary.shape = dataclasses.I3Particle.ParticleShape.Primary
            muon_primary.dir = dataclasses.I3Direction(inj_dir)
            muon_primary.pos = dataclasses.I3Position(inj_pos)
            muon_primary.time = inj_time

            muon = dataclasses.I3Particle()
            muon.dir = dataclasses.I3Direction(inj_dir)
            muon.pos = dataclasses.I3Position(inj_pos)
            muon.time = inj_time
            muon.location_type = dataclasses.I3Particle.LocationType.InIce

            # sample type: MuPlus or MuMinus
            u = self.random_service.uniform(0., 1.)
            if u > 0.5:
                pdg_encoding = 13
            else:
                pdg_encoding = -13
            muon.pdg_encoding = pdg_encoding

            # sample energy
            # muon.energy = self._sample_energy()
            
            muon.energy = muon_energy

            # add muon primary to I3MCTree
            mc_tree.add_primary(muon_primary)

            # add muon as new child
            mc_tree.append_child(muon_primary, muon)

            # add info
            keys=['muon_'+str(ctr)+'_energy','muon_'+str(ctr)+'_pdg_encoding']
            ctr+=1
            injectdict[keys[0]]=muon.energy
            injectdict[keys[1]]=muon.pdg_encoding
            # injectdict={
            #     'muon_energy': muon.energy,
            #     'muon_pdg_encoding': muon.pdg_encoding,
            # }
        injection_info=dataclasses.I3MapStringDouble(injectdict)
        return mc_tree, injection_info,energysample


class CombineMCTrees(icetray.I3ConditionalModule):

    def __init__(self, context):
        """Class to combine two I3MCTrees

        Parameters
        ----------
        context : TYPE
            Description
        """
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddOutBox('OutBox')
        self.AddParameter('tree1', 'The name of the first I3MCTree.')
        self.AddParameter('tree2', 'The name of the first I3MCTree.')
        self.AddParameter(
            'output_key', 'Key to which the combined I3MCTree is written to.')

    def Configure(self):
        """Configures MuonLossProfileFilter.
        """
        self.tree1 = self.GetParameter('tree1')
        self.tree2 = self.GetParameter('tree2')
        self.output_key = self.GetParameter('output_key')

    def DAQ(self, frame):
        """Merge trees
        """
        tree1 = dataclasses.I3MCTree(frame[self.tree1])
        tree2 = dataclasses.I3MCTree(frame[self.tree2])
        tree1.merge(tree2)
        frame[self.output_key] = tree1
        self.PushFrame(frame)
