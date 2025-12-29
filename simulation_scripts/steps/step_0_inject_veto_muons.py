#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/icetray-start
#METAPROJECT icetray/v1.9.2
from __future__ import division
import click
import yaml
import numpy as np
import glob

from icecube.simprod import segments

from I3Tray import I3Tray, I3Units
from icecube import icetray, dataclasses

from utils import create_random_services, get_run_folder
from resources import geometry
from resources.import_events import ImportEvents
from resources.biased_simulation import BaseSimulationBias
# from resources.veto_muon import InjectVetoMuons, CombineMCTrees
from resources.veto_muonbundles import InjectVetoMuons, CombineMCTrees
from scipy.spatial import ConvexHull


@click.command()
@click.argument('cfg', type=click.Path(exists=True))
@click.argument('run_number', type=int)
@click.option('--scratch/--no-scratch', default=True)

def main(cfg, run_number, scratch):
    with open(cfg, 'r') as stream:
        cfg = yaml.full_load(stream)
    cfg['run_number'] = run_number
    cfg['run_folder'] = get_run_folder(run_number)
    if scratch:
        outfile = cfg['scratchfile_pattern'].format(**cfg)
    else:
        outfile = cfg['outfile_pattern'].format(**cfg)
    outfile = outfile.replace(' ', '0')

    # ------------------------------
    # get list of files for this run
    # ------------------------------
    import_cfg = cfg['event_import_settings']
    glob_files = import_cfg['input_file_glob_list']

    import_cfg['run_number'] = cfg['run_number']
    import_cfg['dataset_number'] = cfg['dataset_number']
    import_cfg['folder_num_pre_offset'] = cfg['run_number']//1000
    import_cfg['folder_num'] = (
        import_cfg['folder_offset'] + cfg['run_number']//1000
    )
    import_cfg['folder_pattern'] = import_cfg['folder_pattern'].format(
        **import_cfg)
    import_cfg['run_folder'] = import_cfg['folder_pattern'].format(
        **import_cfg)
    glob_files = import_cfg['input_file_glob_list']
    run_number=import_cfg['run_number'] 
    click.echo('Run: {}'.format(run_number))
    click.echo('Outfile: {}'.format(outfile))
    click.echo('Keys to import: {}'.format(import_cfg['keys_to_import']))
    click.echo('Infile: {}'.format(glob_files))
    # single_file=glob_files[0]
    files = [glob_files[0]+"{:03d}".format(run_number)+".i3.zst"] ###modified!!!

    # sort files
    # files = sorted(files)
    # ------------------------------

    
    click.echo('input Files:')
    for file in files:
        click.echo('\t{}'.format(file))
    click.echo('create random services')
    # crate random services
    if 'random_service_use_gslrng' not in cfg:
        click.echo('random_service_use_gslrng False')
        cfg['random_service_use_gslrng'] = False
    random_services, _ = create_random_services(
        dataset_number=cfg['dataset_number'],
        run_number=cfg['run_number'],
        seed=cfg['seed'],
        n_services=3,
        use_gslrng=cfg['random_service_use_gslrng'])

    # --------------------------------------
    # Build IceTray
    # --------------------------------------
    tray = I3Tray()
    # import events from another I3-file
    if 'num_events' not in import_cfg:
        import_cfg['num_events'] = None

    tray.AddModule(
        ImportEvents,
        'ImportEvents',
        files=files,
        num_events=import_cfg['num_events'],
        keys_to_import=import_cfg['keys_to_import'],
        rename_dict=import_cfg['rename_dict'],
        mctree_name=import_cfg['mctree_name'],
    )

    
    if 'mctree_name' not in cfg['veto_muon_injection_config']:
        cfg['veto_muon_injection_config']['mctree_name'] = (
            import_cfg['mctree_name']
        )
    # Level0.0_nugen_IC86.2012_pass2.022010.000000.i3.bz2
    outfile_runnum=outfile.split('.')[-4][1:]

    nfiles=0
    L5_nfiles=0 
    def EventChecker(frame):
        print('checking events?')
        if FileType==0:
            frame_Run=frame['OldRunNumber'].value
            frame_Event=frame['OldEventID'].value
            
            if (frame_Run in filenums) and (frame.Stop == icetray.I3Frame.Physics):
                indices=np.where(filenums == np.int64(frame_Run))
                # print(filenums[indices])
                if frame_Event in eventnums[indices]:
                    print('filenum:',frame['OldRunNumber'].value,' Event:',frame['OldEventID'].value)
                    global nfiles
                    nfiles+=1
                else:
                    print('Event rejected!')
                    return False
            else:
                    return False 
    

    tray.AddModule(EventChecker)
    tray.AddModule(
        InjectVetoMuons, 'InjectVetoMuons',
        random_service=random_services[0],
        **cfg['veto_muon_injection_config']
    )
#     ##In Vedant's code, this cut removes events that were sampled using single power law or default power law --> I shouldn't need this
#     def uncorr_cut(frame):
#         if frame['MCVetoMuonInjectionInfo']['EnergySample']==False:
#             click.echo('Uncorrelated Events')
#             return False
#         else:
#             return True
#     tray.AddModule(uncorr_cut, "Uncorrelated Cut")

    t_name = import_cfg['mctree_name']

    # # rename I3MCTrees so that we can run PROPOSAL
    # def rename_keys(frame, rename_dict):
    #     for key, new_name in rename_dict.items():
    #         if key in frame:
    #             frame[new_name] = frame[key]
    #             del frame[key]
    #         print('key', key)

    # rename_dict = {
    #     t_name: '_' + t_name,
    #     t_name + '_preMuonProp': '_' + t_name + '_preMuonProp',
    #     t_name + '_preMuonProp_RNGState': (
    #         '_' + t_name + '_preMuonProp_RNGState',
    #     ),
    #     t_name + '_veto_muon': 'I3MCTree_preMuonProp',
    # }
    # tray.AddModule(rename_keys, 'TempRename', rename_dict=rename_dict)

    # temporarily save MMCTrackList
    tray.AddModule('Rename', keys=['MMCTrackList', '_MMCTrackList'])
    # rename other keys
    tray.AddModule('Rename', keys=['cc', 'og_cc'])
    tray.AddModule('Rename', keys=['BDT_bdt_max_depth_4_n_est_1000lr_0.01_seed_3_train_size_50', 'og_BDT_bdt_max_depth_4_n_est_1000lr_0.01_seed_3_train_size_50'])
    tray.AddModule('Rename', keys=['BDT_bdt_max_depth_4_n_est_2000lr_0_02_seed_3_train_size_50', 'og_BDT_bdt_max_depth_4_n_est_2000lr_0_02_seed_3_train_size_50'])

    # propagate injected muon
    cfg['muon_propagation_config'].update({
        'InputMCTreeName': t_name + 'VetoMuon_preMuonProp',
        'OutputMCTreeName': t_name + 'VetoMuon',
    })
    tray.AddSegment(segments.PropagateMuons,
                    'propagate_muons',
                    RandomService=random_services[1],
                    **cfg['muon_propagation_config'])

    # undo renaming
    tray.AddModule('Rename', keys=['MMCTrackList', 'MMCTrackListVetoMuon'])
    tray.AddModule('Rename', keys=['_MMCTrackList', 'MMCTrackList'])

    
    # create combined I3MCTree for CLSIM
    tray.AddModule(
        CombineMCTrees, 'CombineMCTrees',
        tree1=t_name,
        tree2=t_name + 'VetoMuon',
        output_key='CombinedMuonVetoI3MCTree',
    )

    # Bias simulation if desired
    if 'ApplyBaseSimulationBias' in cfg and cfg['ApplyBaseSimulationBias']:
        tray.AddModule(
            BaseSimulationBias,
            'BaseSimulationBias',
            random_service=random_services[2],
            **cfg['BaseSimulationBiasSettings']
        )

    # --------------------------------------
    # Distance Splits
    # --------------------------------------
    if cfg['distance_splits'] is not None:
        click.echo('SplittingDistance: {}'.format(
            cfg['distance_splits']))
        distance_splits = np.atleast_1d(cfg['distance_splits'])
        dom_limits = np.atleast_1d(cfg['threshold_doms'])
        if len(dom_limits) == 1:
            dom_limits = np.ones_like(distance_splits) * cfg['threshold_doms']
        oversize_factors = np.atleast_1d(cfg['oversize_factors'])
        order = np.argsort(distance_splits)

        distance_splits = distance_splits[order]
        dom_limits = dom_limits[order]
        oversize_factors = oversize_factors[order]

        stream_objects = generate_stream_object(distance_splits,
                                                dom_limits,
                                                oversize_factors)
        tray.AddModule(OversizeSplitterNSplits,
                       "OversizeSplitterNSplits",
                       thresholds=distance_splits,
                       thresholds_doms=dom_limits,
                       oversize_factors=oversize_factors)
        for stream_i in stream_objects:
            outfile_i = stream_i.transform_filepath(outfile)
            tray.AddModule("I3Writer",
                           "writer_{}".format(stream_i.stream_name),
                           Filename=outfile_i,
                           Streams=[icetray.I3Frame.DAQ,
                                    icetray.I3Frame.Physics,
                                    icetray.I3Frame.Stream('S'),
                                    icetray.I3Frame.Stream('M')],
                           If=stream_i)
            click.echo('Output ({}): {}'.format(stream_i.stream_name,
                                                outfile_i))
    else:
        click.echo('Output: {}'.format(outfile))
        tray.AddModule("I3Writer", "writer",
                       Filename=outfile,
                       Streams=[icetray.I3Frame.DAQ,
                                icetray.I3Frame.Physics,
                                icetray.I3Frame.Stream('S'),
                                icetray.I3Frame.Stream('M')])
    # --------------------------------------

    click.echo('Scratch: {}'.format(scratch))
    tray.AddModule("TrashCan", "the can")
    tray.Execute()
    tray.Finish()


if __name__ == '__main__':
    main()
