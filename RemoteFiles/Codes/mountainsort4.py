from mountainlab_pytools import mdaio
from mountainlab_pytools import mlproc as mlp
import os
import json
import numpy as np


def sort_dataset(dataset_dir, output_dir, freq_min=300, freq_max=6000, adjacency_radius=100, detect_threshold=3,
                 fs=20000, detect_sign=-1, opts={}):
    
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    # Dataset parameters
    # ds_params = read_dataset_params(dataset_dir)

    # Bandpass filter
    bandpass_filter(
        timeseries=dataset_dir + '/raw.mda',
        timeseries_out=output_dir + '/filt.mda.prv',
        samplerate=fs,
        freq_min=freq_min,
        freq_max=freq_max,
        opts=opts
    )

    # Whiten
    whiten(
        timeseries=output_dir + '/filt.mda.prv',
        timeseries_out=output_dir + '/pre.mda.prv',
        opts=opts
    )

    # Sort
    # detect_sign = -1
    # if 'spike_sign' in ds_params:
    #     detect_sign = ds_params['spike_sign']
    # if 'detect_sign' in ds_params:
    #     detect_sign = ds_params['detect_sign']
    ms4alg_sort(
        timeseries=output_dir + '/pre.mda.prv',
        geom=dataset_dir + '/geom.csv',
        firings_out=output_dir + '/firings_uncurated.mda',
        adjacency_radius=adjacency_radius,
        detect_sign=detect_sign,
        detect_threshold=detect_threshold,
        opts=opts
    )

    # Compute cluster metrics
    compute_cluster_metrics(
        timeseries=output_dir + '/pre.mda.prv',
        firings=output_dir + '/firings_uncurated.mda',
        metrics_out=output_dir + '/cluster_metrics.json',
        samplerate=fs,
        opts=opts
    )

    # Automated curation
    automated_curation(
        firings=output_dir + '/firings_uncurated.mda',
        cluster_metrics=output_dir + '/cluster_metrics.json',
        firings_out=output_dir + '/firings.mda',
        opts=opts
    )


def sort_dataset_segments_part1(dataset_list, output_list, combined_dir, freq_min=300, freq_max=6000,
                                adjacency_radius=100, detect_threshold=3, fs=20000, detect_sign=-1, opts={}):

    timeseries_list=[]
    firings_list=[]

    for dataset_dir, output_dir in zip(dataset_list, output_list):

        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

        # Dataset parameters
        # ds_params = read_dataset_params(dataset_dir)

        # Bandpass filter
        bandpass_filter(
            timeseries=dataset_dir + '/raw.mda',
            timeseries_out=output_dir + '/filt.mda.prv',
            samplerate=fs,
            freq_min=freq_min,
            freq_max=freq_max,
            opts=opts
        )

        # Whiten
        whiten(
            timeseries=output_dir + '/filt.mda.prv',
            timeseries_out=output_dir + '/pre.mda.prv',
            opts=opts
        )

        # Sort
        # detect_sign = -1

        ms4alg_sort(
            timeseries=output_dir + '/pre.mda.prv',
            geom=dataset_dir + '/geom.csv',
            firings_out=output_dir + '/firings_uncurated.mda',
            adjacency_radius=adjacency_radius,
            detect_sign=detect_sign,
            detect_threshold=detect_threshold,
            opts=opts
        )

        firings_list.append(output_dir + '/firings_uncurated.mda')
        timeseries_list.append(output_dir + '/pre.mda.prv')

    offset_list = get_time_offset_list(dataset_list)
    offsets = ','.join(map(str, offset_list))

    # Bandpass filter
    bandpass_filter(
        timeseries=combined_dir + '/combined_raw.mda',
        timeseries_out=combined_dir + '/filt.mda.prv',
        samplerate=fs,
        freq_min=freq_min,
        freq_max=freq_max,
        opts=opts
    )

    # Whiten
    whiten(
        timeseries=combined_dir + '/filt.mda.prv',
        timeseries_out=combined_dir + '/pre.mda.prv',
        opts=opts
    )

    return timeseries_list, firings_list, offsets

    # pyms_anneal_segs(
    #     timeseries_list=timeseries_list,
    #     firings_list=firings_list,
    #     firings_out=combined_dir + '/combined_firings_uncurated.mda',
    #     dmatrix_out=[],
    #     k1_dmatrix_out=[],
    #     k2_dmatrix_out=[],
    #     dmatrix_templates_out=[],
    #     time_offsets=offsets
    # )
    #
    # # Compute cluster metrics
    # compute_cluster_metrics(
    #     timeseries=combined_dir + 'pre.mda.prv',
    #     firings=combined_dir + '/combined_firings_uncurated.mda',
    #     metrics_out=combined_dir + '/cluster_metrics.json',
    #     samplerate=ds_params['samplerate'],
    #     opts=opts
    # )
    #
    # # Automated curation
    # automated_curation(
    #     firings=combined_dir + '/combined_firings_uncurated.mda',
    #     cluster_metrics=combined_dir + '/cluster_metrics.json',
    #     firings_out=combined_dir + '/firings.mda',
    #     opts=opts
    # )

def sort_dataset_segments_part2(combined_dir, timeseries_list, firings_list, offsets):


    pyms_anneal_segs(
        timeseries_list=timeseries_list,
        firings_list=firings_list,
        firings_out=combined_dir + '/firings_uncurated.mda',
        dmatrix_out=[],
        k1_dmatrix_out=[],
        k2_dmatrix_out=[],
        dmatrix_templates_out=[],
        time_offsets=offsets
    )


def sort_dataset_segments_part3(combined_dir, fs=20000, opts={}):


    # Compute cluster metrics
    compute_cluster_metrics(
        timeseries=combined_dir + '/pre.mda.prv',
        firings=combined_dir + '/firings_uncurated.mda',
        metrics_out=combined_dir + '/cluster_metrics.json',
        samplerate=fs,
        opts=opts
    )

    # Automated curation
    automated_curation(
        firings=combined_dir + '/firings_uncurated.mda',
        cluster_metrics=combined_dir + '/cluster_metrics.json',
        firings_out=combined_dir + '/firings.mda',
        opts=opts
    )


def read_dataset_params(dsdir):
    params_fname = mlp.realizeFile(dsdir + '/params.json')
    if not os.path.exists(params_fname):
        raise Exception('Dataset parameter file does not exist: ' + params_fname)
    with open(params_fname) as f:
        return json.load(f)


def bandpass_filter(*, timeseries, timeseries_out, samplerate, freq_min, freq_max, opts={}):
    return mlp.addProcess(
        'ephys.bandpass_filter',
        {
            'timeseries': timeseries
        },
        {
            'timeseries_out': timeseries_out
        },
        {
            'samplerate': samplerate,
            'freq_min': freq_min,
            'freq_max': freq_max
        },
        opts
    )


def whiten(*, timeseries, timeseries_out, opts={}):
    return mlp.addProcess(
        'ephys.whiten',
        {
            'timeseries': timeseries
        },
        {
            'timeseries_out': timeseries_out
        },
        {},
        opts
    )


def ms4alg_sort(*, timeseries, geom, firings_out, detect_sign, adjacency_radius, detect_threshold=3, opts={}):
    pp = {}
    pp['detect_sign'] = detect_sign
    pp['adjacency_radius'] = adjacency_radius
    pp['detect_threshold'] = detect_threshold
    mlp.addProcess(
        'ms4alg.sort',
        {
            'timeseries': timeseries,
            'geom': geom
        },
        {
            'firings_out': firings_out
        },
        pp,
        opts
    )


def compute_cluster_metrics(*, timeseries, firings, metrics_out, samplerate, opts={}):
    metrics1 = mlp.addProcess(
        'ms3.cluster_metrics',
        {
            'timeseries': timeseries,
            'firings': firings
        },
        {
            'cluster_metrics_out': True
        },
        {
            'samplerate': samplerate
        },
        opts
    )['outputs']['cluster_metrics_out']
    metrics2 = mlp.addProcess(
        'ms3.isolation_metrics',
        {
            'timeseries': timeseries,
            'firings': firings
        },
        {
            'metrics_out': True
        },
        {
            'compute_bursting_parents': 'true'
        },
        opts
    )['outputs']['metrics_out']
    return mlp.addProcess(
        'ms3.combine_cluster_metrics',
        {
            'metrics_list': [metrics1, metrics2]
        },
        {
            'metrics_out': metrics_out
        },
        {},
        opts
    )


def automated_curation(*, firings, cluster_metrics, firings_out, opts={}):
    # Automated curation
    label_map = mlp.addProcess(
        'ms4alg.create_label_map',
        {
            'metrics': cluster_metrics
        },
        {
            'label_map_out': True
        },
        {},
        opts
    )['outputs']['label_map_out']
    return mlp.addProcess(
        'ms4alg.apply_label_map',
        {
            'label_map': label_map,
            'firings': firings
        },
        {
            'firings_out': firings_out
        },
        {},
        opts
    )

def pyms_anneal_segs(*,timeseries_list, firings_list, firings_out, dmatrix_out, k1_dmatrix_out, k2_dmatrix_out, dmatrix_templates_out, time_offsets, opts={}):

    return mlp.runProcess(
        'pyms.anneal_segments',
        {
            'timeseries_list':timeseries_list,
            'firings_list':firings_list
        },
        {
            'firings_out':firings_out,
            'dmatrix_out':dmatrix_out,
            'k1_dmatrix_out':k1_dmatrix_out,
            'k2_dmatrix_out':k2_dmatrix_out,
            'dmatrix_templates_out':dmatrix_templates_out
        },
        {
            'time_offsets':time_offsets
        },
        opts
    )


def synthesize_sample_dataset(*, dataset_dir, samplerate=30000, duration=600, num_channels=4, opts={}):
    if not os.path.exists(dataset_dir):
        os.mkdir(dataset_dir)
    M = num_channels
    mlp.addProcess(
        'ephys.synthesize_random_waveforms',
        {},
        {
            'geometry_out': dataset_dir + '/geom.csv',
            'waveforms_out': dataset_dir + '/waveforms_true.mda'
        },
        {
            'upsamplefac': 13,
            'M': M,
            'average_peak_amplitude': 100
        },
        opts
    )
    mlp.addProcess(
        'ephys.synthesize_random_firings',
        {},
        {
            'firings_out': dataset_dir + '/firings_true.mda'
        },
        {
            'duration': duration
        },
        opts
    )
    mlp.addProcess(
        'ephys.synthesize_timeseries',
        {
            'firings': dataset_dir + '/firings_true.mda',
            'waveforms': dataset_dir + '/waveforms_true.mda'
        },
        {
            'timeseries_out': dataset_dir + '/raw.mda.prv'
        }, {
            'duration': duration,
            'waveform_upsamplefac': 13,
            'noise_level': 10
        },
        opts
    )
    params = {
        'samplerate': samplerate,
        'spike_sign': 1
    }
    with open(dataset_dir + '/params.json', 'w') as outfile:
        json.dump(params, outfile, indent=4)


def get_time_offset_list(dataset_list):

    offset_list = [0]

    for dataset_dir in dataset_list:

        with open(dataset_dir + '/raw.mda', 'r') as fid:
            params = np.fromfile(fid, 'int32', 5)
            offset_list.append(params[-1] + offset_list[-1])

    return offset_list[0:-1]