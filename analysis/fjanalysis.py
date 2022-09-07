#!/usr/bin/env python3 

import sys
import os
import numpy as np
import pandas as pd
import scipy.stats
import json

import matplotlib.pyplot as plt
import matplotlib as mpl
from tqdm import tqdm
from tabulate import tabulate

import datapaths

import command 
import matdef
import stats 
import twutils
import plotutils
import shapeplot
import _fj
import astrack
import filesystem

import twanalyse

#
# plotutils.default_style()

# 

#
def step_through(trs, plot_f):
    # plot_f takes a single Track
    stepper = make_stepper(trs, plot_f)
    for _ in stepper:
        plt.show()

def make_stepper(trs, plot_f):
    trs = iter(trs)
    while True:
        tr = next(trs)
        yield plot_f(tr)


#################################################################################
# linearized tracks

@command.defaultsave()
def ltr_vel(): 
    debug = 100
    idx, trs = _fj.slicehelper.load_trs(default_crawling_data_list, debug)
    ltrs = [_fj.linearize(tr) for tr in trs]
    # v_tau_max = 5.0
    # for tr in ltrs:
    #     tr.filter_to_resolution(2*v_tau_max)
    
    allvel = np.concatenate([twanalyse._inst_vel(tr) for tr in ltrs])
    print ('mean ', np.mean(allvel))
    allvel = twutils.trim_tail(allvel, 0.01)
    ax = plt.gca()
    color = iter(plotutils.prop_cycle_color)
    plotutils.ax_kdeplot(ax, allvel, linekw={'label':'all velocities'})
    histstyle = {'alpha':0.25, 'density':True}
    ax.hist(allvel, bins=100, color=next(color), **histstyle)

    track_mean = np.array([np.mean(twanalyse._inst_vel(tr)) for tr in ltrs])
    plotutils.ax_kdeplot(ax, track_mean, linekw={'label':'track mean'})
    ax.hist(track_mean, bins=10, color=next(color), **histstyle)
    ax.set_xlabel(r'velocity ($\mu m s^{-1}$)')
    ax.set_ylabel(r'P')
    vstyle = {'c':'k', 'alpha':0.8, 'linestyle':'--', 'linewidth':2}
    plt.tight_layout()
    ax.legend()


#################################################################################
# experimental data equivalent to twanalyse.summary

# CACHE this data for Fanjin dataset
_cache_at = "fanjin/working_copy/summary/track_{:04d}.json"
def lsummary(ltrs, config={}):
    sd = {}
    #
    if all([len(ltr.step_idx) < 10 for ltr in ltrs]):
        sd['failed'] = True
        sd['failed_condition'] = "step_condition"
        return sd

    sd.update(twanalyse.observables(ltrs))
    sd.update(twanalyse._qaparams(ltrs)) # stats.keep
    sd.update(twanalyse.out_of_plane(ltrs))
    compute_kmsd = config.get('compute_kmsd', True)
    if compute_kmsd:
        good = [tr for tr in ltrs if kmsd_condition(tr)]
        if len(good) > 0:
            sd.update(twanalyse.kmsd_summary(good))
    # sd.update(twanalyse.relinearize(ltrs, 'double_ltrs', 0.24))
    # sd.update(twanalyse.relinearize(ltrs, 'quad_ltrs', 0.48))
    # sd.update(twanalyse.relinearize(ltrs, 'cell_ltrs', 1.0))
    return sd


def kmsd_condition(tr): # condition for computing   msd
    lthreshold = 10.0
    dx = tr.get_step_dx()
    dist = np.sqrt(np.sum(dx**2, axis=1))
    cd = np.sum(dist) > lthreshold
    return cd

import pili
cache_form = os.path.join(pili.root, "../"+_cache_at)
def cache_summary():
    filesystem.safemkdir(cache_form)
    # load fanjin
    print("load fanjin data")
    all_idx, ltrs = _fj.slicehelper.load_linearized_trs('all')
    for idx, ltr in zip(all_idx, ltrs):
        sd = lsummary([ltr])
        file = cache_form.format(idx)
        with open(file, 'w') as f:
            json.dump(sd,f)
 
def load_summary():
    target_idx = _fj.slicehelper.load("all")
    data = []
    for idx in target_idx:
        target = cache_form.format(idx)
        with open(target, 'r') as f:
            summary = json.load(f)
        data.append(summary)
    return data

# cache and load summary

def summary():
    cache_summary()
    return load_summary()

def compute_reference_data(ltrs, reference_idx, objectives):
    # refdata: {subset: {objective: weighted mean value}}
    getter = [twutils.make_get(name) for name in objectives]
    refdata = {}
    reftrs = {}
    localdata = {}
    for key, subidx in reference_idx.items():
        reftrs[key] = [ltrs[idx] for idx in subidx]
        localdata[key] = lsummary(reftrs[key])

    for subset, data in reftrs.items():
        nsteps = [len(ltr.step_idx)-1 for ltr in data]
        ldata = localdata[subset]
        refdata[subset] = {}
        for i, objective in enumerate(objectives):
            ref_value = getter[i](ldata)
            refdata[subset][objective] = ref_value
        
    subsets =  list(reftrs.keys())
    _table = {"subset" : subsets}
    _table.update({name: [refdata[subset][name] for subset in subsets]
        for name in objectives})
    df = pd.DataFrame(_table)
    return df


#################################################################################
# use simulation functions to analyse FJ data

DEBUG = None

default_crawling_data_list = 'default_crawling_list'
default_data_list = 'all'

# NOTE: this "velocity" autocorrelation is actually speed autocorrelation
# speed autocorrelation might be related to pili retraction cycles
# whereas velocity autocorrelation would include directional changes
@stats.keep
def vel_corr():
    trackidx, trackmanager = _fj.slicehelper.load_data(default_data_list, DEBUG)
    trs = trackmanager.get_tracklike()

    # arbitrary time
    fit_maxtime = 4.

    print("Computing velocties")
    vs = [twanalyse._vel(tr) for tr in tqdm(trs)]
    plt.clf()
    deltat = trs[0].tstep
    nsteps = int(np.round(fit_maxtime/deltat))
    # ...
    # Computing Velocity Correlations
    print("Calculating velocity autocorrelations")
    v_corr = [np.correlate(v,v,mode='same') for v in tqdm(vs)]
    assert(all((v_corr_i.size > 2*nsteps for v_corr_i in v_corr)))
    mpts = [int(v_corr_i.size/2) for v_corr_i in v_corr]
    each_corr = np.vstack([v_corr_i[mpt+1:mpt+nsteps+1] for mpt, v_corr_i in zip(mpts, v_corr)])
    time = deltat * np.arange(1,nsteps+1,1)
    tau_decay = np.full(len(trs), np.nan)
    assert(tau_decay.size > 0)
    print("Calculating autocorrelation timescale by fitting to exp function")
    failed_idx = []
    for i, v_corr_i in tqdm(enumerate(each_corr)):
        try:
            a, b, c = twutils.fit_decay(time, v_corr_i, plot=False)
            tau_decay[i] = -c
        except RuntimeError as e:
            failed_idx.append(i)
            print("Warning: scipy failed to fit curve to obtain decay time")
            print(str(e))

    return {'all_velocity_decay_time': list(tau_decay)}

# debugging/data visualistion for velocity correlation
def _show_v_corr(tr):
    v = twanalyse._vel(tr)
    v_corr = np.correlate(v,v,mode='full')
    mpt = int(v_corr.size/2)
    fit_maxtime = 10
    nsteps = int(np.round(fit_maxtime/tr.tstep))
    plt.plot(v_corr[mpt:mpt+nsteps])

#
def _displacement_distribution():
    # lets obtain the displacement distribution used to calculate these velocities ....
    idx, trs = _fj.slicehelper.load_trs(default_data_list)
    ld = stats.load()
    # apply time threshold
    avdt = np.array(ld['all_velocity_decay_time'])
    print('min v_decay time {}s'.format(avdt.min()))
    assert(len(avdt)==len(trs))
    for i, v_decay in enumerate(ld['all_velocity_decay_time']):
        trs[i].filter_to_resolution(2*v_decay)

    disps = np.concatenate([twanalyse._disp(tr) for tr in trs])
    # fj raw threshold
    fj_raw_spatial_threshold = 0.03
    below_threshold_fraction = np.count_nonzero([disps < fj_raw_spatial_threshold]) / disps.size
    print("fraction of displacement events below threshold {}".format(below_threshold_fraction))


@command.defaultsave()
def _plot_velocity_distribution():
    # obtain the mean velocity of each track
    ld = stats.load()
    ax = plt.gca()
    track_mean_vel = twutils.trim(np.array(ld['vel']['each_mean']), 0.01)
    track_median_vel = twutils.trim(np.array(ld['vel']['each_median']), 0.01)
    res = 100
    linedct = plotutils.ax_kdeplot(ax, track_mean_vel, res=res, linekw={'label':'per-track mean velocity'})
    linedct = plotutils.ax_kdeplot(ax, track_median_vel, res=res, linekw={'label':'per-track median velocity'})
    # full velocity distribution
    idx, trs = _fj.slicehelper.load_trs(default_data_list)
    avdt = np.array(ld['all_velocity_decay_time'])
    print('min v_decay time {}s'.format(avdt.min()))
    assert(len(avdt)== len(trs))
    for i, v_decay in enumerate(ld['all_velocity_decay_time']):
        trs[i].filter_to_resolution(2*v_decay)
    vs = np.concatenate([twanalyse._vel(tr) for tr in trs])
    # trim the high velocity tail 
    # TODO: filter by max velocity?
    trimmed_vs = twutils.asytrim(vs, 0.0, 0.05)
    plotutils.ax_kdeplot(ax, trimmed_vs, res=200, linekw={'label':'all crawling tracks'})
    
    ax.set_ylabel('P')
    ax.set_xlabel(r'velocity $\mu m s^{-1}$')
    ax.legend()
    

@stats.keep
def vel():
    # compute velocity distribution
    # idx, trs = _fj.slicehelper.load_trs(default_data_list)
    idx, trs = _fj.slicehelper.load_trs(default_crawling_data_list)
    ld = stats.load()
    tau_decay = np.array(ld["all_velocity_decay_time"])[idx]
    for i, tr in enumerate(trs):
        tr.filter_to_resolution(2*tau_decay[i])
    sd = {}
    sd['vel'] = stats.col_stats(trs, twanalyse._vel)
    _vel(trs)
    return sd

@command.defaultsave()
def _vel(trs):
    dt = np.concatenate([tr['time'][1:] - tr['time'][:-1] for tr in trs])
    disp = np.concatenate([_disp(tr) for tr in trs])
    allvel = disp/dt
    # allvel = np.concatenate([twanalyse._vel(tr) for tr in trs])
    plotutils.kdeplot(allvel, res=40, linekw={'label':'velocity'})
    ax = plt.gca()
    # compare median and mean velocities
    # weighted mean velocity
    dt = np.concatenate([tr['time'][1:] - tr['time'][:-1] for tr in trs])
    disp = np.concatenate([twanalyse._disp(tr) for tr in trs])
    true_mean = np.sum(dt*disp)/np.sum(dt)
    #
    plt.axvline(true_mean, label='mean', c='k')
    plt.axvline(np.median(allvel), label='median', c='g')

    print('true mean', true_mean)
    print('mean', np.mean(allvel))
    print('median', np.median(allvel))
    ax.set_ylabel('P')
    ax.set_xlabel(r'velocity $\mu m s^{-1}$')
    ax.legend()
    return allvel

@command.defaultsave()
def raw_velocity():
    debug = None
    idx, trs = _fj.slicehelper.load_trs(default_crawling_data_list, debug)
    vs = [twanalyse._vel(tr) for tr in trs]
    return _raw_velocity(vs)

def _raw_velocity(vs):
    raw_v_per = np.array([np.mean(v) for v in vs])
    raw_v =  np.concatenate(vs)
    raw_v = np.clip(raw_v, np.quantile(raw_v, 0.05), np.quantile(raw_v, 0.95))
    ax = plt.gca()
    plotutils.ax_kdeplot(ax, raw_v_per, linekw={'label':'track mean raw velocity'})
    plotutils.ax_kdeplot(ax, raw_v, linekw={'label':'raw velocity'})
    ax.set_ylabel('P')
    ax.set_xlabel(r'velocity $\mu m s^{-1}$')
    ax.legend()
    sd = {}
    sd['raw_vel'] = np.mean(raw_v)
    twutils.print_dict(sd)
    return sd

@command.defaultsave()
@stats.keep
def kmsd(use_overlap=False):
    print('Using overlap for kmsd calculation: ', use_overlap)
    varname = 'kmsd_with_overlap' if use_overlap else 'kmsd' 
    # safe to nest modifiers like this?
    @stats.keep
    def compute():
        trackidx, trs = _fj.slicehelper.load_trs(default_data_list, DEBUG)
        return {varname: list(_kmsd(trs, use_overlap))}
    kmsd_arr = np.array(stats.load_or_compute(varname, compute))
    ax = plt.gca()
    ax.hist(kmsd_arr, density=True, alpha=0.2)
    plotutils.ax_kdeplot(ax, kmsd_arr[np.isfinite(kmsd_arr)])
    ax.set_xlabel('P')
    ax.set_ylabel(r'$k_\mathrm{MSD}$')
    sd = {}
    sd[varname] = list(kmsd_arr)
    return sd

def _kmsd(trs, use_overlap=False):
    def msd_no_overlap(tr):
        maxwindow = tr.size//10
        return astrack.msd_no_overlap(tr, maxwindow, use_overlap=use_overlap)
    print("Computing kmsd...")
    # TODO: check goodness of fit by eye // find statistical method for this
    kmsd = [msd_no_overlap(tr) for tr in tqdm(trs)]
    return np.array(kmsd)

#################################################################################
# outliers


# run this from fanjin/working_copy

def make_high_displacement_list():
    crawling_idx = _fj.slicehelper.load('default_crawling_list')
    ld = stats.load()
    disp_length = np.array(ld['displacement_length'])
    high_score = np.quantile(disp_length, 0.90)
    high_displacemnt_idx =  crawling_idx[disp_length[crawling_idx] > high_score]
    print('high displacement crawling data set with {} entries'.format(high_displacemnt_idx.size))
    print("mean displacement_length = ", np.mean(disp_length[high_displacemnt_idx]))
    _fj.slicehelper.save('high_displacement_crawling', high_displacemnt_idx)

# NOTE: this function outputs local.json data to cwd
def outliers():
    """
    An outlier will be identified according to one or more statistics.
    For those statistics, outliers will be to the exteme of some max or min score
    The limiting score will be chosen according to some reasoning about the acceptable values of the statistic 
    or be chosen as a percentile of the data i.e. 1%. 
    Statistics should be computed on the entire dataset so that they can be compared.

    Outlier functions should take the tracking data, load or compute the statistic, the data used to select outliers will be arguments


    External computes: vel_corr, kmsd
    """
    debug = None
    idx, trs = _fj.slicehelper.load_trs('all', debug)
    def maxtime_outliers(idx, trs, minsteps=1000, maxsteps=None):
        outidx = idx[np.array([tr.size < minsteps for tr in trs])]
        return outidx, np.array([]), np.array([])

    def displacement_outliers(idx, trs, minlength=1.0, maxlength=None):
        @stats.keep 
        def compute():
            def displacement_length(tr):
                x, y = tr['x'], tr['y']
                return np.sqrt( (x[0]-x[-1])**2 + (y[0]-y[-1])**2 )
            dlength = np.array([displacement_length(tr) for tr in trs])
            return {'displacement_length': list(dlength)}
        dlength = np.array(stats.load_or_compute('displacement_length', compute))
        outidx = idx[dlength < minlength]
        assert(not np.any(np.isnan(dlength)))
        return outidx, np.array([]), np.array([])

    def aspect_ratio_outliers(idx, trs, min_min_aspect=1.6, max_min_aspect=None):
        # sliding window coarse graining
        def _coarse_graining(arr, framewindow=200):
            # print('Uniform coarse graining with width = {}'.format(framewindow))
            if arr.size <= framewindow:
                return np.nan
            frh = framewindow//2
            sbasis = np.arange(frh, arr.size-frh, 1, dtype=int)
            return np.array([np.mean(arr[j-frh:j+frh-1]) for j in sbasis])
        @stats.keep
        def compute():
            min_lwratio = np.array([np.min(_coarse_graining(tr['length']/tr['width'])) for tr in tqdm(trs)])
            return {'min_aspect_ratio':min_lwratio.tolist()}
        min_lwratio = np.array(stats.load_or_compute('min_aspect_ratio', compute))
        outidx = idx[min_lwratio < min_min_aspect]
        # print('---aspect low quantile ', _quantileofscore(min_lwratio, min_min_aspect))
        assert(not np.any(np.isnan(min_lwratio)))
        return outidx, np.array([]), np.array([])
        
    def vel_corr_outliers(idx, trs, minvcorr=0.5, maxvcorr=5.0):
        v_decay = np.array(stats.load_or_compute('all_velocity_decay_time', vel_corr))[:debug]
        isnan = np.isnan(v_decay)
        nan_idx = idx[isnan]
        safe_v_decay = v_decay[~isnan]
        safe_idx = idx[~isnan]
        left_outidx = safe_idx[safe_v_decay < minvcorr]
        right_outidx = safe_idx[safe_v_decay > maxvcorr]
        return left_outidx, right_outidx, nan_idx
        
    def kmsd_outliers(idx, trs, minkmsd=0.8, maxkmsd=None):
        @stats.keep
        def compute():
            return {'kmsd': list(_kmsd(trs))}
        kmsd_arr = np.array(stats.load_or_compute('kmsd', compute))[:debug]
        # first remove np.nan
        isnan = np.isnan(kmsd_arr)
        nan_idx = idx[isnan]
        safe_kmsd_arr, safe_idx = kmsd_arr[~isnan], idx[~isnan]
        outidx = safe_idx[safe_kmsd_arr < minkmsd]
        return outidx, np.array([]), nan_idx

    methods = [
        ('total_time_outlier', maxtime_outliers, int(100/matdef.TIMESTEP), None),
        ('displacement_outlier', displacement_outliers, 1.0, None),
        ('aspect_ratio_outlier', aspect_ratio_outliers, 1.6, None),
        ('vel_corr_outlier', vel_corr_outliers, 0.5, 5.0),
        ('kmsd_outlier', kmsd_outliers, 0.8, None) 
        ]
        
    sd = {}
    for name, method, low_thresh, high_thresh in methods:
        print(name)
        leftidx, rightidx, nanidx = method(idx, trs, low_thresh, high_thresh)
        print('low_tresh', low_thresh, 'high_thresh', high_thresh)
        count_form = 'left\t{:4d}\tright{:4d}\tnan{:4d}'
        print(count_form.format(len(leftidx), len(rightidx), len(nanidx)))
        sd_part = {
            name:{
                'left':leftidx.tolist(), 'right':rightidx.tolist(), 'nan':nanidx.tolist()
                },
            name+'_threshold': [low_thresh, high_thresh]
            }
        sd.update(sd_part)
    stats.extend(sd)

    to_remove_list = [
        sd['total_time_outlier']['left'],
        sd['aspect_ratio_outlier']['left'],
        sd['vel_corr_outlier']['left'],
        sd['vel_corr_outlier']['right'],
        sd['kmsd_outlier']['left']
    ]
    to_remove_list.extend([sd[k]['nan'] for k in [method[0] for method in methods]])
    to_remove = set()
    for outliers in to_remove_list:
        for i in outliers:
            to_remove.add(i)
    keep_idx = [i for i in idx if i not in to_remove]
    # write out the crawling list
    _fj.slicehelper.save('default_crawling_list', np.array(keep_idx))

    # compute the walking list
    walking_idx = sd['aspect_ratio_outlier']['left']
    _fj.slicehelper.save('default_walking_list', np.array(walking_idx))

    # save low displacement list
    low_disp = sd['displacement_outlier']['left']
    _fj.slicehelper.save('low_displacement_outliers', np.array(low_disp))


def save_median_track():
    # median by mean velocity 
    debug = None
    idx, trs = _fj.slicehelper.load_trs(default_crawling_data_list, debug)
    ld = stats.load()
    def weighted_average(a, t):
        return np.sum(t*a)/np.sum(t)
    tr_mean_vel = [weighted_average(twanalyse._disp(tr), tr.get_dt()) for tr in trs]
    sort_idx = np.argsort(tr_mean_vel)
    cpt = int(sort_idx.size/2)
    dn = 10
    median_idx = sort_idx[cpt-dn:cpt+dn]

    median_trs = [trs[i] for i in median_idx]
    median_disp = np.mean([ld['displacement_length'][i] for i in median_idx])
    print('median track displacement', median_disp)
    print('median track velocity ', np.mean([tr_mean_vel[i] for i in median_idx]))
    # v_decay = ld['all_velocity_decay_time']
    # print('median track v_decay', v_decay[median_idx] )
    # print('median track size', np.sum([tr.size for tr in median_trs]))
    # print('median track by velocity ', median_idx)
    _fj.slicehelper.save('vel_median', median_idx)

    # lets save a high velocity data set aswell
    high_idx = sort_idx[-50:]
    high_trs = [trs[i] for i in high_idx]
    high_disp = np.mean([ld['displacement_length'][i] for i in high_idx])
    print('high track displacement', high_disp)
    print('high track velocity ', np.mean([tr_mean_vel[i] for i in high_idx]))
    _fj.slicehelper.save('vel_high', high_idx)



@command.defaultsave()
def median_track_vel_distribution():
    ld = stats.load()
    v_decay = ld['all_velocity_decay_time']
    idx, trs = _fj.slicehelper.load_trs('vel_median')
    ltrs = [_fj.linearize(median_tr) for median_tr in trs]

    vel_arr = np.concatenate([twanalyse._disp(tr)/tr.get_dt() for tr in trs])
    vel_arr = twutils.trim(vel_arr,0.05)
    ax = plt.gca()
    res = 100
    # lims = (0.0, 0.04)
    lims = (None, None)
    plotutils.ax_kdeplot(ax, vel_arr, res=res, linekw={'label':'base'}, xlims=lims)

    vel_arr = np.concatenate([twanalyse._disp(tr)/tr.get_dt() for tr in ltrs])
    vel_arr = twutils.trim(vel_arr,0.05)
    plotutils.ax_kdeplot(ax, vel_arr, res=res, linekw={'label':'linearized'}, xlims=lims)
    ax.legend()

@command.defaultsave()
def high_track_vel_distribution():
    ld = stats.load()
    v_decay = ld['all_velocity_decay_time']
    idx, trs = _fj.slicehelper.load_trs('vel_high')
    ltrs = [_fj.linearize(high_tr) for high_tr in trs]

    vel_arr = np.concatenate([twanalyse._disp(tr)/tr.get_dt() for tr in trs])
    vel_arr = twutils.trim(vel_arr,0.05)
    ax = plt.gca()
    res = 100
    # lims = (0.0, 0.04)
    lims = (None, None)
    plotutils.ax_kdeplot(ax, vel_arr, res=res, linekw={'label':'base'}, xlims=lims)

    vel_arr = np.concatenate([twanalyse._disp(tr)/tr.get_dt() for tr in ltrs])
    vel_arr = twutils.trim(vel_arr,0.05)
    plotutils.ax_kdeplot(ax, vel_arr, res=res, linekw={'label':'linearized'}, xlims=lims)
    ax.legend()


def _quantileofscore(a, score):
    if score == None:
        return 0.0
    return np.count_nonzero(a<=score)/len(a)

def show_outliers():
    # Natural order to filtering shall be
    # total time, aspect_ratio, vel_corr, kmsd, displacement_length 
    debug = None
    all_examples = False
    idx, trs = _fj.slicehelper.load_trs(default_data_list, debug)
    ld = stats.load()
    outlier_name = ['total_time_outlier', 'displacement_outlier', 'aspect_ratio_outlier', 'vel_corr_outlier', 'kmsd_outlier']
    statistic_name = ['displacement_length', 'min_aspect_ratio', 'all_velocity_decay_time', 'kmsd']
    tr_size = np.array([tr.size for tr in trs])
    statistic = [tr_size] + [np.array(ld[name]) for name in statistic_name]

    # create table of outliers relative to the whole dataset
    table_header = ['variable', 'min', 'max', 'min_q', 'max_q', 'min_outliers', 'max_outliers', 'nan', 'fraction']
    frac_form = '{:d}/{:d}'
    def table_row(name, stat):
        low, high = ld[name+'_threshold']
        isnan = np.isnan(stat)
        nan = np.count_nonzero(isnan)
        safe_stat = stat[~isnan]
        qlow = _quantileofscore(stat, low) 
        qhigh = _quantileofscore(stat, high) if high != None else 1.0
        left = np.count_nonzero(safe_stat < low) if low != None else 0
        right= np.count_nonzero(safe_stat > high) if high != None else 0
        if low == None: low = safe_stat.min()
        if high == None: high = safe_stat.max()
        total = left + right + nan
        return [name, low, high, qlow, qhigh, left, right, nan, frac_form.format(total,len(trs))]
    table = [table_row(name, stat) for name, stat in zip(outlier_name, statistic) ]
    print(tabulate(table, headers=table_header))
    #
    # filter the data using the predefined order 
    order = ['total_time_outlier', 'aspect_ratio_outlier', 'vel_corr_outlier', 'kmsd_outlier',  'displacement_outlier']
    #
    x_label = ['Total Time', 'Displacement Length', 'Min Aspect Ratio', 'Velocity Autocorrelation Time', r'$k_\mathrm{MSD}$']
    #
    presentation = True
    def _distribution(name, safe_stat, x_label, lims):
        low, high = lims
        fig = plt.figure()
        ax = fig.gca()
        # trim_safe_stat = twutils.trim_tail(safe_stat, 0.01)
        trim_safe_stat = safe_stat
        plotutils.ax_kdeplot(ax, trim_safe_stat, linekw={'linewidth':8})
        if presentation:
            vstyle = {'c':'k', 'alpha':0.8, 'linestyle':'--', 'linewidth':4}
        else:
            vstyle = {'c':'k', 'alpha':0.8, 'linestyle':'--'}

        if low != None: 
            ax.axvline(low, **vstyle)
        if high != None:
            ax.axvline(high, **vstyle)
        ax.set_xlim(xmax=1.05 * high)
        fig._plot_name = name

        if presentation:
            ax.set_xticks([low, np.mean(safe_stat), high])
            ax.set_yticklabels([])
            ax.set_title(x_label, fontsize='40')
        else:
            ax.set_xlabel(x_label)
        ax.set_ylabel('P')

        plt.tight_layout()
    #
    def _outlier_tracks(trs, name, low_idx, high_idx, nan_idx, stat):
        # plotutils.default_style()
        nameform = 'track_{:04d}_{:04d}.png'
        fig = plt.figure() # temporary figure
        ax = fig.gca()
        idx_list = [low_idx, high_idx, nan_idx]
        pdir_part_list = ['low/', 'high/', 'nan/']
        for idx, pdir_part in zip(idx_list, pdir_part_list):
            if len(idx):
                pdir = os.path.join(command.pdir, '_'.join([name, pdir_part]))
                datapaths.force_mkdir(pdir)
                for seq_i, score_i in  enumerate(np.argsort(stat[idx])):
                    i = idx[score_i]
                    shapeplot.longtracks(ax, [trs[i]])
                    out = os.path.join(pdir, nameform.format(seq_i, i))
                    print('Plotting outlier to {}'.format(out))
                    plt.savefig(out)
                    ax.clear()
        plt.close()

    def table_generator(order):
        keep_idx = idx
        for name in order:
            print(name)
            i = outlier_name.index(name)
            # stat now depends on keep_idx
            stat = statistic[i][keep_idx]
            # this_trs = [trs[tr_i] for tr_i in keep_idx]
            low, high = ld[name+'_threshold']
            isnan = np.isnan(stat)
            nan_idx = keep_idx[isnan]
            safe_stat = stat[~isnan]
            if low is None: low = safe_stat.min()
            if high is None: high = safe_stat.max()
            safe_keep_idx = keep_idx[~isnan]
            low_idx = safe_keep_idx[safe_stat < low]
            high_idx = safe_keep_idx[safe_stat > high]
            keep_idx = safe_keep_idx[np.logical_and(safe_stat>=low, safe_stat<=high)]
            #
            # (i) plot distributions 
            _distribution(name, safe_stat, x_label[i], (low, high))
            # (ii) plot example tracks
            if all_examples:
                _outlier_tracks(trs, name, low_idx, high_idx, nan_idx, statistic[i])
            #
            qlow = _quantileofscore(safe_stat, low)
            qhigh = _quantileofscore(safe_stat, high) 
            total = len(low_idx) + len(high_idx)
            yield (name, low, high, qlow, qhigh, len(low_idx), len(high_idx), len(nan_idx), frac_form.format(total,len(stat)))

    table_gen = table_generator(order)
    ordered_table = [row for row in table_gen]
    print("---Ordered Filtering---")
    print(tabulate(ordered_table, headers=table_header))
    extra = 'present' if presentation else ''
    command.save_all_figures(show_outliers, svg=True, extra=extra)

if __name__=='__main__':
    ff, thisargs = command.process(locals())
    ff(*thisargs)