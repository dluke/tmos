#!/usr/bin/env python3

import sys, os
import collections 
import itertools
import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import patches
import seaborn as sns

import wrinterval as wr
import parameters

import twutils 
import plotutils
import readtrack
import stats
from tqdm import tqdm

import command
from command import defaultsave

"""
some simple pili analysis methods followed by a method which draws
pili attachment points on hexsphere surface
"""


# remove the +ve tail from exponential data to make it easier to plot
from twutils import trim_tail

#################################################################################
# utilities

def get_attach(pr):
    xyz = np.column_stack(
        [pr['attach_x'], pr['attach_y'], pr['attach_z']]
        )
    return xyz

def get_anchor(pr):
    xyz = np.column_stack(
        [pr['anchor_x'], pr['anchor_y'], pr['anchor_z']]
        )
    return xyz

def get_axis(pr):
    xyz = np.column_stack(
        [pr['paxis_x'], pr['paxis_y'], pr['paxis_z']]
        )
    return xyz


#################################################################################

def wtavg(track, col):
    """compute the time weighted average of a column"""
    time = track['time']
    deltime = time - np.insert(time[:-1], 0, 0.)
    return np.sum(deltime * track[col])/time[-1]

pv_ret = 0.75 # retraction speed
d_bound = 0.004


# excess_length on attachment
@command.defaultsave()
def attach_excess_length():
    evtrs = readtrack.eventset()
    dds = [_excess_length(evtr) for evtr in evtrs]

    att_dx = np.concatenate([dd['att_dx'] for dd in dds])
    iexcess = np.concatenate([dd['inter_excess'] for dd in dds])
    print('mean inter_excess = {}'.format(np.mean(iexcess)))
    print('mean att_dx = {}'.format(np.mean(att_dx)))
    excess_time = np.mean(iexcess)/pv_ret
    print('thats an average excess retraction time of {:10.2f}ms at {}\mu m/s'.format(excess_time*1000., pv_ret))
    one_extend = np.count_nonzero(iexcess < d_bound)
    print('Attached by one extend {}/{}'.format(one_extend, len(att_dx)))

    plt.hist(att_dx, bins=200, alpha=0.70, density=True)
    plt.hist(att_dx, bins=25, alpha=0.35, density=True)

def _excess_length(evtr):
    attidx = evtr['trigger'] == 'attach'
    iexcess = evtr['inter_excess'][attidx]
    pl = evtr['plength'][attidx]
    leq = evtr['pleq'][attidx]
    dd = {}
    dd['att_dx'] = leq-pl
    dd['inter_excess'] = iexcess
    return dd


# check the shortening of pili on detachment events
# TODO move some of this information into summary and save it with stats
def shortening():
    evtrs = readtrack.eventset()
    evtr = evtrs[0]
    relidx = np.nonzero(evtr['process'] == 'release')[0]
    assert(np.all(evtr['trigger'][relidx] == 'release'))
    shortening = evtr['last_shortening'][relidx]
    print('{} release events'.format(relidx.size))

    shortenidx = np.nonzero(shortening > 0)[0]
    print('{} shortening events'.format(shortenidx.size))
    print('thats {:5.2f}%'.format(100 * float(shortenidx.size)/float(relidx.size)))
    mean_short = np.mean(shortening)
    print('mean shortening {}'.format(mean_short))
    print('mean shortening of shortening events {}'.format(np.mean(shortening[shortenidx])))
    _plot_shortening(shortening)

@command.defaultsave()
def _plot_shortening(shortening):
    ax = plt.gca()
    ax.set_xlabel('shortening distance in microns')
    ax.set_ylabel('')
    ax.hist(shortening)

# intersection_excess length during detachment attempts
# once release event probability is fixed there will be no failed release events
# inter_excess is the excess length of the chain inside the surface
# inter_excess already has the d_bound = 0.004 subtracted from the shrink() event
# excess_length is leq - pl for an attached pilus 

@command.defaultsave()
def failed_release():
    evtrs = readtrack.eventset()
    dds = [_failed_release(evtr) for evtr in evtrs]

    relidx_size = sum([dd['relidx_size'] for dd in dds])
    one_shrink_off = sum([dd['one_shrink_off'] for dd in dds])
    f_inter_excess = np.concatenate([dd['f_inter_excess'] for dd in dds])
    print('Failed release events {}/{}'.format(f_inter_excess.size,relidx_size))
    print('Failed by one shrink {}/{}'.format(one_shrink_off,f_inter_excess.size))

    plt.hist(f_inter_excess, bins=200, alpha=0.35, density=True)
    plt.xlabel('intersection excess (microns)')
    plt.hist(f_inter_excess, bins=25, alpha=0.7, density=True)

def _failed_release(evtr):
    # count the failed release events
    relidx = np.nonzero(evtr['process'] == 'release')[0]
    f_relidx = np.nonzero(evtr['trigger'] == 'no_release')[0]
    f_pidx = evtr['pidx'][f_relidx]
    f_inter_excess = evtr['inter_excess'][f_relidx]
    one_shrink_off = np.count_nonzero(evtr['inter_excess'][f_relidx] < d_bound)
    # what can we learn from checking excess_length?
    # f_excess_length = evtr['pleq'][f_relidx] - evtr['plength'][f_relidx]
    dd = {}
    dd['relidx_size'] = relidx.size
    dd['f_inter_excess'] = f_inter_excess
    dd['one_shrink_off'] = one_shrink_off 
    return dd



# types of attachment
def _classify(evtr, debug_release=False):
    attidx = evtr['trigger'] == 'attach'
    # types extension, release, mdstep
    exttype = np.nonzero(evtr['process'][attidx] == 'extension')[0]
    resample_type = np.nonzero(evtr['process'][attidx] == 'resample')[0]
    mdtype = np.nonzero(evtr['process'][attidx] == 'mdstep')[0]
    
    # map back to the whole array
    atts = np.nonzero(attidx)[0]
    mdidx_p = atts[mdtype] - 1
    release_last = evtr['process'][mdidx_p] == 'release'
    samestep = evtr['time'][mdidx_p+1] == evtr['time'][mdidx_p]
    samepilus = evtr['pidx'][mdidx_p+1] == evtr['pidx'][mdidx_p]
    is_release = np.logical_and.reduce((release_last, samestep, samepilus))
    if debug_release and np.count_nonzero(is_release) > 0:
        relidx = atts[mdtype[is_release]]
        for time, p, ext_motor, ret_motor in zip(evtr['time'][relidx], evtr['pidx'][relidx], evtr['ext_motor'][relidx], evtr['ret_motor'][relidx]):
            print('Debug: attach event triggered by release at time {:8.3f} pilus {:6d} motors (ext,ret)=({},{})'.format(time, p, ext_motor, ret_motor))
        
    releasetype = mdtype[is_release]
    mdsteptype = mdtype[~is_release]
    
    return {
        'extension' : atts[exttype],
        'resample' : atts[resample_type],
        'release' : atts[releasetype],
        'mdstep' : atts[mdsteptype]
        }

# note that mdstep has a significant count and large values of inter_excess. This is troubling.
# TODO: move this to a jnotebook for analysis ...
def count_attach():
    # not the most appropriate style
    # plt.style.use(plotutils.get_style('bigax')) 
    evtr = readtrack.eventset()[0]
    at_type = _classify(evtr)
    at_count = {k: len(v) for k, v in at_type.items()}
    print(at_count)
    # plot a bar graph using this dictionary

    fig = plt.figure(figsize=(16,10))
    ax = fig.gca()
    ax.set_ylabel('Number of attachments')
    ax.set_xlabel('Cause of attachment')
    ax.set_yticks([])
    fig._plot_name = 'attach_number'
    plotutils.plot_bar(ax, at_count)
    fig.tight_layout()
    #
    # also compute average intersection length of attachments
    plt.close()
    at_inter = { k : np.mean(evtr['inter_excess'][idx]) for k, idx in at_type.items() }
    fig = plt.figure(figsize=(16,10))
    ax = fig.gca()
    ax.set_ylabel('Avg. Intersecton length ($\mu m$)')
    ax.set_xlabel('Cause of attachment')
    ax.locator_params(nbins=4)
    plt.yticks(fontsize=30)
    fig._plot_name = 'inter_length'
    plotutils.plot_bar(ax,at_inter)
    fig.tight_layout()
    
    # also compute average intersection length of attachments
    plt.close()
    at_excess = { k : np.mean(evtr['pleq'][idx] - evtr['plength'][idx]) for k, idx in at_type.items() }
    fig = plt.figure(figsize=(16,10))
    ax = fig.gca()
    ax.set_ylabel('Avg. Excess Length ($\mu m$)')
    ax.set_xlabel('Cause of attachment')
    ax.locator_params(nbins=4)
    plt.yticks(fontsize=30)
    fig._plot_name = 'excess_length'
    plotutils.plot_bar(ax, at_inter)
    fig.tight_layout()

    
22.2148,
#################################################################################
#

#  ----------------------------------------------------------------
# construct pandas data frame for pilus event data

cols = [
    'track_index', 'pilus_index', 'lifetime', 'n_bindings', 
    'first_binding_length', 'first_binding_contraction', 'first_binding_duration', 'first_taut_duration', 'first_binding_delay',
    'mean_binding_length', 'mean_binding_contraction', 'mean_binding_duration', 'mean_taut_duration', 'mean_taut_delay'
]
# seperate dataframe for pili that don't bind the surface
nbcols = [
    'track_index', 'pilus_index', 'lifetime', 'n_bindings', 'first_extension'
]

# TODO
# reorganised data from columns to {pidx : [data] }
def _reorganize(ptr):
    pidx = ptr['pidx']
    pdata = {idx : [] for idx in pidx}
    for row in ptr._track:
        pdata[row['pidx']].append(row)
    return pdata

def eventdf():
    ptrs = readtrack.eventset()
    _eventdf(ptrs)

def _eventdf(ptrs):
    data = {name: []  for name in cols}
    def _append_stats(track_i, ptr):
        N = ptr
        N_binding
        data.append(np.full(track_i))


    for track_i, ptr in enumerate(ptrs):
        _append_stats(track_i, ptr)


#  ----------------------------------------------------------------

# compute taut time and bound times using event data
def lifetime(plotting=False):
    ptrs = readtrack.eventset()
    _lifetime(ptrs, plotting)

@stats.keep
def _lifetime(ptrs,  plotting=False):
    ldata = [evtr_lifetime(ptr) for ptr in tqdm(ptrs)]
    sdls = [ll[0] for ll in  ldata]

    # collect data
    sd = {}
    arr_data = {}
    arr_keys = [
        'lifetime', 'bound_time', 'taut_delay', 'taut_time', 
        'contract_start', 'contract_end', 'contract_length',
        'effective_contract_length', 'first_extension_length']
    def _arr_stats(arr):
        # Note that these are NOT independent measurements
        # unless we introduce a sampling function
        # TODO:  ...
        st = {}
        st['N'] = len(arr)
        st['mean'] = np.mean(arr)
        st['median'] = np.median(arr)
        # st['std'] = np.std(arr)
        return st
    sd['n_pili_lifetimes'] = sum([d['n_pili_lifetimes'] for d in sdls])
    for k in arr_keys:
        arr_data[k] = np.concatenate([d[k] for d in sdls])
        sd[k] = _arr_stats(arr_data[k])
    sd['bound_pili_participation'] = np.mean([d['bound_pili_participation'] for d in sdls])
    sd['taut_pili_participation'] = np.mean([d['taut_pili_participation'] for d in sdls])
    # sd['taut_pili_extension_ratio'] = np.mean([d['taut_pili_extension_ratio'] for d in sdls])
    twutils.print_dict(sd)

    if plotting:
        taut_lt  = ldata[0][1]
        taut_dl = ldata[0][2]
        fig = plt.figure(figsize=(12,9))
        plot_lifetimes(*[arr_data[name] for name in [
            'lifetime', 'bound_time', 'taut_time', 'taut_time']])
        plt.clf()
        plot_contractions(taut_lt, taut_dl)
        plt.clf()
        plot_contraction_distribution(taut_dl)
        plt.clf()

        # plotting distributions
        effective_taut_dl = np.array([np.sum([ct.dl for ct in pt if ct.dl<0.]) for pt in taut_dl.values() if len(pt)>0])
        sns.histplot(effective_taut_dl)
        plt.savefig("plots/effective_contract_length.png")
        plt.clf()

        lifetime = arr_data['lifetime']
        bound_time = arr_data['bound_time']
        taut_time = arr_data['taut_time']
        pili_lifetime_distribution(trim_tail(lifetime,0.01))
        plt.clf()
        plot_bound_time(trim_tail(bound_time,0.01))
        plt.clf()
        plot_tension_time(trim_tail(taut_time,0.01))
        plt.clf()
        # plotting together

    return sd

Lifetime = collections.namedtuple('Lifetime', ['start', 'end', 'duration'])
Contraction = collections.namedtuple('Contraction', ['start', 'end', 'dl'])

# This analysis is a bit out dated. We now simulate almost exclusivly in the limit
# where pili only interact with the surface for  one retraction during their lifetime
# so for example it would make sense to replace individual taut contractions (taut_dl)
# with the total taut change in length for that pilus
# 
def evtr_lifetime(evtr, plotting=False):
    # sort events by pilus
    pd = {}
    # 
    att_lt = {}
    taut_t = {} # time before becoming taut
    taut_lt = collections.OrderedDict() # taut lifetime
    taut_dl = collections.OrderedDict() # taut contraction of length
    full = {}
    first_extension_length = []
    pili_pidx = set(evtr['pidx'])
    for pidx in pili_pidx:
        this = evtr._track[:][evtr['pidx']==pidx]

        # use the first ext_off event to get the extension length
        ext_off = np.where(this['process'] == 'ext_off')[0]
        if ext_off.any():
            _first = ext_off[0]
            leq = this["pleq"][_first]
            first_extension_length.append(leq)


        spawn = this[this['trigger'] == 'spawn']
        dissolve = this[this['trigger'] == 'dissolve']
        assert(spawn.size == 1), 'Pili can only spawn once'
        assert(dissolve.size <= 1)
        if (dissolve.size != 1):
            continue # 'discard pili which outlive the simulation'

        pd[pidx] = this
        time = this['time']
        full[pidx] = Lifetime(time[0], time[-1], time[-1]-time[0])
        binding = this['isbound']
        taut = this['istaut']
        length = this['pleq'] # for taut pili plength > pleq (very slightly) 
    
        t_transition = np.nonzero(np.diff(taut))[0]
        b_transition = np.nonzero(np.diff(binding))[0]
        if b_transition.size > 1:
            # reduce to even number of binding->unbinding pairs
            if b_transition.size % 2 == 1:
                b_transition = b_transition[:-1]
            btt = time[1:][b_transition]
            lts = [Lifetime(btt[i], btt[i+1], btt[i+1]-btt[i]) for i in range(0, b_transition.size, 2)]
            att_lt[pidx] = lts
            #
            ttt = time[1:][t_transition]
            if t_transition.size % 2 == 1:
                t_transition = t_transition[:-1]
            taut_lt[pidx] = [Lifetime(ttt[i], ttt[i+1], ttt[i+1]-ttt[i]) for i in range(0, t_transition.size, 2)]
            tleq = length[1:][t_transition]
            taut_dl[pidx] = [Contraction(tleq[i], tleq[i+1], tleq[i+1]-tleq[i]) for i in range(0,t_transition.size, 2)]

            # find the tension time (Not the taut lifetime)
            for i in range(0, b_transition.size, 2):
                taut_cut = taut[1:][b_transition[i]:b_transition[i+1]+1]
                taut_idx = np.where(taut_cut==1)[0]
                if taut_idx.size > 0:
                    # only take the first time the pilus becomes taut
                    tidx = taut_idx[0]
                    start = time[1:][b_transition[i]]
                    end = time[1:][b_transition[i] + tidx]
                    # print(this[:][1:][b_transition[i]:b_transition[i+1]+1])
                    # print('Dt' ,end-start)
                    # print()
                    if pidx not in taut_t:
                        taut_t[pidx] = []
                    taut_t[pidx].append( Lifetime(start, end, end-start) )


    # times
    ltall = np.array([lt.duration for lt in full.values()])
    bound_time = np.array([lt.duration for pt in att_lt.values() for lt in pt])
    taut_time = np.array([lt.duration for pt in taut_t.values() for lt in pt])
    taut_lifetime = np.array([lt.duration for pt in taut_lt.values() for lt in pt])
    # lengths
    contract_start_l = np.array([ct.start for pt in taut_dl.values() for ct in pt])
    contract_end_l = np.array([ct.end for pt in taut_dl.values() for ct in pt])
    contract_dl = np.array([ct.dl for pt in taut_dl.values() for ct in pt])

    # per pilus rather than per binding/per becoming taut
    effective_taut_dl = np.array([np.sum([ct.dl for ct in pt if ct.dl<0.]) for pt in taut_dl.values() if len(pt)>0])

    # check whether any of these contraction lengths are in fact positive
    # tper = np.count_nonzero(contract_dl > 0)/contract_dl.size
    #
    sd = {}
    sd['n_pili_lifetimes'] = len(pd) 
    sd['lifetime'] = ltall
    sd['bound_time'] = bound_time
    sd['taut_delay'] = taut_time
    sd['taut_time'] = taut_lifetime
    sd['contract_start'] = contract_start_l
    sd['contract_end'] = contract_end_l
    sd['contract_length'] = contract_dl
    sd['effective_contract_length'] = effective_taut_dl
    sd['bound_pili_participation'] = len(att_lt)/len(pd) if len(pd) > 0 else 0
    sd['taut_pili_participation'] = len(taut_t)/len(pd) if len(pd) > 0 else 0
    # sd['taut_pili_extension_ratio'] = tper
    sd['first_extension_length'] = np.array(first_extension_length)
    return sd, taut_lt, taut_dl


@command.defaultsave()
def plot_lifetimes(ltall, bound_time, taut_time, taut_lifetime):
    ax = plt.gca()
    res = 20
    label = iter([r'Pilus lifetime', r'$\tau_\mathrm{tension}$', r'$\tau_\mathrm{dwell}$'])
    plotutils.ax_kdeplot(ax, trim_tail(ltall,0.01), res=res, linekw={'label':next(label)})
    plotutils.ax_kdeplot(ax, trim_tail(bound_time,0.01), res=res, linekw={'label':next(label)})
    plotutils.ax_kdeplot(ax, trim_tail(taut_time,0.01), res=res, linekw={'label':next(label)})
    ax.set_xlabel('time (s)')
    ax.set_ylabel('P')
    ax.legend()


@command.defaultsave()
def plot_contraction_distribution(taut_dl):
    ax = plt.gca()
    all_taut_start = np.array([dl.start for pt in taut_dl.values() for dl in pt])
    all_taut_end= np.array([dl.end for pt in taut_dl.values() for dl in pt])
    linedata = plotutils.ax_kdeplot(ax, all_taut_start, res=20, linekw={'label':r'$l_\mathrm{tension,on}$'})
    plotutils.ax_kdeplot(ax, all_taut_end, res=20, linekw={'label':r'$l_\mathrm{tension,off}$'})
    basis, y = linedata['space'], linedata['pde']
    vis_idx = np.argwhere(basis > 2.0).ravel()
    y0 = np.interp([2.0], basis, y)[0]

    style = {'alpha':0.5}
    fill_x = np.insert(basis[vis_idx], 0, 2.0)
    fill_fraction = np.count_nonzero(all_taut_start > 2.0)/all_taut_start.size
    ax.fill_between(fill_x, np.zeros(fill_x.size), np.insert(y[vis_idx],0,y0), **style)
    y_offset = 0.2*y0
    xy = (2.5, np.interp([2.5], basis, y)[0] + y_offset)
    ax.annotate(r'{:d}% of data > 2.0$\mu m$'.format(int(100*fill_fraction)), xy=xy, xycoords='data')
    ax.set_xlabel(r'length $\mu m$')
    ax.set_ylabel('P')
    ax.legend()
    # ax.set_ylabel(r'length $\mu m$')

@command.defaultsave()
def plot_contractions(taut_lt, taut_dl):
    mpl.style.use('fivethirtyeight')
    ax = plt.gca()
    style = {'linewidth':1}
    # color = itertools.cycle(plotutils.prop_cycle_color)
    all_taut_lt = [lt for pt in taut_lt.values() for lt in pt]
    all_taut_dl = [dl for pt in taut_dl.values() for dl in pt]
    for lt, dl in zip(all_taut_lt, all_taut_dl):
        # cycle default colors
        # c = next(color)
        # fixed color
        c = 'g'
        ax.plot(lt.start, dl.start, marker='D', color=c)
        ax.plot([lt.start, lt.end], [dl.start,dl.end], color=c, **style)
    ax.set_ylabel(r'length ($\mu m$)')
    ax.set_xlabel('time (s)')
    ax.set_ylim(ymin=0.)


@command.defaultsave()
def plot_bound_time(bound_time): 
    ax = plt.gca()
    ax.set_xlabel(r'$\tau_\mathrm{dwell}$')
    ax.hist(bound_time)

@command.defaultsave()
def plot_tension_time(taut_time): 
    ax = plt.gca()
    ax.set_xlabel(r'$\tau_\mathrm{tension}$')
    ax.hist(taut_time)

@command.defaultsave()
def pili_lifetime_distribution(ltall):
    ax = plt.gca()
    ax.hist(ltall)
    ax.set_xlabel('Pilus lifetime (s)')


#################################################################################
# print some basic averages
# use body tracking data instead
def avgnbound():
    pst = readtrack.eventset()
    nbst = [wtavg(tr, 'nbound') for tr in pst]
    print('avgnbound')
    print(nbst)
    print('mean')
    print(np.mean(nbst))

def avgntaut():
    pst = readtrack.eventset()
    nbst = [wtavg(tr, 'ntaut') for tr in pst]
    print('avgntaut')
    print(nbst)
    print('mean')
    print(np.mean(nbst))

def avgnb(ptr):
    nbst = [wtavg(tr.track, 'nbound') for tr in ptr]
    return np.mean(nbst)

def avgnt(ptr):
    nbst = [wtavg(tr.track, 'ntaut') for tr in ptr]
    return np.mean(nbst)


#################################################################################
# draw pili anchors on hexsphere surface


import tmos
import tmos.base as base
from tmos import surface


# utils
def tonpy(v):
    return np.array([v.x, v.y, v.z])

def tov3d(arr):
    x,y,z = arr
    return base.Vector3d(x,y,z)

@command.defaultsave()
def draw_pili_anchors():
    mapped = compute_pili_anchors()
    ax = plt.gca()
    pili_anchors(ax, mapped, linekw={}, args=parameters.thisread())

def compute_pili_anchors():
    """ draw the pili distribution on a scatter plot. superimpose on one hex """
    args = parameters.thisread()
    assert(args.surface.shape == 'hexspheregrid')
    # construct c++ surface 
    hxsph = surface.HexSphereGrid(args.surface.sphere_radius)

    # read in local pili tracking data
    prs = readtrack.eventsettrack()
    cut = int(prs[0].size//10)
    for pr in prs:
        pr.slice = slice(cut, None, 1)

    # define function to extract attachment points
    def get_attach(pr):
        attidx = np.argwhere(pr['process'] == 'attach')
        xyz = np.column_stack(
                [pr['attach_x'][attidx], 
                pr['attach_y'][attidx], 
                pr['attach_z'][attidx]]
                )
        return xyz
    xyz_lss = list(map(get_attach, prs))
    xyzls = [row for block in xyz_lss for row in block]
    #xyz_block = np.concatenate(xyzls, axis=0)

    # pass to c++ structure to map onto origin hex
    def hexmap(pt):
        grid_pt = tonpy(hxsph.get_xyz( hxsph.get_rq( tov3d(pt) )))
        mapped_pt = pt - grid_pt
        return mapped_pt
    mapped_ls = list(map(hexmap, xyzls))
    mapped = np.stack(mapped_ls, axis=0)
    return mapped

def pili_anchors(ax, mapped, linekw={}, args={}):
    X, Y, Z = np.split(mapped, 3, axis=1)

    scatterstyle = {'alpha':0.3, 'c':'blue'}
    ax.scatter(X, Y, s=1, **scatterstyle)
    # draw the circle
    patchstyle = {'color':'k', 'fill':False, 'linewidth':1, 'alpha':0.5}
    sph = patches.Circle((0,0), args.surface.sphere_radius, **patchstyle)
    ax.add_patch(sph)
    ax.set_aspect('equal')


#################################################################################
# analyse mdevent.dat

## calculate the displacement associated with each pilus binding event
# see notebook/mdevent.py

def _mean(a, default=np.nan): 
    if not(a): return default
    return sum(a)/len(a)


# extract the retractions and corresponding displacements for each pilus displacements 

def _mdreorganize(md):
    pidx = md['pidx']
    pdata = {idx : [] for idx in pidx}
    for row in md._track:
        pdata[row['pidx']].append(row)
    for idx in pdata:
        pdata[idx] = np.stack(pdata[idx])
    return pdata

# construct a dataframe with columns
# pilus_index, no bindings, no retractions, total displacement, first retraction, last retraction

def pilistate(md, trackidx=0):
    norm = np.linalg.norm
    cols = [
        "trackidx", "pidx", "n_bindings", "n_retractions", "xydisplacement", "zdisplacement",
        "first_retraction", "last_retraction"
    ]
    dtypes = [
        int, int,  int, int, float, float, float, float
    ]
    mdpdata = _mdreorganize(md)
    n_pili = len(mdpdata)
    data = {name : np.full(n_pili, np.nan, dtype=dtype) for name, dtype in zip(cols, dtypes)}
    i = 0
    for pidx, block in mdpdata.items():
        process = block["process"]
        retraction = process == "retraction"
        _dat = {}
        _dat['pidx'] = pidx
        _dat['trackidx'] = trackidx 
        _dat['n_bindings'] = max(int(np.sum(process == 'detach')), 1)
        _dat['n_retractions'] = len(block[retraction])
        dxyz = np.column_stack([block["d_x"], block["d_y"], block["d_z"]])
        disp = np.sum(dxyz[:,:2],axis=0)
        _dat["xydisplacement"] = norm(disp)
        zdisp = np.sum(dxyz[:,2],axis=0)
        _dat["zdisplacement"] = norm(zdisp)
        if block[retraction].size > 0:
            _dat["first_retraction"] = block[retraction][0]['time']
            _dat["last_retraction"] = block[retraction][-1]['time']
        for k, v in _dat.items():
            data[k][i] = v
        i = i + 1
    return pd.DataFrame(data)

def ensemble_pilistate(mds):
    return pd.concat([pilistate(md, trackidx=i) for i, md in enumerate(mds)])

def mddistrib():
    mddata = readtrack.mdeventset()
    mddf = ensemble_pilistate(mddata)
    return mddf

def save_xydisp():
    mddf = mddistrib()
    np.save("xydisp.npy", mddf["xydisplacement"].to_numpy())

def mdsummary():
    sd = {}
    mddata = readtrack.mdeventset()
    mddf = ensemble_pilistate(mddata)
    xydisp = mddf["xydisplacement"].to_numpy()
    np.save("xydisp.npy", xydisp)

    disp = mddf["xydisplacement"]
    sd["pdisp"] = {}
    # todo save the displacement distribution
    sd["pdisp"]['mean'] = mddf["xydisplacement"].mean()
    sd["pdisp"]['median'] = mddf["xydisplacement"].median()
    sd["pilus"] = {}
    sd["pilus"]["n_bindings"] =  {"mean" : mddf["n_bindings"].mean()}
    sd["pilus"]["n_retractions"] =  {
        "mean" : mddf["n_retractions"].mean(), 
        "median" : mddf["n_retractions"].median()
        }
    return sd 


# ----------------------------------------------------------------
# sample the md data at interval dt and compute properties for the sampling interval
# properties: displacement, no. retracting pili, nbound pili, pili bound retraction fraction
def movestate(md, dt=0.1, mintime=100):
    """
    mintime is used to allow the system to reach equilibrium
    """
    md['process'] == 'retraction'
    time, nbound, pbrf = md['time'], md['nbound'], md['pbrf'] 
    dxyz = np.column_stack([md['d_x'],md['d_y'],md['d_z']])
    maxtime = time[-1]
    timepoint = np.arange(mintime+dt, maxtime, dt)
    # 
    result = {
        'nret':np.empty(timepoint.size,dtype=int), # count retractions
        'disp':np.empty_like(timepoint),  # sum displacement
        'nbound':np.empty_like(timepoint),
        'pbrf':np.empty_like(timepoint)
        }
    dxy = np.zeros(2)
    # init
    c = 0 
    i = np.searchsorted(time, mintime)
    partial_nb = []
    partial_pbrf = []
    t = time[i]
    print('first timepoint {} and index {}'.format(timepoint[0],i))
    for r_i, tnext in enumerate(timepoint):
        while True:
            if t < tnext:
                dxy += dxyz[i][:2]
                c += 1
                partial_nb.append(nbound[i])
                partial_pbrf.append(pbrf[i])
                i += 1
                t = time[i]
            else:
                d = np.linalg.norm(dxy)
                result['disp'][r_i] = d
                result['nret'][r_i] = c
                result['nbound'][r_i] = _mean(partial_nb)
                result['pbrf'][r_i] = _mean(partial_pbrf, default=0.0)
                dxy = np.zeros(2)
                c = 0 
                partial_nb = []
                partial_pbrf = []
                break
    # assert(len(nret) == timepoint.size)
    return result

# ----------------------------------------------------------------
# plotting
def plot_disp_nret(result):
    disp, nret = result['disp'], result['nret']
    fig, axes = plt.subplots(1,2, figsize=(10,5))
    ax1, ax2 = axes
    ax1.hist(nret, bins=max(nret))
    ax1.set_xlabel('count')
    ax2.hist(disp, bins=max(nret))
    ax2.set_xlabel('displacement')

@command.defaultsave()
def plot_nbound_pbrf(result, tr, par):
    fig, axes = plt.subplots(1,2, figsize=(10,5))
    ax1, ax2 = axes
    nbound = nbound_distrib(tr, par)
    plot_nbound(ax1, nbound, par)
    plot_pbrf(ax2, result, par)

def nbound_distrib(tr, par):
    assert(par.dt == 0.1)
    istart = np.searchsorted(tr['time'],par.mintime)
    time = tr['time'][istart:]
    tr.slice = slice(istart-1,None,None)
    xy = np.column_stack([tr['x'],tr['y']])
    d = np.linalg.norm(xy[1:] - xy[:-1],axis=1)
    tr.slice = slice(istart,None,None)
    nbound = tr['nbound']
    smallidx = d < par.dxthreshold
    par.smallidx = smallidx
    return nbound

def plot_nbound(ax1, nbound, par):
    hstyle = {'alpha': 0.5, 'bins':10}
    small_nb = nbound[par.smallidx]
    ax1.hist(nbound, label='all', **hstyle)
    ax1.hist(small_nb, label='not moving', **hstyle)
    ax1.set_xlabel('nbound')
    ax1.legend()
    
def plot_pbrf(ax2, result, par):
    # round retraction fraction?
    hstyle = {'alpha': 0.5, 'bins':10}
    disp = result['disp']
    smallidx = disp < par.dxthreshold
    ax2.hist(result['pbrf'], label='all', **hstyle)
    ax2.hist(result['pbrf'][smallidx], label='not moving', **hstyle)
    ax2.set_xlabel('pbrf')
    ax2.legend()




if __name__=='__main__':

    ff, thisargs = command.process(locals())
    ff(*thisargs)


