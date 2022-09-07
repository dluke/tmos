#!/usr/bin/env python3

import sys, os
import itertools
import pili # setup environment

import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

import wrinterval as wr

import command
from command import defaultsave

import plotutils
import readtrack




"""
Start with some simple plotting routines for pilitracking data file
"""

# paths
datfile = 'pilitracking.dat'

# global matplotlib configuration
mplstyle = os.path.join(os.getenv('PILI_ROOT'), 'mpl/tracking.mplstyle')
plt.style.use(mplstyle)

# ptdata is reorganised to {pidx : {column : [data] } }
# why don't I just organise it this way when it is read in ???
# usually we want data organised by pilus. By sometimes we may want data organised by time... 
# TODO read and convert data into this form in one function

def reorganise(ptdata):
    pidxmap = {}
    for item in ptdata:
        try:
            context, block = item
        except ValueError:
            print(item)
            raise
        time, npili = context
        pidx_arr = block['pidx']
        for pidx in pidx_arr:
            keys = block.dtype.names
            if pidx not in pidxmap:
                pidxmap[pidx] = {k : [] for k in list(keys) + ['time']}
            pidxmap[pidx]['time'].append(time)
        for key in block.dtype.names: 
            arr = block[key]
            for i, arr_v in enumerate(arr):
                pidxmap[pidx_arr[i]][key].append(arr_v)
    # convert to numpy arrays
    for pidx in pidxmap.keys():
        colmap = pidxmap[pidx]
        for key in colmap:
            colmap[key] = np.array(colmap[key])
    return pidxmap
    
def get_segment(tr, ptr, cand_ev, pidx, ptime):
    tstart, tend = ptime[0], ptime[-1]
    # TODO times not necessarily contiguous
    stidx = np.argwhere(tr['time'] == tstart)[0][0]
    seidx = np.argwhere(tr['time'] == tend)[0][0]
    trwhole = tr[stidx:seidx+1]
    if len(trwhole) != ptime.size:
        print('warning: pilus data not contiguous')
    # cut the segment to the first release event
    i_truerelease = np.logical_and(
        cand_ev['process'] == 'release', 
        cand_ev['trigger'] == 'release')
    release_times = cand_ev['time'][i_truerelease]
    release_time = release_times[0]
    release_idx = np.searchsorted(trwhole['time'], release_time, side='left')
    trseg = trwhole[:release_idx]
    #
    print('pilus {} slice ({},{})'.format(pidx, stidx, seidx))
    prowdata = []
    for tidx in range(stidx, stidx+release_idx):
        ptrstruct = ptr[tidx][1]
        row = ptrstruct[ptrstruct['pidx']==pidx]
        prowdata.append(row)
    return trseg, prowdata

# plot pleg and plenght so they are comparable
# also plot the difference on separate axis
@defaultsave(svg=True)
def length():
    ptdata = readtrack.piliset()[0]
    _length(ptdata)
    
def _length(ptdata):
    pilusdata = reorganise(ptdata)
    color = itertools.cycle(plotutils.colors)
    kwleq = {'linestyle':'-', 'alpha':0.6}
    kwlength = {'linestyle':'--', 'alpha':0.6}
    for pidx, block in pilusdata.items():
        thisc = next(color)
        transition_idx = np.nonzero(np.diff(block['isbound']))[0]
        transition_time = block['time'][transition_idx]
        transition_state= block['isbound'][1:][transition_idx]
        for time, state in zip(transition_time, transition_state):
            if state == 0:
                plt.axvline(time, color='green', linewidth=0.5, alpha=0.5, linestyle='--')
            if state == 1:
                plt.axvline(time, color='k', linewidth=0.5, alpha=0.5)
        plt.plot(block['time'], block['pleq'], color=thisc, **kwleq)
        plt.plot(block['time'], block['plength'], color=thisc, **kwlength)
    plt.xlabel('time (s)')
    plt.ylabel('length ($\mu m$)')
    
#
def anycol(col):
    ptdata = readtrack.piliset()[0]
    print(('extracting column {}'.format(col)))
    pidxmap = {}
    timeline = []
    for item in ptdata:
        context, block = item
        time, npili = context
        timeline.append(time)
        for pidx, colval in zip(block['pidx'], block[col]):
            if pidx not in pidxmap:
                pidxmap[pidx] = (list(),list())
            pidxmap[pidx][0].append(time)
            pidxmap[pidx][1].append(colval)
    # plot pidxmap
    print('plotting...')
    for pidx, data in list(pidxmap.items()):
        time, ydata = data
        plt.plot(time, ydata, label=str(pidx))

    plt.xlabel('time (s)')
    plt.ylabel(col)
    plt.legend()
    command.saveplt(col)


#
visible_ranges = [(0,99.0), (0.3,99.0),(1.0,99.0),(2.0,99.0)]
markers = ['x','^','v']



# plot total number of pili but also number of pili in Tala and Koch visible ranges
@defaultsave(autoshow=False,svg=True)
def npili():
    ptdata = readtrack.piliset()[0]
    _npili(ptdata)

def count_visible(ptdata, visible_ranges=visible_ranges):
    timeline = []
    npy = []
    npy_in_range = {vr: [] for vr in visible_ranges}
    bound_in_range = {vr: [] for vr in visible_ranges}
    taut_in_range = {vr: [] for vr in visible_ranges}
    i = 0
    for item in ptdata:
        context, block = item
        time, npili = context
        timeline.append(time)
        npy.append(npili)

        # if i > 20:
        #     break

        for vr in npy_in_range:
            npy_in_range[vr].append(0)
            bound_in_range[vr].append(0)
            taut_in_range[vr].append(0)

            # print(block['pleq'])

            for i, pleq in enumerate(block['pleq']):
                isbound = block['isbound'][i]
                istaut = block['istaut'][i]
                if pleq >= vr[0] and pleq < vr[1]:
                    npy_in_range[vr][-1] += 1
                    if isbound:
                        bound_in_range[vr][-1] += 1
                    if istaut:
                        taut_in_range[vr][-1] += 1
        # print(npy_in_range[vr][-1])
        i += 1

    return timeline, npy, npy_in_range, bound_in_range, taut_in_range

def _npili(ptdata):

    timeline, npy, npy_in_range, bound_in_range, taut_in_range = count_visible(ptdata)

    fig1 = plt.figure()
    marker = iter(markers)
    plt.plot(timeline, npy, marker=next(marker), label='all')
    npyrform = 'L > {:3.1f}'
    for vr in npy_in_range:
        plt.plot(timeline, npy_in_range[vr], 
                label=npyrform.format(vr[0]))
    plt.xlabel('time (s)')
    plt.ylabel('Number of Pili')
    plt.legend()

    mean_npili = np.mean(npy)
    form = 'Average number of pili in range {} is {:3.1f}/{:3.1f}'
    for vr, nlist in list(npy_in_range.items()):
        print((form.format(vr, np.mean(nlist), mean_npili)))

    fig_hist = plt.figure()
    fig_hist._plot_name = 'distrib'
    offset = iter([-0.2, 0.0, 0.2])
    glob_max = 0
    for vr in npy_in_range:
        #kw = {'density':True, 'label':npyrform.format(vr[0])}
        # +1 for inclusive range, +1 for rightmost bin edge
        local_max = max(npy_in_range[vr])
        glob_max = max(glob_max, local_max)
        binseq = np.array(list(range(local_max+2)))
        hist, bin_edges = np.histogram(npy_in_range[vr], 
                bins=binseq, density=True)
        plt.bar(binseq[:-1] + next(offset), hist, width=0.35, align='center',
                label=npyrform.format(vr[0]))

    plt.xlabel('Number of Pili')
    plt.ylabel('P')
    plt.xticks(list(range(glob_max+1)),list(range(glob_max+1)))
    plt.legend()

#####
# Some experimental use of statistical tests

# use the average number of pili to determine when the simulation reaches a steady state
import mkt
import scipy.ndimage
def init_sim_time():
    ptdata = readtrack.piliset(1)[0]
    _init_sim_time(ptdata)


# statistical tests are not satisfied that there is no trend in the data.
# however slope is very low.
# probably not a good idea to spend too much more time on this

def _init_sim_time(ptdata):
    timeline = []
    num_pili = []
    print('data length ', len(ptdata))
    # for item in ptdata[1000:None:40]:
    for item in ptdata:
        context, _ = item
        time, npili = context
        timeline.append(time)
        num_pili.append(npili)
    timeline = np.array(timeline)
    num_pili = np.array(num_pili)

    # what if we used pl_avg instead of number since its a continuous quantity?

    # smooth_num_pili = scipy.ndimage.gaussian_filter1d(num_pili.astype(float), 400)
    # plt.plot(timeline, num_pili, label='')
    # plt.plot(timeline, smooth_num_pili, label='')

    dt = np.diff(timeline)
    # quick weighted average
    # a_num = np.sum(dt*num_pili[1:])/np.sum(dt)
    print('computing M-K test')
    MK, m, c, p = mkt.test(timeline, smooth_num_pili, eps=1e-6, alpha=0.05, Ha='up')
    print('MK result: {}'.format(MK))
    print('slope: {}'.format(m))
    print('p-value: {}'.format(p))


######

# for debugging wandering pili attachment positions
@defaultsave()
def all_att():
    ptdata = readtrack.piliset()[0]
    _all_att(ptdata)
    
def _all_att(ptdata):
    xls = []
    yls = []
    for item in ptdata:
        context, block = item
        # time, npili = context
        attachedidx = np.nonzero(block['isbound'])[0]
        x = block['attach_x'][attachedidx]
        y = block['attach_y'][attachedidx]
        xls.append(x)
        yls.append(y)
    X = np.concatenate(xls,axis=0)
    Y = np.concatenate(yls,axis=0)
    plt.scatter(X,Y, alpha=0.1)
        
######

# extension/retraction/stalled times, taking into account Tala/Koch visible ranges

@defaultsave()
def motorstate():
    ptdata = readtrack.piliset()[0]
    _motorstate(ptdata)

def _motorstate(ptdata):
    # 
    pili = {}
    for item in ptdata:
        context, block = item;
        time, npili = context
        pidx = block['pidx']
        ext = iter(block['ext_motor'])
        ret = iter(block['ret_motor'])
        pleq = iter(block['pleq'])
        for idx in pidx:
            pili.setdefault(idx, {'time':[], 'ext_motor':[],'ret_motor':[],'pleq':[]})
            pili[idx]['time'].append(time)
            pili[idx]['ext_motor'].append(next(ext))
            pili[idx]['ret_motor'].append(next(ret))
            pili[idx]['pleq'].append(next(pleq))
    # 
    ext_times = []
    ret_times = []
    stall_times = []
    for item in list(pili.items()):
        idx, data = item
        time = data['time']
        ext_motor = np.array(data['ext_motor'],dtype=bool)
        ret_motor = np.array(data['ret_motor'],dtype=bool)
        stall = np.logical_and(~ext_motor,~ret_motor)
        def get_times(barr):
            # 
            tls = []
            diff = np.diff(barr)
            to_state = barr[1:]
            transition_index = np.nonzero(diff)[0]
            tsx = transition_index.size
            for i, tx in enumerate(transition_index):
                if to_state[tx] == 1:
                    n_tx = transition_index[i+1] if i+1 < tsx else -1
                    tdelta = time[n_tx] - time[tx]
                    tls.append(tdelta)
            return tls
        p_ret_times = get_times(ret_motor)
        p_ext_times = get_times(ext_motor)
        p_stall_times = get_times(stall)
        ext_times.extend(p_ext_times)
        ret_times.extend(p_ret_times)
        stall_times.extend(p_stall_times)
    #
    fig_statetime = plt.figure()
    fig_statetime._plot_name = 'extret'
    plotutils.kdeplot(np.asarray(ext_times), res=40, linekw={'label':'extension times'})
    plotutils.kdeplot(np.asarray(ret_times), res=40, 
            linekw={'label':'retraction times'}, xlims=(0, 40))
    plt.xlabel('time (s)')
    plt.ylabel('P')
    plt.ylim(ymin=0)
    plt.legend()

    fig_stalltime = plt.figure()
    fig_stalltime._plot_name = 'stall'
    kw = {'alpha': 0.5, 'density':True}
    plt.hist(stall_times, label='stall times', **kw)
    plt.xlabel('time (s)')
    plt.ylabel('P')
    plt.ylim(ymin=0)
    plt.legend()

    # pili extension/retraction/stalled ratio in visible range
    # NOTE: If this runs too slow, rearrange the loops
    vr_ratio = []
    for vr in [(0.0,99.0)] + visible_ranges:
        ret_count = 0
        ext_count = 0
        stall_count = 0
        for item in ptdata:
            context, block = item
            time, npili = context
            for i, pidx in enumerate(block['pidx']):
                pleq = block['pleq'][i] 
                if pleq > vr[0]:
                    if block['ext_motor'][i] == 1:
                        ext_count += 1
                    elif block['ret_motor'][i] == 1:
                        ret_count += 1
                    else:
                        stall_count += 1

        total = ret_count + ext_count + stall_count
        ratio = (float(ret_count)/total, float(ext_count)/total, float(stall_count)/total)
        vr_ratio.append(ratio)
        form = "For visible range L > {:3.1f} (retraction : extension : stalled) = {}"
        print(form.format(vr[0], ratio))


# for pivot events, see notebook/pivot.py
# summary statistics? 
# min(theta), max(theta) over pili and timesteps, where theta is the angle of pili relative to -e_x

# --- use tmos.base to reconstruct the body and pili geometry
import tmos.base as base
from tmos.base import Vector3d
from tmos.base import Frame

e_z = np.array([0,0,1])

def _construct_frame(trdata):
    ax = Vector3d(trdata['ax_x'], trdata['ax_y'], trdata['ax_z'])
    ax_e1 = Vector3d(trdata['rax_x'], trdata['rax_y'], trdata['rax_z'])
    ax_e2 = ax.cross(ax_e1) # e3 x e1 = e2
    origin = Vector3d(
        (trdata['x']+trdata['trail_x'])/2, 
        (trdata['y']+trdata['trail_y'])/2, 
        (trdata['z']+trdata['trail_z'])/2
        )
    return base.Frame(origin, ax_e1, ax_e2, ax)

def v3d_to_npy(v):
    return np.array([v.x, v.y, v.z])

def _construct_anchor(frame, rowdata):
    anchor = Vector3d(rowdata['anchor_x'], rowdata['anchor_y'], rowdata['anchor_z'])
    lab_anchor = v3d_to_npy(frame.to_lab(anchor))
    return lab_anchor

# this is the attached axis (attachment - anchor)
def _construct_axis(lab_anchor, rowdata):
    lab_attach = np.array([rowdata['attach_x'], rowdata['attach_y'], rowdata['attach_z']]).ravel()
    axis = lab_attach - lab_anchor
    axis_hat = axis/np.linalg.norm(axis)
    return axis_hat

def _pili_reconstruction(ptrs, trs):
    # create the following lists of arrays 
    print("reconstructing frames")
    frames = [ [_construct_frame(row) for row in tr[:]] for tr in tqdm(trs)]
    # min_z_theta = [np.zeros(tr.size, dtype='f8') for tr in trs]
    # max_z_theta = [np.zeros(tr.size, dtype='f8') for tr in trs]
    min_z_theta = []
    max_z_theta = []
    #
    for i, ptdata in enumerate(ptrs):
        frame_list = frames[i]
        min_ = []
        max_ = []
        for idx, item in enumerate(ptdata):
            context, block = item
            time, n = context
            if n == 0:
                min_.append(np.nan)
                max_.append(np.nan)
            else:
                frame = frame_list[idx]
                # TODO, methods updated
                anchors = [_construct_anchor(frame, rowdata) for rowdata in block[:]]
                axes = [_construct_axis(frame, rowdata) for rowdata in block[:]]
                #
                v_theta = [np.arccos(-axis[2]) for axis in axes]
                min_.append(np.min(v_theta))
                max_.append(np.max(v_theta))
        min_z_theta.append(np.array(min_))
        max_z_theta.append(np.array(max_))
    return min_z_theta, max_z_theta

def _check_configuration(ptrs, trs, dd_pivot, dd_action):
    # note we are using pivot

    # use post_processed data -- bound pili only
    summary_stats = ["max_theta", "min_theta", "mean_npili"]
    dd_stats = []
    # 
    for s_tup, e_tup in zip(dd_pivot['s_idx'], dd_pivot['e_idx']):
        # retrieve the start and end indexes of the pivots
        tr_i, s_idx = s_tup
        tr_i, e_idx = e_tup
        ptdata = ptrs[tr_i]
        tr = trs[tr_i]
        # stats for action
        stats = {k: [] for k in summary_stats}
        # for pivot action
        for start_i, end_i in zip(s_idx, e_idx):
            pidx = []
            theta = []
            pcount = []
            # for step in pivot 
            for i in range(start_i, end_i+1):
                trdata = tr[:][i]
                context, block = ptdata[i]
                time, npili = context
                pcount.append(npili)
                for rowdata in block:
                    # reconstruct axis 
                    frame = _construct_frame(trdata)
                    axis = _construct_axis(frame, rowdata)
                    # v_theta = np.arccos(np.dot(axis, -e_z))
                    v_theta = np.arccos(-axis[2])
                    # print('axis', axis, v_theta)
                    r_pidx = rowdata['pidx']
                    pidx.append(r_pidx)
                    theta.append(v_theta)
            # summary statistics for this action
            stats["mean_npili"].append( np.mean(pcount) )
            stats["max_theta"].append( np.max(theta) )
            stats["min_theta"].append( np.min(theta) )
        # collect stats for separately for each track
        dd_stats.append(stats)
    return dd_stats


if __name__=='__main__':

    #anycol('pleq')    

    ff, thisargs = command.process(locals())
    ff(*thisargs)

