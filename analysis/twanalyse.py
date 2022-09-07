#!/usr/bin/env python3

"""
 analysis routines involving simple average / variance

 This is the top level analysis module so don't import it in any other analysis module
"""

import sys, os
import json
join = lambda *x: os.path.abspath(os.path.join(*x))
import numpy as np
import pandas as pd
import scipy.stats
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

import collections
from collections import OrderedDict
from glob import glob
import itertools
from tqdm import tqdm

# 
import pili
import wrinterval as wr
import txtdata
import parameters

# has pdir and datdir for defining where to save data
import command, filesystem 
from command import defaultsave
import stats
import plotutils


# tmos analysis code
import _fj
import twutils
import readtrack
import astrack
import eventanalyse
import pilush 

##################################
# summary method

def filter_failed(trs):
    # read the expected simulation end time 
    simtime = parameters.thisread().system.simtime
    def is_failed(tr):
        # harsh criteria, throw away the whole track if simulation did not complete
        return tr['time'][-1] < simtime
    ok = [tr for tr in trs if not is_failed(tr)]
    total = len(trs)
    print("Filtered out {}/{} simulation tracks".format(total-len(ok),total))
    return ok

# print out summary of useful statistical parameters
def summary():
    trs = readtrack.trackset()
    return _summary(trs)

@stats.keep
def _summary(trs):

    # preprocessing
    tr = trs[0]
    # is_walking = tr['trail_z'][0] >  1.0
    # if is_walking:
    #     print("identified this simulated data as walking and adjusting the leading/trailing poles accordingly")
    # if is_walking or not is_walking:
    trs = [tr.extend_projected_axis_by_radius() for tr in trs]
    trs = [tr.cut_start(10) for tr in trs]
    #

    sd = {}
    if len(trs) == 0:
        sd['failed'] = True
        sd['failed_condition'] = "all simulations crashed"
        return sd

    duration = [tr.get_duration() for tr in trs]
    sd['simtime'] = {}
    sd['simtime']["total"] = np.sum(duration)
    sd['simtime']["each"] = list(duration)
    sd['vel_30'] = _avgvel(trs, step=300)
    sd['vel'] = stats.col_stats(trs, _inst_vel)
    if 'l_total' in trs[0].get_dtype().names:
        sd['l_total'] = {}
        sd['l_total']['mean'] = np.mean([np.mean(tr['l_total']) for tr in trs])

    
    eventanalyse.lifetime()
    sd['nbound'] = stats.col_stats(trs, lambda tr: tr['nbound'])
    sd['ntaut'] = stats.col_stats(trs, lambda tr: tr['ntaut'])
    sd['npili'] = stats.col_stats(trs, lambda tr: tr['npili'])

    # fanjin linearization procedure to make simulation velocities comparable
    ltrs = [_fj.linearize(tr) for tr in trs]
    linearsteps = [ltr.get_nsteps() for ltr in ltrs]
    sd['linearsteps'] = {}
    sd['linearsteps']['total'] = int(np.sum(linearsteps))
    sd['linearsteps']['each'] = list(linearsteps)
    sd.update(out_of_plane(ltrs))

    #
    step_threshold = 10
    if all([len(ltr.step_idx) < step_threshold for ltr in ltrs]):
        sd['failed'] = True
        sd['failed_condition'] = "step_condition"
        return sd

    sd.update(observables(ltrs))
    qasd = _qaparams(ltrs) # stats.keep
    sd.update(qasd)
    sd.update(vel_similarity(ltrs))

    # linearize again  but with d_step = 1.0.  This is to obtain q,a on lengthscales similar to bacteria size
    # sd.update(relinearize(ltrs, 'double_ltrs', 0.24))
    # sd.update(relinearize(ltrs, 'quad_ltrs', 0.48))
    # sd.update(relinearize(ltrs, 'cell_ltrs', 1.0))
    sd.update(kmsd_summary(ltrs))
    # sd.update(__vel_corr(ltrs, 10))

    # mdevent data
    # TODO gracefully check this data exists
    sd.update(eventanalyse.mdsummary())

    # v_tau = stats.load_or_compute('velocity_decay_time', lambda: _vel_corr(trs) )
    # set_resolution(trs, 2*v_tau)
    twutils.print_dict(sd)
    return sd

def lightweight_summary():
    trs = readtrack.trackset()
    return _lightweight_summary(trs, dstep=0.03, localjson="light.json")

def sim_summary():
    trs = readtrack.trackset()

    sd = {}
    
    def deviation(tr):
        dx = tr.get_dx()
        body = tr.get_body_projection()[:-1]
        theta = angle(dx, body)
        # theta is nan if projected body length is zero
        theta = np.nan_to_num(theta)
        return theta
    deviation_data=  np.concatenate([deviation(tr) for tr in trs])

    sd['deviation'] = {}
    sd['deviation']['var'] = np.var(deviation_data)
    sd['deviation']['std'] = np.std(deviation_data)

    lvel = np.concatenate([np.linalg.norm(tr.get_dx(),axis=1)/0.1 for tr in trs])
    sd['lvel'] = {}
    sd['lvel']['mean'] = np.mean(lvel)

    disp_list = [tr.get_dx() for tr in trs]
    sim_ut = []
    sim_up = []
    for svel in disp_list:
        sim_ut.append(svel[1:])
        sim_up.append(svel[:-1])
    totalut = np.concatenate(sim_ut)
    totalup = np.concatenate(sim_up)

    qhat, ahat = mle(totalut, totalup)
    info_qq, info_aa = fisher_information(
        totalut, totalup, qhat, ahat)
    err_q = 1.96/np.sqrt(info_qq)
    err_a = 1.96/np.sqrt(info_aa)
    sd['qhat'] = {'estimate': qhat, 'err': err_q, 'fisher_info':info_qq}
    sd['ahat'] = {'estimate': ahat, 'err': err_a, 'fisher_info':info_aa}

    # long distance coarse graining
    # sd['long'] = {}
    # dstep = 3.0
    # ltrs = [_fj.linearize(tr, dstep) for tr in trs]
    # qasd = _qaparams(ltrs) # stats.keep

    with open("no_coarse.json", 'w') as f:
        json.dump(sd, f, indent='\t')





def _lightweight_summary(trs, dstep=0.03, localjson="light.json"):

    # preprocessing
    tr = trs[0]
    trs = [tr.extend_projected_axis_by_radius() for tr in trs]
    trs = [tr.cut_start(10) for tr in trs]
    #

    sd = {}
    if len(trs) == 0:
        sd['failed'] = True
        sd['failed_condition'] = "all simulations crashed"
        return sd

    duration = [tr.get_duration() for tr in trs]
    sd['simtime'] = {}
    sd['simtime']["total"] = np.sum(duration)
    sd['simtime']["each"] = list(duration)

    # fanjin linearization procedure to make simulation velocities comparable
    sd["dstep"]= dstep
    ltrs = [_fj.linearize(tr, step_d=dstep) for tr in trs]
    linearsteps = [ltr.get_nsteps() for ltr in ltrs]
    sd['linearsteps'] = {}
    sd['linearsteps'] = {}
    sd['linearsteps']['total'] = int(np.sum(linearsteps))
    sd['linearsteps']['each'] = list(linearsteps)

    #
    step_threshold = 10
    if all([len(ltr.step_idx) < step_threshold for ltr in ltrs]):
        sd['failed'] = True
        sd['failed_condition'] = "step_condition"
        return sd

    sd.update(observables(ltrs))
    qasd = _qaparams(ltrs) 
    sd.update(qasd)
    
    with open(local_json, 'w') as f:
        json.dump(sd, f, indent='\t')

    
def relinearize(trs, name, step):
    step_threshold = 10
    sd = {}
    sd[name] = {}
    _ltrs = [_fj.linearize(tr, step_d=step) for tr in trs]
    # check size
    lnstep = [cltr.get_nsteps() for cltr in _ltrs]
    totalstep = int(np.sum(lnstep))
    sd[name]['linearsteps'] = totalstep
    def _lvel(ltrs):
        return stats.col_stats(ltrs, 
            lambda tr: tr.get_step_speed(), 
            weight_f=lambda tr: tr.get_step_dt(), each=False)

    if totalstep < step_threshold:
        sd[name]['failed'] = True
        return sd
    elif totalstep >= step_threshold:
        sd[name]['lvel'] = _lvel(_ltrs)
        cell_qasd = _qaparams(_ltrs)
        sd[name].update(cell_qasd)
    return sd
        

def observables(ltrs):
    sd = {}
    
    pdata = allpolar(ltrs)
    sd['deviation'] = {}
    sd['deviation']['var'] = np.var(pdata.deviation)
    sd['deviation']['std'] = np.std(pdata.deviation)

    sd['lvel'] = stats.col_stats(ltrs,
        lambda tr: tr.get_step_speed(), 
        weight_f=lambda tr: tr.get_step_dt())
    return sd

def kmsd_summary(ltrs):
    sd = {}
    sd['kmsd'] = {}
    _kmsd = np.array([kmsd(tr) for tr in ltrs])
    duration = np.array([tr.size for tr in ltrs])
    mean, var = twutils.wmeanvar(_kmsd, ws=duration)
    sd['kmsd']['mean'] = mean
    sd['kmsd']['var'] = var
    return sd

def kmsd(tr):
    _ptr = kmsd_one(tr)
    if _ptr is None:
        return np.nan
    p, cov, scaling, msd_n = _ptr
    kmsd, _intercept = p
    return kmsd

# read tracks and linearise
# convenient method
def linearise(trs): 
    return [_fj.linear]

### velocity profile similarity metrics
# NOTE: these metric actually require a reference dataset!
# that means we have to choose what to compare to in advance
# if we want to choose at a later date, all we can do is save the
# step velocities and distributions to some additional meta data file
# which we try to avoid

def load_simulation_target():
    rundir = join(pili.root, '../run')
    path = join(rundir, "825bd8f/target/t0/data/")
    return readtrack.trackset(ddir=path)

def get_walking_target_dir():
    rundir = join(pili.root, '../run')
    path = join(rundir, "825bd8f/target/t2/")
    return path

def load_walking_target():
    path = join(get_walking_target_dir(), 'data/')
    return readtrack.trackset(ddir=path)

@stats.keep
def run_vel_similarity(fanjin=False):
    trs = readtrack.trackset()
    ltrs = [_fj.linearize(tr) for tr in trs]
    return vel_similarity(ltrs, fanjin)

def vel_similarity(ltrs, fanjin=True):
    # 1. load speed distributions for fanjin data subsets
    # 2. compute two similarity metrics for each subset
    sd = {}
    sim_vel = np.concatenate([ltr.get_step_speed() for ltr in ltrs])
    # save this simulation velocity distribution as numpy locally 
    np.save("lvel.npy", sim_vel)
    
    def lvelscore(sim_vel, ref_vel):
        simscore = {}
        v1  = sim_vel- np.mean(sim_vel)
        v2  = ref_vel - np.mean(ref_vel)
        simscore["chi"] = chisquare(v1, v2)
        simscore["ks_statistic"] = ks_stat(sim_vel, ref_vel)
        # simscore["ks_nomean"] = ks_stat(v1, v2)
        return simscore

    if fanjin:
        sd["fanjin"] = {}
        distrib = _fj.load_subset_speed()
        for name, ref_vel in distrib.items():
            sd["fanjin"][name] = lvelscore(sim_vel, ref_vel)


    # same for simulation target
    sd["sim"] = {}
    #
    # target = load_simulation_target()
    # refltrs = [_fj.linearize(tr) for tr in target]
    # ref_vel = np.concatenate([_ltr.get_step_speed() for _ltr in refltrs])
    # sd["sim"]["t0"] = lvelscore(sim_vel, ref_vel)

    # # same for walking target
    # target = load_walking_target()
    # refltrs = [_fj.linearize(tr) for tr in target]
    # ref_vel = np.concatenate([_ltr.get_step_speed() for _ltr in refltrs])
    # sd["sim"]["t2"] = lvelscore(sim_vel, ref_vel)

    targets = ["t0", "t2", "t6", "t7"]
    for t in targets:
        path = join(pili.root, "../run/825bd8f/target/", t)
        trs = readtrack.trackset(ddir=join(path, 'data'))
        refltrs = [_fj.linearize(tr) for tr in trs]
        ref_vel = np.concatenate([_ltr.get_step_speed() for _ltr in refltrs])
        sd["sim"][t] = lvelscore(sim_vel, ref_vel)

    return sd

# similarity using kde and weighted squared difference
def chisquare(v1, v2, res=100):
    xn1, xm1 =  np.quantile(v1, 0.025), np.quantile(v1, 0.975)
    xn2, xm2 =  np.quantile(v2, 0.025), np.quantile(v2, 0.975)
    xn, xm = min(xn1, xn2), max(xm1, xm2)
    mspace = np.linspace(xn, xm, res)
    method = "scott"
    kde1 = scipy.stats.gaussian_kde(v1, bw_method=method)
    kde2 = scipy.stats.gaussian_kde(v2, bw_method=method)
    pde1 = kde1.evaluate(mspace)
    pde2 = kde2.evaluate(mspace)
    chisquared = np.sum((pde2 - pde1)**2/(pde1 + pde2))
    return np.sqrt(chisquared)

def ks_stat(v1, v2):
    # doesn't really make sense to subtract means for the ks_statistic
    # but changing it now means rerunning a lot of summary statistics
    # TODO: compare results with these lines uncommented
    # v1  = v1 - np.mean(v1)
    # v2  = v2 - np.mean(v2)
    ks_statistic, pvalue = scipy.stats.ks_2samp(v1, v2)
    return ks_statistic

#################################################################################
### persistence and activity parameters

# call with step velocity u[1:], u[:-1]
def q_estimator(u_t,u_p, sample=None):
    if sample is None:
        sample = np.array(range(u_t.shape[0]))
    return np.sum( (u_t*u_p).sum(axis=1)[sample] )/np.sum( (u_p*u_p).sum(axis=1)[sample] )

def mle(u_t,u_p, sample=None):
    if sample is None:
        sample = np.array(range(u_t.shape[0]),dtype=int)
    qhat = q_estimator(u_t,u_p, sample)
    upart = (u_t - qhat*u_p)[sample]
    # sum over both dimensions here is ok
    ahat = np.sqrt(np.sum(upart*upart)/(2*(u_t.shape[0])))
    return qhat, ahat

def fisher_information(u_t, u_p, qhat, ahat):
    n = u_t.shape[0]
    qq = np.sum( (u_p*u_p).sum(1) )/ahat**2
    _num = u_t - qhat*u_p
    aa = 3*np.sum(_num*_num)/ahat**4 - 2*(n-1)/ahat**2
    return qq, aa
# observed_qq, observed_aa = fisher_information(u[1:],u[:-1],qhat,ahat)

@stats.keep
def qaparams():
    trs = readtrack.trackset()
    ltrs = [_fj.linearize(tr) for tr in trs]
    return _qaparams(ltrs)

def _qaparams(ltrs, vquantile=0.99):
    sd = {}
    step_vel = [ltr.get_step_velocity() for ltr in ltrs]
    sim_ut = []
    sim_up = []
    for svel in step_vel:
        sim_ut.append(svel[1:])
        sim_up.append(svel[:-1])
    totalut = np.concatenate(sim_ut)
    totalup = np.concatenate(sim_up)
    # velocity threshold is needed for fanjin data since 
    # there are some huge velocities (probably errors) and these effect the estimate significantly
    th = None
    if vquantile:
        # note: not very efficient
        norm = np.linalg.norm
        s_t = norm(totalut,axis=1)
        s_p = norm(totalup,axis=1)
        th = np.quantile(s_t, vquantile)
        select = np.logical_and(s_t < th, s_p < th)
        totalut = totalut[select]
        totalup = totalup[select]
        sd['mle_excluded'] = np.sum(~select)/select.size
    sd['mle_vquantile'] = vquantile
    sd['mle_vthreshold'] = th

    qhat, ahat = mle(totalut, totalup)
    info_qq, info_aa = fisher_information(
        totalut, totalup, qhat, ahat)
    err_q = 1.96/np.sqrt(info_qq)
    err_a = 1.96/np.sqrt(info_aa)
    sd['qhat'] = {'estimate': qhat, 'err': err_q, 'fisher_info':info_qq}
    sd['ahat'] = {'estimate': ahat, 'err': err_a, 'fisher_info':info_aa}
    return sd


#################################################################################
### Pivot actions
# action frequency
def freq(trs, dd):
    fq = sum([len(actions) for actions in dd['dx_all']])/sum(
        [tr['time'][-1]-tr['time'][0] for tr in  trs])
    return fq

def _pivot_summary_stats(trs, dd_pivot, dd_action):
    sd = {}
    sd["pivot_freq"] = freq(trs, dd_pivot)
    sd["median_pivot_duration"] = np.median(np.concatenate(dd_pivot['dt_all']))
    pivot_trail_dx = np.concatenate(dd_pivot['dx_trail_all'])
    sd["median_pivot_trail_dx"] = np.median(pivot_trail_dx)
    sd["mean_pivot_trail_dx"] = np.mean(pivot_trail_dx)
    #
    sd["median_action_duration"] = np.median(np.concatenate(dd_action['dt_all']))
    action_trail_dx = np.concatenate(dd_action['dx_trail_all'])
    sd["median_action_trail_dx"] = np.median(action_trail_dx)
    sd["mean_action_trail_dx"] = np.mean(action_trail_dx)
    return sd

@stats.keep
def pivot_summary_stats():
    trs = readtrack.trackset()
    dd_pivot, dd_action = pivot_action('all', trs)
    sd = _pivot_summary_stats(trs, dd_pivot, dd_action)
    twutils.print_dict(sd)
    return sd


#################################################################################
# Look for "Pivot" action in tracking data

# first it is useful to check the width and length
def body_dimensions(trs):
    basis_time = [tr['time'] for tr in trs]
    widths = [tr['width'] for tr in trs]
    lengths = [tr['length'] for tr in trs]
    # compute mean width/length distributions
    dd = {}
    dd['mean_width'] = np.array([np.mean(w) for w in widths])
    dd['mean_length'] = np.array([np.mean(l) for l in lengths])
    # linear regression gradient distribution
    import numpy.polynomial.polynomial as poly
    dd['w_gradient'] = np.array([poly.polyfit(time, width, 1)[1] for time, width in zip(basis_time, widths)])
    dd['l_gradient'] = np.array([poly.polyfit(time, length, 1)[1] for time, length in zip(basis_time, lengths)])
    return dd


def pivot_action(tr_idx, trs, max_delta=200):
    """
    max_delta is the max number of simulation steps that we compute displacements over
    hence max_delta = 200 means that the longest pivot can be 20 seconds. 
    """
    if tr_idx == 'all':
        tr_idx = list(range(len(trs)))
    def _evaluate(tr, start_idx, end_idx):
        dx = np.sqrt( (tr['x'][start_idx] - tr['x'][end_idx])**2 
            + (tr['y'][start_idx] - tr['y'][end_idx])**2 )
        dx_trail = np.sqrt( (tr['trail_x'][start_idx] - tr['trail_x'][end_idx])**2 
            + (tr['trail_y'][start_idx] - tr['trail_y'][end_idx])**2 )
        if 'length' in tr._track.dtype.fields:
            dl = tr['length'][end_idx] - tr['length'][start_idx] 
        else:
            # assume simulation body length does not change
            dl = np.zeros(start_idx.size)
        dt = 0.1 * (end_idx-start_idx)
        return dx, dx_trail, dl, dt
    items = ['s_idx', 'e_idx', 'dx_all', 'dx_trail_all', 'dl_all', 'dt_all']
    dd_pivot = {k : [] for k in items}
    dd_action = {k : [] for k in items}
    print("Computing pivot actions.")
    for i, tr in zip(tr_idx, trs):
        pivot_start_idx, pivot_end_idx, action_start_idx, action_end_idx = _pivot_action(tr, max_delta)
        dx, dx_trail, dl, dt = _evaluate(tr, pivot_start_idx, pivot_end_idx)
        dd_pivot['s_idx'].append((i, pivot_start_idx))
        dd_pivot['e_idx'].append((i, pivot_end_idx))
        dd_pivot['dx_all'].append(dx)
        dd_pivot['dx_trail_all'].append(dx_trail)
        dd_pivot['dl_all'].append(dl)
        dd_pivot['dt_all'].append(dt)
        dx, dx_trail, dl, dt = _evaluate(tr, action_start_idx, action_end_idx)
        dd_action['s_idx'].append((i, action_start_idx))
        dd_action['e_idx'].append((i, action_end_idx))
        dd_action['dx_all'].append(dx)
        dd_action['dx_trail_all'].append(dx_trail)
        dd_action['dl_all'].append(dl)
        dd_action['dt_all'].append(dt)
    return dd_pivot, dd_action


# look for displacements of the trailing pole relative to the leading pole
# we choose to do this on unlinearised data
# TODO clean up and separate the definitons of "pivots" and "actions" 
d_err = 0.03
d_x = 2 * d_err
d_step = 0.12
d_trail_x = d_step + d_x
def _pivot_action(tr, max_delta=None):
    time = tr._track['time']
    x, y = tr._track['x'], tr._track['y']
    trail_x, trail_y = tr._track['trail_x'], tr._track['trail_y']
    #
    # compute the largest number of steps that the leading pole does not displace by more than d_x
    if max_delta is None:
        step_idx, step_displacement = _fj.step_to_distance(tr, step_d=d_x)
        max_delta = np.diff(step_idx).max()
    # use this value to compute all possible displacements
    delta_x = np.zeros((max_delta, time.size))
    delta_trail_x = np.zeros((max_delta, time.size))
    # print("computing all displacements up to \Delta steps = {}".format(max_delta))
    for d_i in range(1, max_delta):
        delta_x[d_i][d_i:] = np.sqrt( (x[d_i:] - x[:-d_i])**2 + (y[d_i:] - y[:-d_i])**2 )
        delta_trail_x[d_i][d_i:] = np.sqrt( (trail_x[d_i:] - trail_x[:-d_i])**2 + (trail_y[d_i:] - trail_y[:-d_i])**2 )
    # true/false arrays
    x_c = delta_x < d_x
    trailx_c = delta_trail_x > d_trail_x
    #
    pivot_sq = np.logical_and(x_c, trailx_c)
    # still need to extract the individual actions, taking the longest displacment in each case
    pivot = np.zeros(time.size)
    for i, row in enumerate(pivot_sq):
        for j in np.nonzero(row == True)[0]:
            pivot[j-i:j] = 1
    pdiff = np.diff(np.insert(np.append(pivot,0),0,0))
    pivot_start_idx = np.nonzero(pdiff == 1)[0]
    pivot_end_idx = np.nonzero(pdiff == -1)[0]
    assert(pivot_start_idx.size == pivot_end_idx.size)
    # now instead find the fastest action which still gives a displacement (d_trail_x - d_x) > d_step
    # TODO clean this up, find the fastest actions from the original distances
    start_idx = list(pivot_start_idx)
    end_idx = list(pivot_end_idx)
    action_start_idx = []
    action_end_idx = []
    for i in range(1, max_delta):
        dx_chunks = [delta_x[i][start_i:end_i] for start_i, end_i in zip(start_idx, end_idx)]
        dxt_chunks = [delta_trail_x[i][start_i:end_i] for start_i, end_i in zip(start_idx, end_idx)]
        to_remove = []
        for chk_i, (xchunk, xtchunk) in enumerate(zip(dx_chunks, dxt_chunks)):
            piv_idx = np.nonzero( (xtchunk - xchunk) > d_step )[0]
            # tmp
            if piv_idx.size > 0:
                idx = piv_idx[0]
                action_start_idx.append(pivot_start_idx[chk_i]+idx-i)
                action_end_idx.append(pivot_start_idx[chk_i]+idx)
                to_remove.append(chk_i)
        for chk_i in to_remove[::-1]:
            del start_idx[chk_i]
            del end_idx[chk_i]
    action_start_idx = np.array(action_start_idx, dtype=int)
    action_end_idx = np.array(action_end_idx, dtype=int)
    return pivot_start_idx, pivot_end_idx, action_start_idx, action_end_idx


########################################################
# analysis methods

def _concat(ll):
    return [x for l in ll for x in l]

def try_out_of_plane():
    trs = readtrack.trackset()
    ltrs = [_fj.linearize(tr) for tr in trs]
    sd = out_of_plane(ltrs)
    twutils.print_dict(sd)

def out_of_plane(ltrs):
    # compute some stats relevant to walking 
    # i.e. body angle, aspect ratio, pole travel 
    # compute aspect ratio variance 
    # fanjin ellipse tracking is on the tips of the poles where as simulated tracking is sphere centers
    # -- fanjin style tracking pretty much relies on the face that exactly vertical orientation is rare
    # in addition to aspect ratio, record the variance of (projected length / max length)
    #  -- theta angle can be computed from projected length and real length (but we don't know the real length)
    sd = {}
    is_simulated = ltrs[0].source == "simulated"
    data = {name:[] for name in 
        ["ztheta", "aspect", "pole_travel_score", "lratio"]
    }
    for ltr in ltrs:
        head = ltr.get_head()
        tail = ltr.get_trail()

        body_vec = head - tail
        body_ax = body_vec/np.linalg.norm(body_vec,axis=1)[:,np.newaxis]
        plane_project = np.sqrt(body_ax[:,0]**2 + body_ax[:,1]**2)
        _theta = np.abs(np.arctan2(body_ax[:,2], plane_project))
        # CAREFUL -- we should really modify the experimental data so emulate tracking of the same part of the cell 
        #  but it leads to some nasty edge cases
        if is_simulated:
            # we know maxl = 3.0 for simulated data by apply the same procedure anyway
            maxl = np.quantile(plane_project+1.0, 0.99) 
            data["ztheta"].append(_theta) #  only makes sense for simulated data
            data["aspect"].append( (plane_project+1.0)/1.0)
            data["lratio"].append( (plane_project+1.0)/maxl )
        else:
            maxl = np.quantile(ltr['length'], 0.99) # we never use the actual maximum of experimental data, always a quantile like .99 or .95
            width = ltr['width']
            data["aspect"].append( plane_project/width )
            data["lratio"].append( ltr['length'][ltr.step_idx]/maxl )
        data["pole_travel_score"] = _fj.pole_travel_score(ltr)

    def _stats(arr):
        return {"mean": np.mean(arr), "var": np.var(arr)}
    sd["ztheta"] = _stats(_concat(data["ztheta"]))
    sd["aspect"] = _stats(_concat(data["aspect"]))
    sd["lratio"] = _stats(_concat(data["lratio"]))
    sd["pole_travel_score"] = np.mean(data["pole_travel_score"])
    return sd

def get_candidate_distributions():
    candidate = 2924
    candidate_track = _fj.lintrackload([2924])[0]
    pdata = allpolar([candidate_track])
    velocities = candidate_track.get_step_velocity()
    return velocities, pdata.deviation

def get_linearised_data(ddir=None):
    trs = readtrack.trackset(ddir=ddir)
    ltrs = [_fj.linearize(tr) for tr in trs]
    return ltrs

def velocity_profiles(linearised=True):
    _velocity_profiles(get_linearised_data(), linearised)

# save the linearised velocity profiles
def _velocity_profiles(ltrs, linearised=True):
    prep='linearised' if linearised else '' 
    vdir = 'plots/vel/'
    filesystem.safemkdir(vdir)
    rule = '_'.join([prep, 'vel_{:04d}.png'])
    fig, ax = plt.subplots(1, 1, figsize=(100,6))
    for i, ltr in enumerate(ltrs):
        ax.set_xlabel('time')
        ax.set_ylabel('$\mu m/s$')
        if linearised:
            dt = ltr.get_step_dt()
            vel = ltr.get_step_velocity()
        else:
            dt = ltr['time'][1:] - ltr['time'][:-1]
            vel = ltr.get_head_v()
        speed = np.linalg.norm(vel,axis=1)
        ax.plot(np.cumsum(dt), speed)
        out = os.path.join(vdir, rule.format(i))
        print('save velocity profile to ', out)
        plt.savefig(out)
        plt.clf()

# for printing a single variable to std out
# use with seqwalk.py utility
def stat(var):
    value = stats.load().get(var, None)
    twutils.print_dict({var:value})

def anylims(col):
    assert col in txtdata.keys
    track = readtrack.Track()
    coldata = track[col]
    return np.min(coldata), np.max(coldata)

@defaultsave()
def anycol(col):
    assert col in txtdata.keys
    _anycol(readtrack.Track(), col)

_NUMERIC_KINDS = set('buifc')
def _anycol(track, col):
    data = track[col]
    plt.plot(track['time'], data)
    plt.xlabel('time')
    plt.ylabel(col)

    if data.dtype.kind in _NUMERIC_KINDS:
        print("mean {}".format(col), np.mean(data))
    plt.gcf()._plot_name = col

def anymeanvar(col):
    assert col in txtdata.keys
    track = readtrack.Track().track
    coldata = track[col]
    return np.mean(coldata), np.var(coldata)
    

################################################################################
# 

# 
def _sample_time(track, step=1):
    time = track['time']
    return time[np.arange(0, time.size, step)]

#
def cxy_vel(tr, step=1):
    time = tr.track['time']
    if time.size < step:
        return np.array([])
    cxy = tr.get_cxy()
    x, y = cxy[:,0], cxy[:,1]
    deltat = time[step] - time[0]
    idx = np.arange(0, time.size, step)
    disp = np.sqrt( (x[idx[1:]] - x[idx[:-1]])**2 + (y[idx[1:]] - y[idx[:-1]])**2 )
    return disp/deltat

def _vel(track, step=1):
    time = track['time']
    if time.size < step:
        return np.array([])
    x, y = track['x'], track['y']
    idx = np.arange(0, time.size, step)
    deltat = time[idx[1:]] - time[idx[:-1]]
    disp = np.sqrt( (x[idx[1:]] - x[idx[:-1]])**2 + (y[idx[1:]] - y[idx[:-1]])**2 )
    nonzero = np.nonzero(deltat)[0]
    zeroargs = np.argwhere(deltat== 0.)
    if not (np.all(disp[zeroargs]==0.)):
        problem = np.argwhere(np.logical_and(deltat==0., disp!=0.))
        print('Warning: In 0 time we are expected to have moved 0 distance. problem.size = {}'.format(problem.size))
    # avoid division by zero
    return disp[nonzero]/deltat[nonzero]

# step mean speed
def lvel(ltr, trail=False):
    dt = ltr.get_step_dt(trail)
    dx = ltr.get_step_dx(trail)
    vel = np.nan_to_num( dx / dt[:,np.newaxis] )
    speed = np.sqrt(np.sum(vel**2, axis=1))
    #weightd mean
    return np.sum(dt * speed)/np.sum(dt)

# instantaneous velocity  (actually speed)
# use this for simulation
def _inst_vel(track, step=1):
    time = track['time']
    x, y = track['x'], track['y']
    idx = np.arange(0, time.size, 1)
    deltat = time[idx[1:]] - time[idx[1:]-step]
    disp = np.sqrt( (x[idx[1:]] - x[idx[1:]-step])**2 + (y[idx[1:]] - y[idx[1:]-step])**2 )

    safety_check = False
    if safety_check:
        # avoid division by zero
        nonzero = np.nonzero(deltat)[0]
        zeroargs = np.argwhere(deltat== 0.)
        if not (np.all(disp[zeroargs]==0.)):
            problem = np.argwhere(np.logical_and(deltat==0., disp!=0.))
            print('Warning: In 0 time we are expected to have moved 0 distance. problem.size = {}'.format(problem.size))
        return disp[nonzero]/deltat[nonzero]
    else:
        return disp/deltat


def _disp(track, step=1):
    time = track['time']
    if time.size < step:
        return np.array([])
    x, y = track['x'], track['y']
    return __disp(time, x, y, step=step)
    
def _trail_disp(track, step=1):
    time = track['time']
    if time.size < step:
        return np.array([])
    x, y = track['trail_x'], track['trail_y']
    return __disp(time, x, y, step=step)

def __disp(time, x, y, step=1):
    idx = np.arange(0, time.size, step, dtype=int)
    deltat = time[idx[1:]] - time[idx[:-1]]
    disp = np.sqrt( (x[idx[1:]] - x[idx[:-1]])**2 + (y[idx[1:]] - y[idx[:-1]])**2 )
    return disp

@defaultsave(True)
def disp(step=1):
    # step is the number of indices to skip in the data, not a time
    # i.e. if data is collected at 0.1 second intervals
    #  step = 10 will give 1 second time resolution,
    #  step = 300 will give 30 second time resolution
    for tr in readtrack.trackset():
        vt = _disp(tr.track, step)
        tx = _sample_time(tr.track, step)
        plt.plot(tx[1:], vt)

"""
Carefully compute the decay time of velocity autocorrelation function 
"""
def vel_corr():
    trs = readtrack.trackset()
    return _vel_corr(trs)

@stats.keep
def _vel_corr(trs):
    # m_lifetime = stats.load_or_compute('lifetime', eventanalyse.lifetime)['mean']
    # pili_timescale = np.ceil(m_lifetime) # round up to whole number of seconds
    # # tmp
    # pili_timescale = min(pili_timescale, 4.)
    pili_timescale = 6.
    return __vel_corr(trs, pili_timescale)

def __vel_corr(trs, fit_maxtime):
    # velocity correlation function
    mpl.style.use('fivethirtyeight')
    vs = [_vel(tr) for tr in trs]
    """
    The autocorrelation is an exponetial decay followed by a long tail which is expected to be noise with constant mean. 
    fit_maxtime is used to cut the tail of the correlation
    """
    deltat = trs[0].tstep
    nsteps = int(np.round(fit_maxtime/deltat))
    # ...
    v_corr = [np.correlate(v,v,mode='full') for v in vs]
    assert(all((v_corr_i.size > 2*nsteps for v_corr_i in v_corr)))
    # slice and take the mean
    mpts = [int(v_corr_i.size/2) for v_corr_i in v_corr]
    each_corr = np.vstack([v_corr_i[mpt+1:mpt+nsteps+1] for mpt, v_corr_i in zip(mpts, v_corr)])

    time = deltat * np.arange(1,nsteps+1,1)
    # why take the mean of correlation functions rather than the mean of the fitted parameters?
    m_v_corr = np.mean(each_corr, axis=0)
    # a + b * exp(c*t)
    a, b, c = twutils.fit_decay(time, m_v_corr, plot=True)
    tau_decay = -c
    sd = {'velocity_decay_time':tau_decay}
    twutils.print_dict(sd)
    return sd

def set_resolution(trs, cut, init=1000):
    # cut is the time in seconds to cut the track with
    # deltat = 0.1 then init=1000 is 100 seconds
    for tr in trs:
        tr.filter_to_resolution(cut, initial=init)
    # cut out the initial relaxation period
    # TODO:
    # I tried to use (statistical) K-M test on the number of pili to determine where to cut
    # but could not reliably satisfy the test (\alpha=0.5)
    return trs


# plot velocity distribution
@command.defaultsave()
def vel():
    trs = readtrack.trackset()
    # v_tau = stats.load_or_compute('velocity_decay_time', lambda: _vel_corr(trs) )
    # set_resolution(trs, 2*v_tau)
    _vel_distribution(trs, res=200)
    ax = plt.gca()
    ax.set_xlim(0,3.0)

def trim_small(x, small):
    return x[x>small]

# instantaneous velocity distribution
def _vel_distribution(trs, label=None, hist=True , res=40):
    allvel = np.concatenate([_inst_vel(tr) for tr in trs])
    # allvel = twutils.trim_tail(allvel, 0.05)
    # allvel = trim_small(allvel, 0.01)
    lims = (None,None)
    ax = plt.gca()
    plotutils.ax_kdeplot(ax, allvel, res=res, linekw={'label':label},xlims=lims, hist=hist)
    ax.set_ylabel('P')
    ax.set_xlabel(r'velocity $\mu m s^{-1}$')
    ax.legend()
    plt.tight_layout()


# linearized track velocity distribution NOTES

# This is a mess because there are two procedures 
# (i) FJ linearization procedure used on FJ data to segment the track based on distance travelled 
# (ii) segmentation of the track using autocorrelation time to get approximately statistically independent segments
# For descriptive statistics I may choose to use one or neither
# For comparision with FJ data it makes sense to use (i) and if I want to compute std error then in addition I should use (ii)

# The way velocity calculation for (ii) is currently implemented is incorrect 
# I mixed up two different things
# (i) The timescale on which to compute velocity ( coarse graining timescale)
# (ii) The track segmentation timescale
# Because of our KMC procedure we should expert a reliable instantaneous velocity only by taking a timescale >> the relevent kmc timescale
# In this case \tau_{retraction} = 0.0053 s >> 0.1 seconds
# On the other hand obtaining a reliable instanteous velocity at 0.1 second time resolution for FJ data is generally not possible
# and that is the reason for the linearization procedure

# The general procedure for obtaining comparable velocity distribution then is to apply the linearization procedure and then take instantaneous velocity


@command.defaultsave()
def raw_vel():
    trs = readtrack.trackset()
    _vel_distribution(trs)

@command.defaultsave()
def raw_ltr_vel():
    trs = readtrack.trackset()
    ltrs = [_fj.linearize(tr) for tr in trs]
    _vel_distribution(ltrs)

@command.defaultsave()
def ltr_vel():
    # plotutils.default_style()
    trs = readtrack.trackset()
    ltrs = [_fj.linearize(tr) for tr in trs]
    ax = plt.gca()
    _ltr_vel(ax, ltrs)

def _ltr_vel(ax, ltrs, kdestyle={}):
    linvel = np.concatenate([np.linalg.norm(ltr.get_step_velocity(),axis=1) for ltr in ltrs])
    outd = plotutils.ax_kdeplot(ax, linvel, res=100, hist=True, **kdestyle)
    ax.set_ylabel('P')
    ax.set_xlabel(r'velocity $\mu m s^{-1}$')
    plt.tight_layout()




"""
Again. M-K test indicates a trend but there shouldn't be one. should run simulations for longer?
"""
import scipy.ndimage
import mkt
def pl_avg_mkt():
    trs = readtrack.trackset()
    tr = trs[0]
    pl_avg = tr['pl_avg']
    time = tr['time']
    plt.plot(time, pl_avg)
    smooth_pl_avg = scipy.ndimage.gaussian_filter1d(pl_avg.astype(float), 400)
    plt.plot(time, smooth_pl_avg)

    init = 1000
    Ha = 'up'
    MK, m, c, p = mkt.test(time[init:], pl_avg[init:], eps=1e-6, alpha=0.05, Ha=Ha)
    print('raw data ---')
    print('MK result: {}'.format(MK))
    print('slope: {}'.format(m))
    print('p-value: {}'.format(p))

    MK, m, c, p = mkt.test(time[init:], smooth_pl_avg[init:], eps=1e-6, alpha=0.05, Ha=Ha)
    print('smooth data ---')
    print('MK result: {}'.format(MK))
    print('slope: {}'.format(m))
    print('p-value: {}'.format(p))

    plt.show()


@defaultsave(True)
def all_vel(step=1):
    for tr in readtrack.trackset():
        vt = _vel(tr.track, step)
        tx = _sample_time(tr.track, step)
        plt.plot(tx[1:], vt)


def avgvel(step=1):
    m_vel = _avgvel(readtrack.trackset())
    print('avg velocity {} \mu/s'.format(m_vel))
    
def _avgvel(trs, step=1):
    """mean velocity for trackset"""
    return np.mean([np.mean(_vel(tr, step=step)) for tr in trs])

# average the number of steps used for md dynamics
def ansteps():
    stepcol = readtrack.Track().track['nsteps']
    nst = np.mean( stepcol[np.nonzero(stepcol)] )
    print('average nsteps ', nst)


#############################################################################
# comparable stats

import _analysefj, readmat
import fjstats
def getstats():
    trs = readtrack.trackset()
    stats = fjstats.getstats(trs)
    for name, stat in list(stats.items()):
        plt.clf()
        #plotutils.kdeplot(stat)
        plt.hist(stat, bins=5)
        dst_out = os.path.join(command.pdir, name+'_dst.png')
        plt.savefig(dst_out)


# reorientation
# ?
def dotax():
    trs = readtrack.trackset()
    ft = readmat.ExpTracks.usesimulated(trs)
    ft.compute_cxy()
    command.default_result_dirs()
    _analysefj.plot_frmodes(ft, fpercent=99)


###
# Compare distributions

# this directory
analysis_dir = os.path.dirname(os.path.abspath(__file__))

#https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test
def ks_measure(expvels, trs):
    """
    Comparing with the experimental displacement distribution
    """
    # recompute the simulated velocities
    disps = np.concatenate(_disp(trs))
    vels = 10*disps
    vels = fjstats.cut_for_rescaled(vels, 90)

    # load the fj distribution
    print('compute ks statistic')
    return scipy.stats.ks_2samp(expvels, vels)

def ks(trs):
    print('loading fjdata/disps.npy')
    expdisps = np.loadtxt(os.path.join(analysis_dir, 'fjdata/disps.npy'))
    measure = ks_measure(10*expdisps, trs)
    print(measure)
    return measure

def one_ks():
    trs = readtrack.trackset()
    
import simstats

def stime(step_d):
    print('computing step times with d  = ', step_d)
    st = readmat.this_fjtrack()
    simstats._steptime(st, step_d)

# --------------------------------------------------------------------------------

@stats.keep
def persistence():
    trs = readtrack.trackset()
    ltrs = [_fj.linearize(tr) for tr in trs]
    return mle_persistence(ltrs)

def mle_persistence(ltrs):
    # combine data
    u_t = []
    u_p = []
    for tr in ltrs:
        upart = tr.get_head_v()
        u_t.append(upart[1:])
        u_p.append(upart[:-1])
    u_t = np.concatenate(u_t)
    u_p = np.concatenate(u_p)

    qhat ,ahat = mle(u_t, u_p)
    sd = {}
    sd['qhat'] = {}
    sd['ahat'] = {} 
    sd['qhat']['estimate'] = qhat
    sd['ahat']['estimate'] = ahat 
    return sd

# --------------------------------------------------------------------------------
# periogram
import scipy.signal

@command.defaultsave()
def first_signal(period_threshold=20.0, n_annotate=4):
    trs = readtrack.trackset()
    ltrs = [_fj.linearize(tr) for tr in trs]
    plot_signal(ltrs[0])

def plot_signal(track, period_threshold=20.0, n_annotate=4):
    velocity = track.get_head_v()
    track_speed = np.linalg.norm(velocity, axis=1)

    window = 'blackmanharris'
    f, Pxx = scipy.signal.periodogram(track_speed, 10.0,
        window=window, scaling='density')

    peaks, properties = scipy.signal.find_peaks(Pxx)
    peaks_value = [Pxx[pid] for pid in peaks]
    descending = sorted(list(zip(peaks, peaks_value)), key= lambda t:t[1], reverse=True)

    fig = plt.figure(figsize=(10,5))
    ax = fig.add_axes([0,0,1,1]) 
    ax.plot(f, Pxx)
    ax.set_title('periodogram (blackmanharris)')
    ax.set_xlabel('frequency (Hz)')
    delta = (0.0,0.0)
    peak_count = 0
    for pv in descending:
        if peak_count >= n_annotate:
            break
        idx, value = pv
        period = 1./f[idx]
        if period > period_threshold:
            continue
        else:
            peak_count += 1
        xy = (f[idx], Pxx[idx])
        xyd = (xy[0]+delta[0], xy[1]+delta[1])
        x, y = xy
        plt.plot(x, y, marker='o', c='k', markersize=4.0)
        ax.annotate('{:.3f}s'.format(period), xy=xy, xycoords='data',
                    xytext=xyd 
                    # arrowprops=dict(facecolor='black', shrink=0.05)
                    # horizontalalignment='right', verticalalignment='top',
                    )
    return fig

# ----------------------------------------------------------------
# new polar plots and deviation angle calculation for linearised data

pdatafields = ['fastidx', 'dx_theta', 'body_theta', 'deviation', 'velocity']
# Polardata = collections.namedtuple('polardata', pdatafields)
class Polardata():
    def __init__(self, fastidx, dx_theta, body_theta, deviation, velocity):
        self.fastidx = fastidx
        self.dx_theta = dx_theta
        self.body_theta = body_theta
        self.deviation = deviation
        self.velocity = velocity
    
    def __iter__(self):
        return (getattr(self, attr) for attr in pdatafields)

def disp_angle(disp):
    n_disp = np.linalg.norm(disp, axis=1)[:,np.newaxis]
    da = disp[1:]/n_disp[1:]
    db = disp[:-1]/n_disp[:-1]
    theta = np.arccos(np.sum(da*db,axis=1))
    theta *= np.sign(np.cross(da, db))
    return theta

def angle(dx, dy):
    norm_dx = dx/np.linalg.norm(dx, axis=1)[:,np.newaxis]
    norm_dy = dy/np.linalg.norm(dy, axis=1)[:,np.newaxis]
    theta = np.arccos(np.sum(norm_dx * norm_dy,axis=1))
    theta *= np.sign(np.cross(norm_dx, norm_dy))
    return theta

def polar_hist(ax, arr_theta):
    nbins = 18
    barstyle = {'alpha':1.0, 'width':0.3,
        'linewidth':2.0, 'edgecolor':'k'}
    bins = np.linspace(-np.pi,np.pi,nbins+1,True)
    line, edges = np.histogram(arr_theta, bins=bins)
    centers = (edges[1:] + edges[:-1])/2
    handle = ax.bar(centers, line, **barstyle)
    return handle

def _allpolar(ltr, vthreshold):
    dx = ltr.get_step_dx()
    velocity = ltr.get_step_velocity()
    speed = np.linalg.norm(velocity, axis=1)
    fastidx = speed > vthreshold
    dx_theta = np.arctan2(dx[:,1],dx[:,0])
    lead = np.stack([ltr['x'],ltr['y']],axis=1)
    trail = np.stack([ltr['trail_x'],ltr['trail_y']],axis=1)
    b_ax = lead[ltr.step_idx] - trail[ltr.step_idx]
    body_theta = np.arctan2(b_ax[:,1], b_ax[:,0])

    # !
    theta = angle(dx, b_ax[1:])
    # theta = angle(dx, b_ax[:-1])
    
    # theta is nan if projected body length is zero
    theta = np.nan_to_num(theta)
    pdata = Polardata(fastidx, dx_theta, body_theta[1:], theta, speed)
    pdata.vthreshold = vthreshold
    return pdata

def allpolar(ltrs, vthreshold=0.3):
    _pdata = [_allpolar(ltr, vthreshold) for ltr in ltrs]
    data = []
    for attr in pdatafields:
        data.append( np.concatenate([getattr(p,attr) for p in _pdata]) )
    pdata = Polardata(*data)
    pdata.vthreshold = vthreshold
    return pdata

def polar_dev(axes, theta, fastidx):
    ax1, ax2 = axes
    ax1.set_title('deviation angle')
    sbar = polar_hist(ax1, theta[~fastidx])
    fbar = polar_hist(ax2, theta[fastidx])
    ax1.legend([sbar],['slow'])
    ax2.legend([fbar],['fast'])

def plotpolar(axes, polardata):
    fastidx, dx_theta, body_theta, theta, vel = polardata
    # --------------------
    ax1, ax2 = axes[0]
    ax1.set_title('head displacement')
    sbar = polar_hist(ax1, dx_theta[~fastidx])
    fbar = polar_hist(ax2, dx_theta[fastidx])
    ax1.legend([sbar],['slow'])
    ax2.legend([fbar],['fast'])

    # ---------------------
    ax1, ax2 = axes[1]
    ax1.set_title('body axis')
    sbar = polar_hist(ax1, body_theta[~fastidx])
    fbar = polar_hist(ax2, body_theta[fastidx])
    ax1.legend([sbar],['slow'])
    ax2.legend([fbar],['fast'])

    # ---------------------
    polar_dev(axes[2], theta, fastidx)

    vthreshold = polardata.vthreshold
    q = np.count_nonzero(vel > vthreshold)/vel.size
    axes[2][0].set_title(r'velocity $>$ {:3.1f} ({:d}\%)'.format(
        vthreshold, int(100*q)
    ))
    plt.gcf().tight_layout()

@command.defaultsave()
def polar(vthreshold=0.3,usepercentile=False):
    trs = readtrack.trackset()
    ltrs = [_fj.linearize(tr) for tr in trs]
    velocity = ltrs[0].get_step_velocity()
    speed = np.linalg.norm(velocity, axis=1)
    if usepercentile:
        vthreshold = np.quantile(speed, float(vthreshold)/100)
        print('using vthreshold ', vthreshold)
    
    polardata = allpolar(ltrs, vthreshold)
    fig, axes = plt.subplots(3, 2, subplot_kw=dict(polar=True), 
        figsize=(8,12))
    plotpolar(axes, polardata)

# ----------------------------------------------------------------

def unit(a):
    return a/np.linalg.norm(a)

# return dictionary with fields [dt, dx, velocity]
# here we take "actions" as steps.  
def _actions(ltrack):
    dt = ltrack.get_step_dt()
    dx = ltrack.get_step_dx()
    vel = dx / dt[:,np.newaxis]
    s = {
        'dt': dt,
        'dx': np.linalg.norm(dx,axis=1),
        'velocity': np.linalg.norm(vel,axis=1)
    }
    return s

def actions(ltrs):
    conc = {}
    data = [_actions(ltr) for ltr in tqdm(ltrs)]
    for k in data[0].keys():
        conc[k] = np.concatenate([d[k] for d in data])
    return conc
    
# fanjin uses log scale on histograms and its helpful to trim the data
npcount = np.count_nonzero
def trim_limits(arr, lim):
    left, right = lim
    left_out = arr<left
    right_out = arr>right
    N = arr.size
    form = 'trimming {} elements < {} && {} elements > {}'
    print(form.format(npcount(left_out)/N, left, npcount(right_out)/N, right))
    return arr[np.logical_and(~left_out, ~right_out)]

def plot_actiondata(actions):
    # A: step time 
    fig, axes = plt.subplots(1,3, figsize=(15,5))
    ax1, ax2, ax3 = axes
    vthreshold = 0.3
    fastidx = actions['velocity'] > vthreshold
    fast_dt = actions['dt'][fastidx]
    slow_dt = actions['dt'][~fastidx]
    fast_vel = actions['velocity'][fastidx]
    slow_vel = actions['velocity'][~fastidx]
    fast_dx = actions['dx'][fastidx]
    slow_dx = actions['dx'][~fastidx]

    hstyle = {'rwidth':0.9,'alpha':0.8,'histtype':'bar','stacked':True}
    dtlim = [0,20]
    color = ['#1f77b4', '#ff7f0e']
    def _dt_hist(ax, fast_dt, slow_dt):
        ax.set_xscale('log')
        edges = np.geomspace(0.09,dtlim[1],num=12)
        # print(edges)
        ax.hist([fast_dt,slow_dt], label=['fast','slow'], 
            bins=edges, **hstyle, color=color[::-1])
        ticks = [0.1,0.2,0.4,1.6,6.4,dtlim[1]]
        ax.set_xticks(ticks)
        ax.set_xticklabels(ticks)
        ax.set_xlabel('time (s)')
        ax.legend()
    _dt_hist(ax1, fast_dt, slow_dt)
    
    vlim = [0,100]
    def _dv_hist(ax, slow_vel, fast_vel):
        ax.set_xscale('log')
        edges = np.geomspace(0.0001,100,num=16)
        # ax.hist(vel, bins=edges, **hstyle)
        trim_slow_vel = trim_limits(slow_vel, [0,100])
        trim_fast_vel = trim_limits(fast_vel, [0,100])
        ax.hist([trim_fast_vel,trim_slow_vel], 
             label=['fast','slow'], bins=edges, color=color[::-1], **hstyle)
        # ax.hist(trim_fast_vel, bins=edges, **hstyle)
        ticks = [0.001,0.01,0.1,1.0,10.0]
        ax.set_xticks(ticks)
        ax.set_xticklabels(ticks)
        ax.set_xlabel(r'velocity $\mu m s^{-1}$')
        ax.set_xlim([0.0001,100])
        ax.axvline(0.3,linestyle='--', color='k')
        ax.legend()
    _dv_hist(ax2, slow_vel, fast_vel)

    def _d_hist(ax, slow_dx, fast_dx):
        ax.set_xscale('log')
        edges = np.geomspace(0.1,0.25,num=16)
        ax.hist([fast_dx,slow_dx], bins=edges, 
            label=['fast','slow'], color=color[::-1], **hstyle)
        ticks = [0.06,0.12,0.24]
        ax.set_xlim(0.06,0.24)
        ax.set_xticks(ticks,minor=False)
        ax.set_xticklabels(['%.2f' % t for t in ticks])
        ax.set_xticks([],minor=True)
        ax.set_xlabel(r'displacement ($\mu m$)')
        ax.legend()
    _d_hist(ax3, slow_dx, fast_dx)


# plot linearised velocity profiles with latex table
default = ["k_ext_off", "dwell_time", "pilivar",  "anchor_angle_smoothing_fraction", "k_spawn"]
xlabel = r"linearised speed $(\mu m\, s^{-1})$"
def plot_sim_with_table(ax, simdata_example,  par=default):
    # parameters
    args = parameters.thisread(directory=simdata_example)
    row = [str(args.pget(name)) for name in par]
    col = [[x] for x in row]
    df = pd.DataFrame({"parameter":par, "value":row})
    #
    simtrs = readtrack.trackset(ddir=join(simdata_example, "data/"))
    simltrs = [_fj.linearize(tr) for tr in simtrs]

    sim_vel = np.concatenate([ltr.get_step_speed() for ltr in simltrs])
    xlim = (0, np.quantile(sim_vel, 0.98))
    sns.histplot(sim_vel, binrange=xlim, ax=ax, stat="density")
    ax.grid(False)
    # ax.text(.5,.5, df.to_latex().replace('\n', ' '),
    #     transform=ax.transAxes, fontsize=20)
    ax.set_xlabel(xlabel)


###########################
# kmsd -- as best as i can compute it
# check the size of the trajectories first
# does it really make sense to compute MSD for  trajectories that displace by only a few cell widths

def fit_kmsd(basis, msd, weights=None):
    l_basis, l_msd = np.log(basis), np.log(msd)
    p, cov = np.polyfit(l_basis, l_msd, deg=1, cov=True)
    return p, cov

def kmsd_one(tr):
    low = 20
    ff = 2
    scaling = np.arange(low, int(tr['time'].size/ff)+1, 1).astype(int)
    if tr['time'].size < ff*low or scaling.size < 10:
        return None
    xyz = tr.get_head()
    x_, y_, _ = tr.get_origin()
    x, y = xyz[:,0] - x_, xyz[:,1] - y_ # subtract the starting position
    msd_n = np.full(scaling.size, np.nan)
    for i, window in enumerate(scaling):
        msd_n[i] = np.mean( 
                (x[window:] - x[:-window])**2 +  (y[window:] - y[:-window])**2 
                )
    # linear fit
    try:
        p, cov = fit_kmsd(scaling, msd_n)
    except ValueError:
        return None
    kmsd, _intercept = p
    if kmsd < 0.:
        return None
    return p, cov, scaling, msd_n


if __name__=='__main__':

    ff, thisargs = command.process(locals())
    ff(*thisargs)


