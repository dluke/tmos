#!/usr/bin/env python3

import sys, os
import numpy as np
import matplotlib.pyplot as plt
import command
from tqdm import  tqdm
from command import defaultsave

import filesystem
import parameters
import plotutils 
import shapeplot

import stats
import twutils
import readtrack
from readtrack import Track, trackset

# Reserve this file for short tracking functions

"""
Write functions to operate on lists of Tracks and then add auto-generate 
convenience methods which load the original Tracks on fly.

Set <method>.trs = True  for methods which can operate on lists of Tracks
Set <method>.pili = True for methods which can operator on binding/unbinding data

Consider that in general we want to split loading, analysis and plotting
"""


from shapeplot import apply_surface_annotation

        
###
# plotutils.default_style()

@defaultsave(svg=True)
def lt(n=None):
    ax = plt.gca()
    shapeplot.longtracks(ax, trackset()[:n])

@defaultsave(svg=True)
def surfacelt():
    ax = plt.gca()
    args = parameters.thisread()
    apply_surface_annotation(args)(shapeplot.longtracks)(ax, trackset())

@defaultsave(svg=True)
def st():
    ax = plt.gca()
    shapeplot.shapetracks(ax, trackset())

@defaultsave()
def shapelt():
    ax = plt.gca()
    shapeplot.ltdraw(ax, trackset() )

# testing
def test_lt():
    import _fj
    trs = _fj.npyload(2538).get_tracklike()
    shapeplot.longtracks(plt.gca(), trs)
    plt.show()

# radius of gyration

def radgy_dst(trs):
    rgy = np.zeros(len(trs))
    for i, tr in enumerate(trs):
        cxy = tr.get_cxy()
        M = np.mean(cxy, axis=0)
        r = np.subtract(cxy, M) # subtract centre of mass
        rgy[i] =  np.sqrt(np.mean(r**2))
    return rgy


def radgy(trs):
    return np.mean(radgy_dst(trs))

# kmsd for simulated trajectory
# we assume that each trajectory is a sample of some hypothertical population

# TODO: is it better to compute kmsd for each trajectory and take mean or to take means of msd(\tau) first?

@stats.keep
def kmsd():
    def kmsd_one(tr):
        use_overlap = False
        maxwindow = tr.size//10
        return msd_no_overlap(tr, maxwindow, use_overlap=use_overlap)
    trs = readtrack.trackset()
    kmsd_arr = [kmsd_one(tr) for tr in trs]
    sd = {}
    sd['kmsd'] = stats.stats(kmsd_arr)
    twutils.print_dict(sd)
    return sd
    
# any reason for doing this?
def msd_no_overlap(tr, maxwindow=200, use_overlap=False, c_mindata=1000):
    # FJ data max length is 2000 seconds, 2e5 data points, but many tracks are shorter

    # TODO: check the goodness of these linear fits by eye
    # how to set minwindow?
    minwindow = 10
    scaling = np.arange(minwindow, maxwindow, 1).astype(int)
    def msd_one(tr):
        msd_n = np.full(scaling.size, np.nan)
        weights = np.full(scaling.size, np.nan)
        if tr.size < c_mindata:
            return np.nan
        time = tr['time']
        xyz = tr.get_head()
        x_, y_, _ = tr.get_origin()
        x, y = xyz[:,0] - x_, xyz[:,1] - y_ # subtract the starting position
        
        for i, window in enumerate(scaling):
            if not use_overlap:
                end_to_end = np.arange(0, tr.size, window)
                start_window = end_to_end[:-1]
                end_window = end_to_end[1:]
                msd_n[i] = np.mean(
                    (x[start_window]-x[end_window])**2 + (y[start_window]-y[end_window])**2
                    )
                weights[i] = end_to_end.size - 1
            else:
                msd_n[i] = np.mean(
                    (x[window:] - x[:-window]) ** 2 + (y[window:] - y[:-window]) ** 2
                )
                weights[i] = x.size - window + 1
        # TODO we use sqrt(N) for weightings, is this the right 
        kmsd, _ = fit_kmsd(scaling, msd_n, np.sqrt(weights))
        return kmsd 
    return msd_one(tr)

def fit_kmsd(basis, msd, weights=None):
    l_basis, l_msd = np.log(basis), np.log(msd)
    p, cov = np.polyfit(l_basis, l_msd, deg=1, cov=True)
    return p, cov

def kmsd_one(tr):
    low = 20
    ff = 2
    scaling = np.arange(low, int(tr['time'].size/ff)+1, 1).astype(int)
    if tr['time'].size < ff*low:
        return None
    assert(scaling.size > 1)
    xyz = tr.get_head()
    x_, y_, _ = tr.get_origin()
    x, y = xyz[:,0] - x_, xyz[:,1] - y_ # subtract the starting position
    msd_n = np.full(scaling.size, np.nan)
    for i, window in enumerate(scaling):
        msd_n[i] = np.mean( 
                (x[window:] - x[:-window])**2 +  (y[window:] - y[:-window])**2 
                )
    # linear fit
    p, cov = fit_kmsd(scaling, msd_n)
    kmsd, _intercept = p
    if kmsd < 0.:
        return None
    return p, cov, scaling, msd_n


def kmsd_dst(trs):
    """
    return an msd value for each track by using a sliding window
    """
     # averaging window size range
    dst = np.array([kmsd_one(tr) for tr in  tqdm(trs)])
    print()
    # some values are invalid because the data is insufficient 
    validdst = dst[~np.isnan(dst)]
    print('cut {}/{} tracks from kmsd distribution'.format(dst.size-validdst.size, dst.size))
    return validdst

# for step, x and y msd

def _msdx(trs):
    def msdx_one(tr):
        x =  tr.track['x']
        return (x[-1] - x[0])**2
    return np.mean(list(map(msdx_one, trs)))

def _msdy(trs):
    def msdx_one(tr):
        x =  tr.track['y']
        return (x[-1] - x[0])**2
    return np.mean(list(map(msdx_one, trs)))


# simple pili average length
def _plengths(trs):
    return np.mean([np.mean(tr.track['pl_avg']) for  tr in trs])

## tilt angle

# one track
def tiltangle(track):
    axx, axy = track['ax_x'], track['ax_y']
    tilt = np.arccos(np.clip(np.sqrt(axx*axx + axy*axy), 0, 1))
    sgn = np.where(track['ax_z'] >= 0., -1, 1.)
    return tilt * sgn


# all tracks
def ptilt_f(ax, trs, label=False):
    """
    tilt angles vs time for every track
    """
    for tr in trs:
        tilt = tiltangle(tr.track)
        ax.plot(tr.track['time'], tilt)
    return ax
# default save 
@defaultsave(True)
def ptilt():
    ax = ptilt_f(plt.gca(), trackset())
    ax.set_ylabel('tilt angle')
# average 
def avgtilt(trs):
    # wtf?
    def first(a): return [a[0]] 
    return np.mean([tiltangle(tr.track).mean() for tr in first(trs)])

def avgncontact(trs):
    return np.mean([tr['ncontacts'].mean() for tr in trs])

def _pncontact(ax, trs):
    for tr in trs:
        ax.plot(tr['ncontacts'])


@defaultsave(True)
def pncontact():
    _pncontact(plt.gca(), trackset())


# plot track using scatter
@defaultsave(True)
def plottrack(idx=None, tmeslice=None):
    # get x, y coordinates
    #tf = get_track_file(idx)
    tr= Track()
    track = tr.slice(tmeslice)

    x = track['x']
    y = track['y']
    events = track['process']

    # search out the release events
    releases = np.argwhere(events == 'release').flatten()

    xrel = x[releases]
    yrel = y[releases]

    plt.plot(x, y, '-')
    plt.plot(xrel, yrel, 'ro') 
#plottrack.defaults=[0, (100.,104.)]

def rate_count(tmeslice=None):
    """
    print the number of times a process happened
    mainly for debugging
    """
    tr = Track()
    track = tr.slice(tmeslice)

    # sum up the frequencies of each process
    process = track['process']
    freq = {}
    for proc in process:
        if proc not in freq:
            freq[proc] = 0
        else:
            freq[proc] += 1
    print(freq)


def list_ftime():
    """How much simulation did we run before crashing or raising exception"""
    tfs = readtrack.trackset()
    ftime = []
    for i, tr in enumerate(tfs):
        ft = tr.track['time'][-1]
        ftime.append(ft)
        #print '{} {}'.format(i, ft)

    args = parameters.thisread()
    mxt = args.system.simtime
    # binary sucess measure
    def issuccess(mxt):
        return 1 if np.all(np.array(ftime) > mxt) else 0
    # analogy success measure
    def success_value(mxt):
        return np.mean(ftime)/mxt
    return success_value(mxt)
# convenience method for printing results ot stdout
def ft(): print(list_ftime())


if __name__=='__main__':
    ff, thisargs = command.process(locals())
    ff(*thisargs)
