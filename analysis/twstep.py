#!/usr/bin/env python3

"""
Need some methods to operate on tracks and identify crossings and other on-step 
and off-step behaviour.
"""

import sys, os, math
import numpy as np
import matplotlib.pyplot as plt

import datapaths
import command
import parameters
from command import defaultsave
import readtrack

import astrack


"""
We need data on the step surface to do this correctly so run analysis from the simulation
directory and read the local config.txt 
"""

# This parameter is the distance inside the cell surface that pili anchors are located
# Use it as the maximum distance that may intersect the surface.
# Naturally I << args.cell.R
inside = I = 0.04

try:
    args = parameters.thisread()
    surface = args.surface
    sep = surface.sep
    height = surface.height
    R = args.cell.R
    length = args.cell.length
except IOError:
    sep = 10.
    height = 5.
    R = .5
    length = 2.

# for testing purpose only
# extrfile = '/data/dan/run/dec_18/stepped_surface/run_low_step_again/height_001.000_eps_00025.0000_pilivar_00.78540/data/bacterium_00000.dat'
# tr = readtrack.Track(extrfile)

import shapeplot
@command.defaultsave()
def steplt():
    ax = plt.gca()
    shapeplot.longtracks(ax, readtrack.trackset())
    xlim = ax.get_xlim()
    draw_step_lines(ax, xlim)

##################################################################################
# utilities


def transform(x):
    #return floating point value in 'step-space' where numbers close to integers are close to steps
    return x/sep + 0.5

def int_close(x, tol=R/sep):
    return np.abs(np.round(x) - x) < tol
    
def locate(x):
    """
    Given an x coordinate which x-y plane surface are we close to.
    x-y surfaces are numbered as integers starting at 0.
    """
    return int(math.floor(x/sep + 0.5))

def nplocate(x):
    # locate function for numpy arrays
    return np.array(np.floor(x/sep + 0.5), dtype=int)

def xextents(tr):
    xh = tr.track['x']
    xt = tr.track['trail_x']
    # separate into xmin and xmax
    flipidx = xh < xt
    xmax = np.where(flipidx, xt, xh)
    xmin = np.where(flipidx, xh, xt)
    return xmin-I, xmax+I


guidestyle = {'color':'k', 'linestyle':'--', 'alpha':0.4}
spanstyle = {'alpha':0.06, 'color':'k'}

def draw_step_lines(ax, xlims):
    args = parameters.thisread()
    sep = args.surface.sep
    """surface separation"""
    xn, xm = xlims
    xdist = xm - xn
    nsteps = math.ceil(xdist/sep) + 1
    # draw this many steps in both directions
    xgrid = np.arange(-nsteps*sep -5., nsteps*sep+5.+sep/2., sep)
    draw_const_lines(ax, xgrid, guidestyle)
    draw_span(ax, xgrid, spanstyle)

def draw_const_lines(ax, xgrid, style):
    lines = []
    for i, x in enumerate(xgrid):
        line = ax.axvline(x, **style)
        lines.append(line)
    return lines

def draw_span(ax, xgrid, style):
    for i, x in enumerate(xgrid):
        # draw shaded regions for the low step
        if locate(x) % 2 == 0 and i+1 != xgrid.size:
            ax.axvspan(x, xgrid[i+1], **style)



#################################################################################

"""
Let functions have signature f(tr) where tr is a Track object 
and data is accessed like tr.track['x']

Compile a list of useful data
1. currently bridging step, locate(xmin) != locate(xmax)
2. head/tail currently touching step, int_close(transform(xmin||xmax))
3. Crossings and non-crossings 
4. visit map integrated over y-direction 
5. Pili surface attachments which land across steps 
6. plot the rotation of the body during crossing events
"""
#... identity crossings and non-crossings ...

# helper methods

def _headtail_touching(tr):
    xh = tr.track['x']
    xt = tr.track['trail_x']
    ith = int_close(transform(xh))
    itt = int_close(transform(xt))
    return ith, itt

def _is_touching(tr):
    ith, itt = _headtail_touching(tr)
    return np.asarray(np.logical_or(ith, itt), dtype=int)

# data holder class

class CrossingData(object):
    """
    fields are:
    to_touching
    from_touching
    crossingidx
    cross_upstep
    cross_downstep
    """
    def __init__(self, to_touching, from_touching, tr):
        """
        to_touching/from_touching are the indices where we start/stop touching
        x is just a reference to x position of head
        """
        # arrays with size = no. 'crossing attempts'
        self.to_touching = to_touching
        self.from_touching = from_touching
        self.tr = tr # reference to Track object
        x = self.tr.track['x']
        pre_locate = nplocate(x[to_touching]) 
        post_locate = nplocate(x[from_touching])

        # arrays with count_nonzero() = no. crossings
        crossingidx = pre_locate != post_locate
        self.crossingidx = crossingidx
        cross_upstep = np.logical_and(pre_locate % 2 == 0, crossingidx)
        cross_downstep = np.logical_and(pre_locate % 2 == 1, crossingidx)
        self.cross_upstep = cross_upstep
        self.cross_downstep = cross_downstep

    def listxcross(self, x):
        """ For a 'crossing attempt' 
        record every time the value x actually passes the line
        """
        listc = []
        for tidx, fidx in zip(self.to_touching, self.from_touching):
            window = x[tidx:fidx]
            ll = nplocate(window)
            # neighbour diff
            crossed = np.argwhere( (ll[1:]-ll[:-1]) != 0 ).flatten()
            ct = crossed + tidx
            listc.append(ct)
        return listc

    def get_n_crossings(self):
        return np.count_nonzero(self.crossingidx)

    def get_n_upsteps(self):
        return np.count_nonzero(self.cross_upstep)

    def get_n_downsteps(self):
        return np.count_nonzero(self.cross_downstep)

    def get_upsteps(self):
        touchidx = self.cross_upstep
        start, end = self.to_touching[touchidx], self.from_touching[touchidx]
        return start, end

    def get_downsteps(self):
        touchidx = self.cross_downstep
        start, end = self.to_touching[touchidx], self.from_touching[touchidx]
        return start, end

    @classmethod
    def touchstate(cls, tr):
        """ Construct touching data object from track
        """
        # find the indices of the track immediately before the touching state changes
        tarr = _is_touching(tr)
        tdiff = tarr[1:] - tarr[:-1]
        # no idea why np.argwhere makes column vectors
        to_touching = np.argwhere(tdiff == 1).flatten()
        from_touching = np.argwhere(tdiff == -1).flatten()
        # the size of these arrays should not be different by more than 1
        assert(abs(to_touching.size - from_touching.size) <= 1)
        #set the array sizes equal
        to_touching = to_touching[:from_touching.size]

        return cls(to_touching, from_touching, tr)

# methods for rtw.py
def _utouchstate(trs):
    return list(map(CrossingData.touchstate, trs))

def _crosscount(crs):
    return [cr.get_n_crossings() for cr in crs]

def _upcrosscount(crs):
    return [cr.get_n_upsteps() for cr in crs]

def _downcrosscount(crs):
    return [cr.get_n_downsteps() for cr in crs]



################################################################################
# Apply functions to set of tracks

def where():
    trs = readtrack.trackset()
    for i, tr in enumerate(trs):
        cdata = CrossingData.touchstate(tr)
        print('{:05d} no. crossings {}'.format(i, cdata.get_n_crossings()))

def crossing():
    trs = readtrack.trackset()

    tr = trs[1]
    cdata = CrossingData.touchstate(tr)
    ax = plt.gca()
    _crossing_axz(ax, cdata, tr)

    #... fill in with plotting methods for crossing data ...

def crossing_ncontacts():
    ax = plt.gca()
    trs = readtrack.trackset()
    cdatas = [CrossingData.touchstate(tr) for tr in trs]
    _crossing_ncontacts(ax, cdatas, trs)

sfdir = 'plots/crossing_contacts/'
sform_up = os.path.join(sfdir, 'cell_{:05d}_upcrossing_{:04d}.png')
sform_down = os.path.join(sfdir, 'cell_{:05d}_downcrossing_{:04d}.png')

def _crossing_ncontacts(ax, cdatas, trs):
    datapaths.force_mkdir(sfdir)
    for i, tr in enumerate(trs):
        cdata = cdatas[i]
        start, end = cdata.get_upsteps()
        for j, rsd in enumerate(zip(start, end)):
            st, ed = rsd
            fname = sform_up.format(i, j)
            ax.clear()
            print('writing to ', fname)
            _one_crossing_ncontacts(ax, tr, st, ed, fname)

        start, end = cdata.get_downsteps()
        for j, rsd in enumerate(zip(start, end)):
            st, ed = rsd
            fname = sform_down.format(i, j)
            ax.clear()
            print('writing to ', fname)
            _one_crossing_ncontacts(ax, tr, st, ed, fname)


def _one_crossing_ncontacts(ax, tr, start, end, fname, line_kw={}):
    ax.plot(np.arange(end-start)/10, tr.track['ncontacts'][start:end], **line_kw)
    ax.set_ylabel(r'No. Contacts')
    ax.set_xlabel(r'time (s)')
    plt.savefig(fname)

def _crossing_moments(cdata):
    listc = cdata.listxcross(cdata.tr.track['x'])
    #print np.argwhere(cdata.cross_upstep).fl
    upmoments = [listc[i] for i in np.argwhere(cdata.cross_upstep).flatten()]
    downmoments = [listc[i] for i in np.argwhere(cdata.cross_downstep).flatten()]
    return upmoments, downmoments

def all_crossing_axz():
    trs = readtrack.trackset()
    cdatas = [CrossingData.touchstate(tr) for tr in trs]

    formdir = 'plots/crossing_axz/'
    datapaths.force_mkdir(formdir)
    ax = plt.gca()
    for i, (tr, cdata) in enumerate(zip(trs, cdatas)):
        ax.clear()
        count = _crossing_axz(ax, cdata, tr)
        if count:
            form = os.path.join(formdir, 'crtilt_{:02d}.png'.format(i))
            plt.savefig(form)

def _crossing_axz(ax, cdata, tr):
    ax.set_ylim(-np.pi/2., np.pi/2.)
    theta = astrack.tiltangle(tr.track)
    upmoments, downmoments = _crossing_moments(cdata)
    count = 0
    for i, (start, end) in enumerate(zip(*cdata.get_upsteps())):
        count += 1
        ax.plot(np.arange(end-start)*0.1, theta[start:end], c='b')
        ax.plot((upmoments[i]-start)*0.1, theta[upmoments[i]], 
                linestyle=None, marker='|', c='k', markersize=40)
    for i, (start, end) in enumerate(zip(*cdata.get_downsteps())):
        count += 1
        ax.plot(np.arange(end-start)*0.1, theta[start:end], c='g')
        print('dm', downmoments[i])
        ax.plot((downmoments[i]-start)*0.1, theta[downmoments[i]], 
                linestyle=None, marker='|', c='k', markersize=40)
    ax.set_ylabel(r'tilt angle')
    ax.set_xlabel(r'time (s)')
    return count

#@defaultsave(True)
def touching():
    trs = [readtrack.trackset()[0]]
    touching_data = [_is_touching(tr) for tr in trs]
    for trdata in touching_data:
        plt.plot(trdata)
    plt.show()


################################################################################
# Visit Map

def integrated_visit(ax, trs):
    import scipy.stats
    xh = np.concatenate([tr.track['x'] for tr in trs])
    xlims = min(xh), max(xh)
    ax.set_xlim(*xlims)

    lines = []
    for tr in trs:
        xh = tr.track['x']
        kernel = scipy.stats.gaussian_kde(xh)
        xspace = np.linspace(xlims[0], xlims[1], 1000)
        vmap = list(map(kernel.evaluate, xspace))
        line = ax.plot(xspace, vmap)
        lines.append(line)

    draw_step_lines(ax, xlims)
    return lines


@defaultsave()
def xvisit():
    trs = readtrack.trackset()
    #fig = plt.figure(figsize=(10.5, 5.))
    ax = plt.gca()
    ax.set_xlabel('x', fontsize=20.)
    ax.set_ylabel('Visit density', fontsize=20.)
    ax.axes.yaxis.set_ticks([])
    integrated_visit(ax, trs)

if __name__=='__main__':
    ff, thisargs = command.process(locals())
    ff(*thisargs)

    # testing
    #xextents(tr)




