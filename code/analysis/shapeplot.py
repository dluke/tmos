

"""
Use matplotlib.patches to draw the shape of the cell from trs.track
"""
import math
import itertools

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

import matdef
import parameters
import plotutils

verbose = False

##args = parameters.read()
#R = args.cell.R
#D = 2*R
#l = args.cell.length

R = 0.5
D = 2*R 
l = 2


### the same utilities

def tonpy(v):
    return np.array([v.x, v.y, v.z])


#################################################################################
# surface annotations

def draw_sine_crest(ax, args, xlims):
    guidestyle = {'color':'k', 'linestyle':'--', 'alpha':0.4}
    assert(args.surface.shape == 'sineplane')
    sep = 2 * np.pi * args.surface.B
    print('surface.B = ', args.surface.B)
    print('shapeplot annotation period', sep)
    offset = sep/4.
    xn, xm = xlims
    #nsteps = math.ceil((xm-xn)/sep) + 1
    nleft = math.ceil(abs(xn)/sep) 
    nright = math.ceil(abs(xm)/sep) 
    xgrid = np.arange(-nleft*sep + offset, nright*sep + offset + offset/2., sep)
    def draw_const_lines(ax, xgrid, style):
        for i, x in enumerate(xgrid):
            line = ax.axvline(x, **style)
    draw_const_lines(ax, xgrid, guidestyle)
    return ax

def draw_hex_grid(ax, args, xlims, ylims, usepatches=False):
    guidestyle = {'color':'k', 's':1, 'alpha':0.5}
    assert(args.surface.shape == 'hexspheregrid')
    N = 10
    from tmos import surface
    # create object
    hxsph = surface.HexSphereGrid(args.surface.sphere_radius)
    # list of Hexc objects
    hxrange = hxsph.coordinate_range(surface.Hexc(0,0), N)
        
    def convert(hx):
        # convert Hexc position to numpy point
        trialpoint = tonpy( hxsph.get_xyz(hx) )[:2]
        # filter out vertices which are out of scope
        x, y = trialpoint
        if x < xlims[0] or x > xlims[1] or y < ylims[0] or y > ylims[1]:
            return None
        return trialpoint
    def apply_convert(hxrange):
        ls = []
        for hx in hxrange:
            pt = convert(hx)
            if pt is not None:
                ls.append(pt)
        return ls
    np_pts = np.array(apply_convert(hxrange))
    x, y = np_pts[:,0], np_pts[:,1]

    # draw scatter plot
    ax.scatter(x, y, **guidestyle)
    
    # also draw circular patches with no fill
    if usepatches:
        patchstyle = {'color':'k', 'fill':False, 'linewidth':1, 'alpha':0.5}
        for pt in np_pts:
            sph = patches.Circle(pt, args.surface.sphere_radius, **patchstyle)
            ax.add_patch(sph)
    return ax

from functools import wraps
def apply_surface_annotation(params=None):
    def annotate(f):
        """ args needs to be an optional argument of f or else passed to decorator """
        @wraps(f)
        def f_mod(*a, **kw):
            # two chances to find args != None 
            args = kw.pop('args', None)
            if args is None:
                args = params
            
            #
            ax = f(*a, **kw)
            if args and args.surface.shape == 'sineplane':
                draw_sine_crest(ax, args, ax.get_xlim())
            if args and args.surface.shape == 'hexspheregrid':
                draw_hex_grid(ax, args, ax.get_xlim(), ax.get_ylim())
            # add more surface annotations 
            # ...
            return ax
        return f_mod
    return annotate

#################################################################################

def shapetracks(ax, trs, linekw={}):
    longtracks(ax, trs, linekw)
    args = parameters.thisread()
    if args.surface.shape == 'sineplane':
        draw_sine_crest(ax, args, ax.get_xlim())
    return ax

# plot track using scatter
def longtracks(ax, trs, linekw={}):
    for tr in trs:
        track = tr.track
        x = track['x']
        y = track['y']
        ax.plot(x, y, '-', **linekw)
    plotutils.prettyax_track(ax, trs)
    ax.set_aspect('equal')
    return ax
longtracks.trs = True


colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
def ltdraw(ax, trs, sample=10):
    citer = itertools.cycle(colors)
    for n, tr in enumerate(trs):
        if verbose:
            print('drawing shapes for track', n)
        sim_capsdraw(ax, tr, sample, color=next(citer))
    plotutils.prettyax_track(ax, trs)
    longtracks(ax, trs)
    return ax

def _project_xy(ax_x, ax_y):
    # tails to project vertical axis onto 0,0
    nm = np.sqrt(ax_x**2+ ax_y**2)
    return ax_x/nm, ax_y/nm

def sim_capsdraw(ax, tr, sample, color='k'):
    ax.set_aspect('equal')
    track = tr.track
    step = float(sample)/float(tr.track['time'][-1]) * tr.size
    idx = np.array(np.arange(0, tr.size, step), dtype=int)
    if verbose:
        print("attempt to draw {} cell outlines".format(idx.size))
    x, y = tr.track['x'][idx], tr.track['y'][idx]
    tx, ty = tr.track['trail_x'][idx], tr.track['trail_y'][idx]
    # ax_x, ax_y = _project_xy(tr.track['ax_x'][idx], tr.track['ax_y'][idx])
    # why not compute axes from head,tail positions ...
    ax_x, ax_y = x-tx, y-ty
    norm = np.sqrt(ax_x**2 + ax_y**2)
    ax_x, ax_y = ax_x/norm, ax_y/norm
    # todo
    return capsdraw(ax, (x,y), (tx,ty), (ax_x,ax_y), 
            color=color, hashead=True)

def exp_capsdraw(ax, yr, ts=30., tnum=None):
    """Draw cell body shapesof experimental data using axis and centroid and dimensions"""
    ctracks = yr.get_columns(matdef.CENTROID, timestep=ts, tracknum=tnum)
    lw = yr.get_columns(['length', 'width'], timestep=ts, tracknum=tnum)
    yr.compute_ax()
    allaxes = yr.get_columns(['ax_x', 'ax_y'], timestep=ts, tracknum=tnum)
    citer = itertools.cycle(colors)
    for i, ctrack in enumerate(ctracks):
        # unpack
        cx, cy = ctrack
        length, width  = lw[i]
        ax_x, ax_y = allaxes[i]
        #
        lhalf = (length-width)/2.
        x = cx + ax_x * lhalf
        y = cy + ax_y * lhalf
        tx = cx - ax_x * lhalf
        ty = cy - ax_y * lhalf
        capsdraw(ax, (x,y), (tx,ty), (ax_x,ax_y), 
                R=width/2., color=next(citer), hashead=False)
    return ax


# by returning a list of Matplotlib Artists we can animate the drawing 
def capsdraw(ax, xy, txy, ax_xy, R=R, color='b', hashead=True, style={}):
    D = 2*R
    headalpha = 0.2
    # unpack
    x, y = xy
    tx, ty = txy
    ax_x, ax_y = ax_xy
    ax_x, ax_y = R * ax_x, R * ax_y

    #
    theta = np.degrees(np.arctan2(ax_y, ax_x))
    side_x = x - ax_y
    side_y = y + ax_x
    uside_x = x + ax_y
    uside_y = y - ax_x
    tside_x = tx - ax_y
    tside_y = ty + ax_x
    tuside_x = tx + ax_y
    tuside_y = ty - ax_x

    #artists = []
    for i in range(x.size):
        l1 = ax.plot([side_x[i], tside_x[i]], [side_y[i], tside_y[i]], color=color, **style)
        l2 = ax.plot([uside_x[i], tuside_x[i]], [uside_y[i], tuside_y[i]], color=color, **style)
        p1 = patches.Arc((x[i],y[i]),D,D, theta[i],
                #theta[i]-90., theta[i]+90.,
                color=color, **style)
        p2 = patches.Arc((tx[i],ty[i]),D,D, theta[i]+180,
                #theta[i]+90., theta[i]-90.,
                color=color, **style)
        arc1 = ax.add_patch(p1)
        arc2 = ax.add_patch(p2)
        artists = [l1, l2, arc1, arc2]
        if hashead:
            headp = patches.Circle((x[i],y[i]),R,fill=True, alpha=headalpha, color=color)
            head = ax.add_patch(headp)
            artists.append(head)
    #return artists
    


if __name__=='__main__':
    import readtrack
    trs = readtrack.trackset()
    ax = plt.gca()
    ax = capsdraw(ax, trs[0], 300)
    plt.show()



