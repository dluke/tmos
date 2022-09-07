import os, sys
import numpy as np
import math, itertools
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.patches as patches

from tqdm import tqdm
import parameters

import shapeplot
import twanalyse
import astrack
import twstep

""" Animations """
import readtrack 


"""
In addition to matplotlib animations
It might be time to either
1. Improve the way we store vtk surface data 
2. Reconstruct vtk files from the text output and the config.txt
"""


# 

R = 0.5
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
norm = np.linalg.norm

class CellArtist(object):
    """
    The cell shapes we draw are composed of two lines, two circles/arcs and a shaded circle.
    Want to hide the implementation to update all this matplotlib.patches objects.

    Note: Not derived from matplotlib artist. This could be called 'Cell Artist Manager'.
    """
    def __init__(self, l1, l2, p1, p2, phead):
        self.l1 = l1
        self.l2 = l2
        self.p1 = p1
        self.p2 = p2
        self.phead = phead

    @property
    def artists(self):
        return [self.l1, self.l2, self.p1, self.p2, self.phead]

    def set_data(self, triter_data):
        xy, txy, ax_xy = triter_data
        x, y = xy
        tx, ty = txy
        norm = np.sqrt( (x-tx)**2 + (y-ty)**2 )
        ax_x, ax_y = (x-tx)/norm, (y-ty)/norm
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

        #
        self.l1.set_data([side_x, tside_x], [side_y, tside_y])
        self.l2.set_data([uside_x, tuside_x], [uside_y, tuside_y])
        self.p1.center = (x,y)
        self.p2.center = (tx,ty)
        self.phead.center = (x,y)
        # the positions are updated by changing the objects directly but here we return the artist list anyway
        return self.artists 
    
    def set(self, **kw):
        for artist in [self.l1,self.l2,self.p1,self.p2,self.phead]:
            artist.set(**kw)

    @classmethod
    def init(cls, ax, color='b'):
        """
        Consruct the artist list for the first time.
        """
        style = {"linewidth": 2}
        l1, = ax.plot([], [], color=color, **style)
        l2, = ax.plot([], [], color=color, **style)

        p1 = patches.Arc((0,0), 2*R, 2*R, 0, color=color, **style)
        p2 = patches.Arc((0,0), 2*R, 2*R, 0, color=color, **style)
        ax.add_patch(p1)
        ax.add_patch(p2)
        
        headp = patches.Circle((0,0),R,fill=True, alpha=0.3, color=color)
        ax.add_patch(headp)
        new = cls(l1,l2,p1,p2,headp)
        new.set(zorder=2)
        return new

# for experimental data
class ExpCellArtist(CellArtist):
    """
    override CellArtist so that we take data like
    data = xy, ax_xy, length, width
    where xy is now the tip of the head so with should subtract width/2. along axis
    """

    def _set_radius(self, r):
        self.p1.width = r
        self.p1.height = r
        self.p2.width = r
        self.p2.height = r
        self.phead.width = r
        self.phead.height = r

    def set_data(self, triter_data):
        xy, ax_xy, length, width = triter_data
        x, y = xy
        ax_x, ax_y = ax_xy
        rlength = length-width
        hax_x, hax_y = width/2. * ax_x, width/2. * ax_y
        ax_x, ax_y = rlength * ax_x, rlength * ax_y
        # 
        x -= hax_x,
        y -= hax_y
        #
        side_x = x - hax_y
        side_y = y + hax_x
        uside_x = x + hax_y
        uside_y = y - hax_x
        tside_x = side_x - ax_x
        tside_y = side_y - ax_y
        tuside_x = uside_x - ax_x
        tuside_y = uside_y - ax_y
        tx = x - ax_x
        ty = y - ax_y
        #
        self.l1.set_data([side_x, tside_x], [side_y, tside_y])
        self.l2.set_data([uside_x, tuside_x], [uside_y, tuside_y])
        self.p1.center = (x,y)
        self.p2.center = (tx,ty)
        self.phead.center = (x,y)
        self._set_radius(width)
        # the positions are updated by changing the objects directly but here we return the artist list anyway
        return self.artists 

def shape_iter(tr, idx):
    x, y = tr.track['x'][idx], tr.track['y'][idx]
    tx, ty = tr.track['trail_x'][idx], tr.track['trail_y'][idx]
    ax_x, ax_y = shapeplot._project_xy(tr.track['ax_x'][idx], tr.track['ax_y'][idx])
    for i, _  in enumerate(idx):
        yield (x[i],y[i]), (tx[i], ty[i]), (ax_x[i], ax_y[i])

def exp_shape_iter(tr, idx):
    x, y = tr.track['x'][idx], tr.track['y'][idx]
    ax_x, ax_y = shapeplot._project_xy(tr.track['ax_x'][idx], tr.track['ax_y'][idx])
    length = tr.track['length'][idx]
    width = tr.track['width'][idx]
    for i, _  in enumerate(idx):
        yield (x[i],y[i]), (ax_x[i], ax_y[i]), length[i], width[i]


def _surface_annotation(ax):
    args = parameters.thisread()
    # surface annotation
    if args.surface.shape == 'infsteps':
        twstep.draw_step_lines(ax, ax.get_xlim())
    elif args.surface.shape == 'hexspheregrid':
        # condition to avoid too many small circles
        usepatches = True if args.surface.sphere_radius > 1. else False
        shapeplot.draw_hex_grid(ax, args, ax.get_xlim(), ax.get_ylim(), usepatches)
    elif args.surface.shape == 'plane':
        xlim, ylim = ax.get_xlim(), ax.get_ylim()
        xl, xr = xlim
        yl, yr = ylim
        dl = 5
        xmin, xmax = xl*np.floor(xl/dl), xr*np.ceil(xr/dl)+dl/2
        ymin, ymax = yl*np.floor(yl/dl), yr*np.ceil(yr/dl)+dl/2,
        xspace = np.linspace(xmin, xmax, (int(abs(xmax-xmin))//dl)+1, True)
        yspace = np.linspace(ymin, ymax, (int(abs(ymax-ymin))//dl)+1, True)
        style = {"alpha":0.4,"linewidth":1,'color':'k'}
        for i, x in enumerate(xspace):
            ax.plot([x,x],[yspace[0],yspace[-1]], **style)
        for i, y in enumerate(yspace):
            ax.plot([xspace[0],xspace[-1]],[y,y], **style)

# animate the body 
def outline(fig, trs, sample=1, camera='fixed', savefile=False, fps=None):
    ax = fig.gca()
    citer = itertools.cycle(colors)

    print("Constructing Iterators")
    # get the longest track
    ltrack = trs[ np.argmax(np.array([tr.track['time'][-1] for tr in trs])) ]
    mxtime = ltrack.track['time'][-1]
    # step = int(round(float(sample)/float(mxtime) * ltrack.size))
    # use idx array to step the Animation
    idx = np.array(np.arange(0, ltrack.size-1, sample), dtype=int)
    # print(idx)

    # Need to set axis limits in advance
    print("Computing Axes Limits")
    if camera=='follow':
        vsize = 10.
        x0, y0, _ = trs[0].get_origin()
        xlims = (x0-vsize/2,x0+vsize/2)
        ylims = (y0-vsize/2,y0+vsize/2)
    elif camera=='fixed':
        xlims = readtrack.trxmin(trs, 'x'), readtrack.trxmax(trs, 'x')
        ylims = readtrack.trxmin(trs, 'y'), readtrack.trxmax(trs, 'y')
    else:
        print("undefined camera type ", camera)
        sys.exit()
    # set initial axis limits
    ax.set_aspect('equal')
    ax.set_xlim(*xlims)
    ax.set_ylim(*ylims)
    

    # linedata artists:
    lstyle = {'linewidth':4, 'linestyle':'--', 'alpha':0.6, 'color':'blue'}
    trstyle = {'linewidth':4, 'linestyle':'--', 'alpha':0.6, 'color':'red'}
    ld_lineartists = [ax.plot([], [], **lstyle)[0] for tr in trs]
    tr_lineartists = [ax.plot([], [], **trstyle)[0] for tr in trs]

    # body artist
    data_getter = [shape_iter(tr, idx) for tr in trs]
    cellartists = [CellArtist.init(ax, color=next(citer)) for i, _ in enumerate(trs)]

    # timer
    textartist = ax.text(0.10,0.90,"{:04.1f}s".format(0), fontsize=36, transform=ax.transAxes)

    #
    _surface_annotation(ax)
    """Need a CellArtist for every cell in the animation"""
    # setup
    cxy = trs[0].get_cxy()
    # driving function for the animation
    def update_shapes(n, tq):
        tq.update(1)
        artists = []
        time = float(n)/10
        if camera=='follow':
            xn, yn = cxy[n]
            ax.set_xlim(xn-vsize/2,xn+vsize/2)
            ax.set_ylim(yn-vsize/2,yn+vsize/2)
        textartist.set_text("{:04.1f}s".format(time))
        for i, cartist in enumerate(cellartists):
            artists.extend( cartist.set_data(next(data_getter[i])) )
        for i, lartist in enumerate(ld_lineartists):
            tr = trs[i]
            lartist.set_data(tr['x'][:n], tr['y'][:n])
            artists.append( lartist )
        for i, trartist in enumerate(tr_lineartists):
            tr = trs[i]
            trartist.set_data(tr['trail_x'][:n], tr['trail_y'][:n])
            artists.append( trartist )

        return artists
    
    print('num. frames ', idx.size)
    # inter = 10 if savefile else 2
    inter = sample
    tq = tqdm(total=idx.size)
    track_ani = animation.FuncAnimation(fig, update_shapes, iter(idx),
            interval=inter, blit=False, repeat=False, save_count=idx.size, fargs=[tq])

    if savefile:
        print('saving animation to {}'.format(savefile))
        track_ani.save(savefile, fps=fps)
    else:
        print("Display")
        plt.show()
    tq.close()

# animate the line-track
# plot track using scatter

def ani_longtracks(trs, to_save=False):
    # trs = rt.trackset()
    fig = plt.gcf()
    ntracks = len(trs)

    step = 10
    grid = (2, 3)

    print("Creat Grid Layout")
    ax = plt.subplot2grid(grid, (0,0), rowspan=2)
    rx = plt.subplot2grid(grid, (0,1), colspan=2)
    #zx = plt.subplot2grid(grid, (1,1), colspan=2)

    ## setup right axes 
    rx.set_xlim(0, readtrack.maxtime(trs))
    rx.set_ylim(0, np.pi)
    # rx.set_ylabel(r'$\theta_{\text{tilt}}$')
    rx.set_ylabel(r'$\theta$')
    print("Computing tilt angles")
    ttilt = list(map(astrack.tiltangle, readtrack.nptracks(trs)))
    tiltlines = [rx.plot([], [], '-')[0] for _ in range(ntracks)]

    print("Setting up Axes")
    ## setup main axis
    ax.set_aspect('equal')
    delta = 0.5
    xlim = readtrack.trmin(trs, 'x')-delta, readtrack.trmax(trs, 'x')+delta
    ylim = readtrack.trmin(trs, 'y')-delta, readtrack.trmax(trs, 'y')+delta
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    # ax.plot() returns list of line objects but we create one so take [0] element
    lines = [ax.plot([], [], '-')[0] for _ in range(ntracks)]

    plt.tight_layout()

    """
    # construct iterators
    def make_tilt_iter(tr, tiltline):
        return (line.set_data(tr.track['time'][:n], ttilt[:n]) for n in range(0, tr.size,step))
    def make_iter(tr, line):
        return (line.set_data(tr.track['x'][:n], tr.track['y'][:n]) for n in )
    iterators = [make_iter(tr, line) for tr, line in zip(trs, lines)]
    t_iterators = [make_tilt_iter(tr, line) for tr, line in zip(trs, tiltlines)]
    applines = lines + tiltlines
    appiters = iterators + t_iterators
    """

    def update_track(n):
        for tr, line in zip(trs, lines):
            line.set_data(tr.track['x'][:n], tr.track['y'][:n])
        for tr, tilt, line in zip(trs, ttilt, tiltlines):
            line.set_data(tr.track['time'][:n], tilt[:n])
        return lines + tiltlines
        #return lines 

    N = int(max([tr.size for tr in trs])/step)
    print("Construct FuncAnimation")
    base = list(range(0,trs[0].size,step))
    track_ani = animation.FuncAnimation(fig, update_track, base,
            interval=5, blit=False, repeat=False)

    if to_save:
        savename = 'ani_lt.mp4'
        print('saving animation to {}'.format(savename))
        track_ani.save(savename)
    print("Display")
    plt.show()

