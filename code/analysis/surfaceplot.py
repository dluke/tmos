
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

# some relevant debugging/plotting functions for testing convergence

def plothistory(history):
    hsx, hsy, hsteta = np.transpose(history[50:])
    # should cut off the zero elements

    #
    x = np.arange(hsx.size) + 1.
    plt.plot(x, hsx, label='fx')
    plt.plot(x, hsy, label='fy')
    plt.plot(x, hsteta, label='ftheta')
    plt.legend()
    plt.show()

#https://stackoverflow.com/a/24567352/2144720
def circles(x, y, s, c='b', vmin=None, vmax=None, **kwargs):
    """
    Make a scatter of circles plot of x vs y, where x and y are sequence 
    like objects of the same lengths. The size of circles are in data scale.

    Parameters
    ----------
    x,y : scalar or array_like, shape (n, )
        Input data
    s : scalar or array_like, shape (n, ) 
        Radius of circle in data unit.
    c : color or sequence of color, optional, default : 'b'
        `c` can be a single color format string, or a sequence of color
        specifications of length `N`, or a sequence of `N` numbers to be
        mapped to colors using the `cmap` and `norm` specified via kwargs.
        Note that `c` should not be a single numeric RGB or RGBA sequence 
        because that is indistinguishable from an array of values
        to be colormapped. (If you insist, use `color` instead.)  
        `c` can be a 2-D array in which the rows are RGB or RGBA, however. 
    vmin, vmax : scalar, optional, default: None
        `vmin` and `vmax` are used in conjunction with `norm` to normalize
        luminance data.  If either are `None`, the min and max of the
        color array is used.
    kwargs : `~matplotlib.collections.Collection` properties
        Eg. alpha, edgecolor(ec), facecolor(fc), linewidth(lw), linestyle(ls), 
        norm, cmap, transform, etc.

    Returns
    -------
    paths : `~matplotlib.collections.PathCollection`

    Examples
    --------
    a = np.arange(11)
    circles(a, a, a*0.2, c=a, alpha=0.5, edgecolor='none')
    plt.colorbar()

    License
    --------
    This code is under [The BSD 3-Clause License]
    (http://opensource.org/licenses/BSD-3-Clause)
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle
    from matplotlib.collections import PatchCollection

    if np.isscalar(c):
        kwargs.setdefault('color', c)
        c = None
    if 'fc' in kwargs: kwargs.setdefault('facecolor', kwargs.pop('fc'))
    if 'ec' in kwargs: kwargs.setdefault('edgecolor', kwargs.pop('ec'))
    if 'ls' in kwargs: kwargs.setdefault('linestyle', kwargs.pop('ls'))
    if 'lw' in kwargs: kwargs.setdefault('linewidth', kwargs.pop('lw'))

    patches = [Circle((x_, y_), s_) for x_, y_, s_ in np.broadcast(x, y, s)]
    collection = PatchCollection(patches, **kwargs)
    if c is not None:
        collection.set_array(np.asarray(c))
        collection.set_clim(vmin, vmax)

    ax = plt.gca()
    ax.add_collection(collection)
    ax.autoscale_view()
    if c is not None:
        plt.sci(collection)
    return collection

#surface is a hexsurface
from tmos.vector3d import Vector3d
from tmos.surface import HexSurface

def spawn_hex_grid(origin, surface, n):
    ex = surface.ex
    ey = surface.ey
    ez = surface.ez
    initial = [ey, ex, -ey, ez, -ex]
    order = [ey, -ez, ex, -ey, ez, -ex]
    count = [1, 0, 1, 1, 1, 1] 
    pts = [Vector3d()] 
    trace = Vector3d()
    for init in initial:
        trace = trace + init
        pts.append(trace)
    ntotal = len(pts)
    i = 0
    def nxt(i): return 0 if i+1 == 6 else i+1
    while ntotal < n:
        count[i] += 1
        for _ in range(count[i]):
            trace = trace + order[i]
            pts.append(trace)
            ntotal += 1
        i = nxt(i)
    return pts

def to_xy(vls):
    x = []
    y = []
    for v in vls:
        x.append(v.x)
        y.append(v.y)
    return x, y

def plot_surface(ax, args):
    surface = HexSurface(args.R, args.r)
    # configuration
    #color
    fill = False

    # decide this based on the longtracks distance...
    n = 50

    hex_grid = spawn_hex_grid(Vector3d(), surface, n)
    x, y = to_xy(hex_grid)
    out = circles(x, y, surface.plane_r, c='g', alpha=0.5, fc='none')
    ax.add_collection(out)
    return ax
    
#############################################################

def test_plot_surface():
    ax = plt.gca()
    hs = HexSurface(1., 0.5)
    plt.axes().set_aspect('equal', 'datalim')
    ax = plot_surface(ax, hs)
    plt.show()

if __name__=='__main__':
    test_plot_surface()
    pass

