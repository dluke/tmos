
import os
import itertools

from command import defaultsave

import matplotlib as mpl
from matplotlib import pyplot as plt
import numpy as np

import readtrack

from collections import OrderedDict
import scipy.stats

import pili
mpl_stylelib_path = os.path.join(os.getenv('PILI_ROOT'), 'mpl/')
def get_style(stylename):
    stylefile = stylename
    if not stylefile.endswith('.mplstyle'):
        stylefile += '.mplstyle'
    return os.path.join(mpl_stylelib_path, stylefile)

def default_style():
    plt.style.use(['fivethirtyeight', get_style('ft8')])

# useful 
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
prop_cycle_color = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

# write to stdout
def print_dict(sd, cwidth=40):
    form = '{:<%ds} {}' % cwidth
    for key, val in sd.items():
        print(form.format(key, val))

# plot linearised velocity appropriately
def plot_lvel(ax, tr, xcut=None):
    dt = tr.get_step_dt()
    lvel = tr.get_step_speed()
    time = np.cumsum(dt)
    if xcut:
        icut = np.searchsorted(time, xcut)
        time = time[:icut]
        lvel = lvel[:icut]
    l =  ax.plot(np.repeat(time,2)[1:], np.repeat(lvel,2)[:-1])
    return l



# plot a dictionary as a bar graph
@defaultsave()
def plot_bar(ax, dct):
    i = 0
    names = []
    for name, height in dct.items():
        ax.bar(i, height)
        i += 1
        names.append(name)
    ax.set_xticks(range(len(dct)))
    ax.set_xticklabels(names)

#
def kde_eval(stat, res):
    xn = np.min(stat)
    xm = np.max(stat)
    kde = scipy.stats.gaussian_kde(stat)
    mspace = np.linspace(xn, xm, res)
    pde = kde.evaluate(mspace)
    return mspace, pde

#
verbose = False
def ax_kdeplot(ax, stat, linekw={},  res=200, xlims=(None,None), hist=False):
    xn, xm = xlims
    if xn is None:
        xn = np.min(stat)
    if xm is None:
        xm = np.max(stat)
    if verbose:
        print("Computing pdf with {} data points at resolution {}".format(stat.size, res))
    kde = scipy.stats.gaussian_kde(stat)
    mspace = np.linspace(xn, xm, res)
    pde = kde.evaluate(mspace)
    p, = ax.plot(mspace, pde, **linekw)
    # hist = False
    if hist:
        histstyle = {'alpha':0.25}
        ax.hist(stat, bins=res//2, range=(xn,xm), density=True, color=p.get_color(), **histstyle)
    outd = OrderedDict()
    # should return plot handle instead of this dictionary
    outd['space'] = mspace
    outd['pde'] = pde
    outd['handle'] = p
    return outd


def kdeplot(stat, linekw={},  res=200, xlims=(None,None)):
    return ax_kdeplot(plt.gca(), stat, linekw, res, xlims)

#################################################################################

def fix_tex_label(st):
    return st.replace('_', '\_')

def mid90(arr):
    """cut the middle 90 values out of a sorted array"""
    return np.clip(arr, np.percentile(arr, 5), np.percentile(arr, 95))

@defaultsave()
def kmsd_vdiff_scatter(grads, meanvdiff, varvdiff):
    plt.scatter(grads, mid90(meanvdiff), c=mid90(varvdiff))
    plt.xlabel('k_msd')
    plt.ylabel('<|vlead| - |vtrail|>')

from operator import sub
def get_aspect(ax):
    # Total figure size
    figW, figH = ax.get_figure().get_size_inches()
    # Axis size on figure
    _, _, w, h = ax.get_position().bounds
    # Ratio of display units
    disp_ratio = (figH * h) / (figW * w)
    # Ratio of data units
    # Negative over negative because of the order of subtraction
    data_ratio = sub(*ax.get_ylim()) / sub(*ax.get_xlim())
    return disp_ratio / data_ratio


def abline(axes, slope, intercept):
    """Plot a line from slope and intercept"""
    x_vals = np.array(axes.get_xlim())
    y_vals = intercept + slope * x_vals
    plt.plot(x_vals, y_vals, '--')


#################################################################################


def ltall():
    """
    plot x,y tracks but write them all to individual image files
    """
    trs = readtrack.trackset()
    for i,tf in enumerate(trs):
        track = tr.track
        plt.clf()
        plt.axes().set_aspect('equal')
        x = track['x']
        y = track['y']

        plt.plot(x, y, '-', c=colors[i%len(colors)])
        form = 'plots/single_track_{:02d}'.format(i)
        plt.savefig(form)



##################################################################################
# pretty ax parameters

def _short_form(value, tick_number):
    return "{:.1f}".format(value)

def prettyax_track(ax, trs):

    def _round_away(a,b):
        return np.floor(a)-0.1, np.ceil(b)+0.1
    left, right = _round_away(*ax.get_xlim())
    ax.set_xticks([left+0.1,right-0.1])
    ax.set_xlim(left, right)
    left, right = _round_away(*ax.get_ylim())
    ax.set_yticks([left+0.1,right-0.1])
    ax.set_ylim(left, right)

    ax.xaxis.set_major_formatter(plt.FuncFormatter(_short_form))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(_short_form))
    ax.tick_params(axis='both', which='major', labelsize=16)
    # ax.tick_params(axis='both', which='minor', labelsize=16)
    return ax

def invisax(ax):
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
    return ax


#################################################################################
# generate a pill geometry in 2d

e_x = np.array([1,0])
e_y = np.array([0,1])
def pillgeometry(R, length, res=60):
    origin = np.array([0,0])
    pill = [origin - length/2. * e_x - R*e_y]
    right_semi = np.vstack(
            [np.array([np.cos(theta), np.sin(theta)])
        for theta in np.linspace(-np.pi/2, np.pi/2, res+1, True)]
        )
    right_semi_at = (R*right_semi) + (length/2. * e_x).T
    left_semi = np.vstack(
            [np.array([np.cos(theta), np.sin(theta)])
        for theta in np.linspace(np.pi/2, 3*np.pi/2, res+1, True)]
        )
    left_semi_at =  (R*left_semi) - (length/2. * e_x).T
    fpt = right_semi_at[0].reshape((1,2))
    return np.concatenate([right_semi_at, left_semi_at, fpt], axis=0)


#################################################################################
# for plotting heatmaps of high dimensional data projected onto 2d

class KernelEstimator(object):
# with reference to
# https://github.com/mwaskom/seaborn/blob/b13dfdb6570b943391f5c230b038b44d600a19a7/seaborn/external/kde.py#L82
    def __init__(self, xy, value, bw_method="scott"):
        self.xy = xy
        self.value = value
        self.d = 2
        nx, ny = self.xy.shape[:2]
        self.neff = nx * ny
        self.set_bandwidth(bw_method=bw_method)
    
    def evaluate(self, points):
        n, d = points.shape
        result = np.zeros(n)
        for i, pt in enumerate(points):
            diff = self.xy - pt
            kern = np.exp( -np.sum(diff*diff,axis=1)/(self.bw * 2.0) ) / self.bw
            denom = np.sum(kern)
            if denom == 0.:
                result[i] == 0.
            else:
                result[i] = np.sum(self.value * kern) / denom
        return result

    def set_bandwidth(self, bw_method):
        if bw_method == "scott":
            scott = np.power(self.neff, -1./(self.d+4))
            bw = scott
        elif np.isscalar(bw_method) and not isinstance(bw_method, str):
            bw = bw_method
        self.bw = bw
        return self.bw





