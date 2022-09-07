

# these methods are a mess, go through all of them, todo

# the analysis code for FJ data

from collections import OrderedDict
from glob import glob
import sys, os, random

import numpy as np
import scipy.stats

import plotutils # formerly as utils
import twutils

import readmat
import matdef
import shapeplot

import command
from command import defaultsave
from matplotlib import pyplot as plt
import matplotlib



##########################################################################
# basic analysis

def track_disps(ft, timestep, cols=matdef.CENTROID):
    acxy = ft.get_columns(cols, timestep=timestep)
    return [np.linalg.norm(cxy[:,1:] - cxy[:,:-1], axis=0) for cxy in acxy]

def disps(ft, timestep=None, cols=matdef.CENTROID):
    if timestep is None:
        timestep = ft.timestep
    disp = np.concatenate(track_disps(ft, timestep, cols))
    return disp

def meanvel(ft, cols=matdef.LEAD):
    avgvel = np.mean(disps(ft, cols=cols)/ft.timestep)
    return avgvel

# end to end distance of the track
def disp_length(ft, cols=matdef.LEAD):
    tracks = ft.get_columns(cols)
    def displ(track):
        x, y = track
        xd = np.array([x[-1] - x[0], y[-1] - y[0]])
        return np.linalg.norm(xd)
    return list(map(displ, tracks))

def track_vel_dot_ax(Ftracking, timestep, point=matdef.CENTROID):
    """Reorientations. compute angles between displacements and body axis"""
    if 'ax_x' not in Ftracking.colmap:
        Ftracking.compute_ax()
    if point==matdef.CENTROID and 'center_x' not in Ftracking.colmap:
        Ftracking.compute_cxy()

    def dispv(track):
        cx, cy = track
        dnorm = np.sqrt( (cx[1:]- cx[:-1])**2 + (cy[1:] - cy[:-1])**2 )
        vstack = np.stack(
                [(cx[1:]- cx[:-1])/dnorm,
                (cy[1:] - cy[:-1])/dnorm],
                axis=1)
        return vstack
    tracks = Ftracking.get_columns(point, timestep=timestep)
    allaxes = Ftracking.get_columns(['ax_x', 'ax_y'], timestep=timestep)
    ax_x = Ftracking.get_whole_col('ax_x')
    ax_y = Ftracking.get_whole_col('ax_y')

    pertrack = []
    for track, axes in zip(tracks, allaxes):
        dv = dispv(track)
        axblock = axes.T[1:]
        costh = np.sum(dv * axblock, axis=1)
        theta = np.arccos(costh)
        #np.nan_to_num(theta, copy=False) # didn't fix problem with nan values?
        theta = np.sign(np.cross(dv, axblock)) * theta
        pertrack.append(theta)
    return pertrack

def vel_dot_ax(Ftracking, timestep, point=matdef.CENTROID):
    return np.concatenate(track_vel_dot_ax(Ftracking, timestep, point))


def plot_disps(Ftracking, extra=''):
    Ftracking.filter_short(mtime=100.)
    #Ftracking.filter_by_method(filter_kmsd_1)
    nbins = 40
    stepv = [1, 10, 30]
    upper_xlim = [0.5, 2.0, 2.0]

    fig, axes = plt.subplots(1,len(stepv))
    fig.set_size_inches(18, 8)
    for i, (ax, step) in enumerate(zip(axes, stepv)):
        dsp = disps(Ftracking, timestep=step)
        """normalise bins to represent probabilities"""
        weights = np.full(len(dsp), 1./len(dsp))

        ax.hist(dsp, range=(0,upper_xlim[i]), bins=nbins, weights=weights)
        #ax.hist(dsp, bins=nbins, weights=weights)

        ax.set_xlabel(r'$\Delta r$ ($\mu m$) in {:n} sec'.format(step),
                fontsize = 24)
        ax.set_ylabel('probability')
        ax.set_ylim(0, 0.6)
        ax.axvline(0.4, c='r')
    tosave = os.path.join('plots/', '_'.join(['disps', extra]))
    tosave += '.png'
    print('saving to ', tosave)
    #plt.suptitle("Displacments distribution for original data.")
    plt.tight_layout()
    command.default_result_dirs()
    plt.savefig(tosave)

##########################################################################
# the methods we actually use for analysing reorientation distribution

def frmodes(Ftracking, fastpercent, timestep=0.1, col=matdef.CENTROID): 
    # compute the velocity of the head
    pu = disps(Ftracking, timestep=timestep, cols=col)/timestep

    # compute FanJin reorientation angles
    thetas = vel_dot_ax(Ftracking, timestep=timestep, point=col)

    fast = np.percentile(pu, fastpercent)
    lines = linehistogram(thetas, fast, pu)
    return (fast, lines)

# thetas - an array of angles
# fast - a fast speed, calculate with np.percentile
# pu - an array of speeds/displacements
def linehistogram(thetas, fast, pu):
    print('thetas.size', thetas.size)
    finite = np.isfinite(thetas)
    n_nonfinite = np.count_nonzero(~np.isfinite(thetas))
    # why to we get infinite values
    print('portion of infinite values of theta {}/{}'.format(n_nonfinite, pu.size))

    # exclude infinite
    pu = pu[finite]
    thetas= thetas[finite]

    fast_th= thetas[pu>fast]
    slow_th= thetas[pu<fast]

    nbins = 36
    bins = np.linspace(-np.pi,np.pi,nbins+1,True)
    bin_centres = bins[:-1] + (bins[1]-bins[0])/2.

    fastline, _ = np.histogram(fast_th, bins=bins, density=True)
    if not slow_th.size:
        slowpair = None
    else:
        slowline, _ = np.histogram(slow_th, bins=bins, density=True)
        slowpair = (bin_centres, slowline)
    fastpair = (bin_centres, fastline)
    return fastpair, slowpair


########################################
# The plotting we actually use for simulated data

# plotting method for grids of polar reorientation angle plots
def polar(ax, fast_lines, leg_fontsize=5, svalue=''):
    (fast, lines) = fast_lines 
    fastpair, slowpair = lines
    config = {'width':0.15, 'alpha':0.6}
    if fastpair:
        bin_centres, fastline = fastpair
        ax.bar(bin_centres, fastline, color='b', label='$> {:3.2f}$'.format(fast), 
                **config)
    if slowpair:
        bin_centres, slowline = slowpair
        ax.bar(bin_centres, slowline, color='r', label='$< {:3.2f}$'.format(fast),
                **config)

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_yticks([0.2, 0.5])
    ax.set_title(svalue)

    if slowpair and fastpair:
        ax.legend(fontsize=leg_fontsize)
        #ax.legend(fontsize=leg_fontsize, loc='upper center', bbox_to_anchor=(0.5,-0.05))
    return ax

########################################## # plot the average distribution of reorientation angles 
def plot_frmodes(Ftracking, fpercent=95, timestep=0.1, extra=''):
    fast, lines = frmodes(Ftracking, fpercent, timestep)

    fig = plt.figure(figsize=(10,10))
    ax = plt.subplot(111, projection='polar')
    ax = polar(ax, (fast, lines), leg_fontsize=30)
    plt.tight_layout()

    tosave = 'plots/total_rmodes_' + Ftracking.source + '.svg'
    plt.savefig(tosave)

##########################################
# correlation functions

def corr_orient(exp):
    allaxes = exp.get_columns(['ax_x', 'ax_y'])
    ts = 1.
    maxtime = 1000.
    tau = np.arange(ts, maxtime*ts , ts, dtype=int)
    def corr_or(tau, axes):
        """return numpy array containing the correlation function for one track"""
        ax = axes.T
        corr = np.full_like(tau, np.nan, dtype=float)
        for i, tauv in enumerate(tau):
            step = exp.get_step(tauv)
            if step > ax.shape[0]-1:
                break
            dot = np.sum(ax[step:] * ax[:-step], axis=1)
            #nem_dot = np.abs(dot)
            #"""nematic like"""
            #if not np.any(np.isnan(dot)):
                #print i
            # why nan is possible?
            corr[i] = np.nanmean(dot)
            #corr[i] = np.mean(dot)
        return corr
    corrblock = np.vstack([corr_or(tau, axes) for axes in allaxes])
    # average over tracks
    return tau, np.mean(corrblock, axis=0)

@command.defaultsave(True)
def plot_corr_orient(exp):
    plt.clf()
    tau, corr = corr_orient(exp)
    #print corr
    assert not np.any(np.isnan(corr))
    plt.plot(tau, corr)
    #plt.loglog(tau, corr)


### simple

# what is the mean of R angle
def check_assumptions(ft):
    ft.override_orientation_using_lead_trail()
    thetas = vel_dot_ax(ft, timestep=ft.timestep, point=matdef.LEAD)
    mor = np.mean(thetas)
    def to_degrees(a):
        return 180./np.pi * a
    print('mean R angle', to_degrees(mor))



if __name__=='__main__':

    ft = init(trial, alldirs)
    # check_assumptions(ft)




