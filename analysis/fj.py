#!/usr/bin/env python3

# the analysis code for FJ data

import sys, os
import numpy as np
import collections

import datapaths
import readmat
import matdef
import command
#import shapeplot
import twanimation

import matplotlib.pyplot as plt
import scipy.stats
from tqdm import tqdm

import _fj


"""
Analysis of raw FanJin data with goal to see whether walking <-> crawling transitions
can be detected by aspect ratio. 
"""


def anycol(col, idx):
    ft = _fj.npyload(idx)
    ft.compute_aspect()
    coldata = ft.get_col(ft.tracks[0], col)
    print('mean/var', np.mean(coldata), np.var(coldata))
    plt.plot(coldata)
    plt.show()

#################################################################################
# try step detection

# using methods from
#Generalized methods and solvers for noise removal from piecewise constant signals. II. New methods
sys.path.append('/home/dan/tools/step_detection')

# robust jump penalization

# from pwc_medfiltit import pwc_medfiltit
# from pwc_cluster import pwc_cluster
# from pwc_jumppenalty import pwc_jumppenalty
# from pwc_bilateral import pwc_bilateral
# from pwc_tvdip import pwc_tvdip

# from pwc_tvdrobust import pwc_tvdrobust
# def tvdr():
#     ft = _fj.npyload(1279)
#     width = ft.get_col(ft.tracks[0], 'width')
#     length = ft.get_col(ft.tracks[0], 'length')
#     #result = pwc_tvdrobust(width, lamb=100.)
#     #result = pwc_bilateral(width)
#     result = pwc_jumppenalty(width, square=False)
#     plt.plot(width)
#     plt.plot(result)
#     result = pwc_jumppenalty(length, square=False)
#     plt.plot(length)
#     plt.plot(result)

#     plt.show()

# simple convolution method

def step_detect():
    ft = _fj.npyload(1279)
    width = ft.get_col(ft.tracks[0], 'width')
    length = ft.get_col(ft.tracks[0], 'length')
    #_step_detect(width)
    _step_detect(length)
    plt.show()


#https://stackoverflow.com/questions/48000663/step-detection-in-one-dimensional-data
def _step_detect(wh):
    wh -= np.average(wh)
    up_step = np.hstack((np.ones(len(wh)), -1*np.ones(len(wh))))
    down_step = -up_step
    wh_upstep = np.convolve(wh, up_step, mode='valid')
    wh_downstep = np.convolve(wh, down_step, mode='valid')

    # get the peak of the convolution, its index
    # cleaner than np.where(dary_step == dary_step.max())[0][0]
    plt.plot(wh_upstep, linestyle='--')
    plt.plot(wh_downstep, linestyle='--')
    return
    upstep_indx = np.argmax(wh_upstep)  
    downstep_indx = np.argmax(wh_downstep)  
    # plots
    plt.plot(wh)
    plt.axvline(upstep_indx, c='b')
    plt.axvline(downstep_indx, c='r')


###################################################################################

def order_width_var():
    _order_width_var(_fj.npyloadall())

def _order_width_var(ft):
    print("getting width data")
    width = np.array(ft.get_columns(['width']), dtype=object)
    # compute variances
    vwidth = list(map(np.var, width))
    # sort
    wsort = sorted(enumerate(vwidth), key=lambda t:t[1])
    # get key and value lists
    dwsort = collections.OrderedDict(wsort)
    keys, vals = list(dwsort.keys()), list(dwsort.values())
    #pof = scipy.stats.percentileofscore(vals, 0.01) # 92%
    cutoff = 0.002
    ok = np.argmax(np.array(vals)>cutoff)
    goodidx = np.array(sorted(dict(wsort[:ok]).keys()), dtype=int)
    slicefile= 'fjslice/width_ok.npy'
    print("writing to", slicefile)
    np.savetxt(slicefile, goodidx.astype(int))

    #_wv_plot_everything(width, keys, vals)

def _wv_plot_everything(width, keys, vals):
    # plot everything
    def imap(i): return keys[i]
    print("Constructing defaultsaveall function")
    wh_iter = ((i, wh) for i, wh in enumerate(width[keys]))
    @command.defaultsaveall(wh_iter, imap)
    def ordered_vwidth(wh):
        plt.ylim(0.5, 2.5)
        plt.title('var(width) = {}'.format(wh.var()))
        plt.plot(wh)
    
    ordered_vwidth()

###################################################################################
# filtering crawling/walking data

# this is superceded by fjanalysis.py 

def aspect_switching():
    whs = np.load('plots/data/whs.npy')
    _aspect_switching(whs)

def coarse_graining():
    whs = np.load('plots/data/whs.npy')
    goodidx = np.load('plots/data/whs.meta.npy')
    framewindow = 200
    smwhs = _coarse_graining(framewindow, whs, goodidx)
    _cg_plot_everything(framewindow, whs, smwhs, goodidx)

def _coarse_graining(framewindow, whs, goodidx):
    print('computing coarse grained tracks with framewindow = {}'.format(framewindow))
    frh = framewindow/2
    smwhs =[]
    for i, wh in tqdm(enumerate(whs)):
        sbasis = np.arange(frh, wh.size-frh, 1)
        smoothwh = np.array([np.mean(wh[j-frh:j+frh-1]) for j in sbasis])
        smwhs.append(smoothwh)

    return smwhs

def _cg_plot_everything(framewindow, whs, smwhs, goodidx):
    frh = framewindow/2
    ax = plt.gca()
    formdir = 'plots/fj_coarse'
    datapaths.force_mkdir(formdir)
    form = os.path.join(formdir, 'coarse_wh_{:05d}.png')

    for i, idx in enumerate(goodidx):
        wh = whs[i]
        smoothwh = smwhs[i]
        # plot side by side
        ax.clear()
        ax.plot(wh)
        sbasis = np.arange(frh, wh.size-frh, 1)
        ax.plot(sbasis, smoothwh)
        print('writing to', form.format(idx))
        plt.savefig(form.format(idx))


def _aspect_switching(whs):
    # look for tracks with aspect ratio switches and plot the ratio vs time to plots/fj_switching
    cutoff = 1.6
    # coarse graining window ( in no. frames)
    framewindow = 100
    """
    I think switching happens over a period of ~100s frames, 10s seconds
    This is consistent with the retraction rate of pili.
    """
    # look for changes larger than delta
    delta = 0.5
    def search_switch(wh):
        """
        Should combine several measures to build up a confidence value in having found a
        walking <-> crawling transition
        """
        # thin to window
        diff = np.abs(wh[framewindow:] - wh[:-framewindow])
        didx = diff > delta
        _didx = np.argwhere(didx).ravel()
        if _didx.size is 0:
            return _didx
        idx = np.argwhere(np.logical_or(wh[_didx] < cutoff, wh[_didx+framewindow] < cutoff)).ravel()
        _didx = _didx[idx]
        _didx = _fj._discard_consecutive(_didx)
        return _didx

    alldir = 'plots/fj_switching/'
    datapaths.force_mkdir(alldir)
    form = os.path.join(alldir, 'track_{:05d}.png')
    sw_whs = list(map(search_switch, whs))
    swidx = [(i, sw) for i,sw in enumerate(sw_whs) if sw.size>0]
    for (i, sw) in swidx:
        plt.clf()
        plt.plot(whs[i])
        for s in sw:
            plt.axvline(s, lw=0.5, c='k', linestyle='--')
        target = form.format(i)
        print("writing file to ", target)
        plt.savefig(target)


@command.defaultsave(True)
def fj_mean_aspect(whs):
    means = np.array([np.mean(wh) for wh in whs])
    kde = scipy.stats.gaussian_kde(means)
    mspace = np.linspace(means.min(), means.max(), 100)
    plt.plot(mspace, list(map(kde.evaluate, mspace)))

###################################################################################
# animation

def animate_one(idx, sample=10):
    track = _fj.trackload([idx])[0]
    outline([track], sample)

def outline(trs, sample=10):
    twanimation.outline(plt.gcf(), trs, sample)

if __name__=='__main__':
    ff, thisargs = command.process(locals())
    ff(*thisargs)
