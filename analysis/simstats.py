
"""
fjstats and yrstats contain plotting functions for exp data
here we apply that too simulated data
"""

import os, sys
import numpy as np
import scipy.stats
import command, filesystem
import plotutils

import _analysefj

import readtrack
import readmat
import fjstats
import readyr
import yrstats 

import matplotlib.pyplot as plt


### Compare with simulated data

pdir = os.path.join(os.path.dirname(__file__), 'simplots/')
filesystem.safemkdir(pdir)

simdir = "/data/dan/run/"
# 
# This simulation file has kmsd close to experiment 
simdirs = {}
simdirs['kmsdsim'] = "/data/dan/run/dec_18/plane_set_parameters/crawling_4x4x5/tau_00.500_pilivar_00.50000_free_tau_eovers_0001.000/data"

simdirs['highpvar'] = os.path.join(simdir, "dec_18/plane_set_parameters/crawling_4x4x5/tau_02.000_pilivar_20.00000_free_tau_eovers_0001.000/data")

simdirs['highpersist'] = "/data/dan/run/dec_18/plane_set_parameters/crawling_4x4x5/tau_00.500_pilivar_20.00000_free_tau_eovers_0001.000/data"

simdirs['rangle_by_eye'] = "/data/dan/run/dec_18/plane_set_parameters/crawling_4x4x5/tau_00.500_pilivar_00.50000_free_tau_eovers_0000.800/data"

def load_simulated(dset='kmsdsim'):
    ddir = simdirs[dset]
    trs  = readtrack.trackset(ddir=ddir)
    st = readmat.ExpTracks.usesimulated(trs)
    return st

def sim_steptime(st, step_d=0.4, scaling=1., fit=False, linekw={}):
    step_d = step_d
    fjstats.fj_steptime(st, step_d, scaling, fit, linekw)
    plt.title('')
    plt.savefig(os.path.join(pdir, '_stime_distribution.png'))

def _steptime(st, step_d=0.4):
    plt.clf()
    meta = fjstats.fj_steptime(st, step_d)
    out = 'plots/stime/_stime_distribution.png'
    filesystem.safemkdir(out)
    plt.savefig(out)
    command.writemeta(meta, out)

def sim_fjr_corr(st, form):
    st.compute_orientation()
    # yaxis limits
    axislims = [(None, None)]*3 
    yrstats.fjr_corr(st, form, axislims=axislims)


#### comparing distributions

# use ksmeasure compare deviations angle for YR and SIM data
@command.defaultsave(pdir=pdir)
def simtheta(st):
    fastpercent = 99
    sim_theta = yrstats.frtheta(st, fastpercent)
    np.nan_to_num(sim_theta, copy=False)
    linkw = {}
    histstyle = {'density':True, 'bins':30, 'alpha':0.4}
    plt.hist(sim_theta, **histstyle)
    plotutils.kdeplot(sim_theta, linkw)

def ks_against_yr(st):
    fastpercent = 95
    sim_theta = yrstats.frtheta(st, fastpercent)
    np.nan_to_num(sim_theta, copy=False)
    # 
    yt = readyr.inithres()
    yr_theta = yrstats.frtheta(yt, fastpercent)

    linkw = {'label':'yr'}
    plotutils.kdeplot(yr_theta, linkw)
    linkw = {'label':'sim'}
    plotutils.kdeplot(sim_theta, linkw)
    plt.legend()
    plt.show()
    
    ksd, pv = scipy.stats.ks_2samp(sim_theta, yr_theta)
    print('ksd', ksd)


#

def wlsq_samples(X, Y, rescale=False):
    xout = plotutils.kdeplot(X)
    xspace, xpde = xout['space'], xout['pde']
    yspace, ypde = yout['space'], yout['pde']
    # stub ...

def wlsq(A, B):
    # weighted sum of squares between arrays A and B
    # weights could be A or B for similar distributions, lets take A+B/2
    #omega = np.sqrt((A+B)/2.)
    omega = np.sqrt(A)
    return 1./np.sum(omega) * np.sum(omega*(A-B)**2)



@command.defaultsave(pdir=pdir)
def divider_dimension(ft):
    #dbasis = np.linspace(0.03, 0.6, 21, True)
    dbasis = np.logspace(-2, 1, 31, True)
    dbasis, divider_count = fjstats._divider_dimension(ft, dbasis)
    cut = slice(10,20)

    norm_dim, intercept = np.polyfit(dbasis[cut], divider_count[cut], 1)
    plotutils.abline(plt.gca(), norm_dim, intercept)
    plt.title(r'{} D = {:.04f}'.format(ft.source, norm_dim))
    print('norm divider dimension', norm_dim) 

    grad2 = fjstats.finite_grad(dbasis, divider_count)
    plt.plot(dbasis, grad2)
    print('2nd div error', np.mean(np.abs(grad2[cut])))
    #grad2 = fjstats.convolve_grad(dbasis, divider_count)
    #plt.plot(dbasis[:-1], grad2)


    dl = _analysefj.disp_length(ft)
    xguide = (np.min(dl)/10.,)
    for x in xguide:
        plt.axvline(np.log(x), alpha=0.6)
    
    print('end to end distance', np.min(dl), np.max(dl))



if __name__=='__main__':

    st = load_simulated('rangle_by_eye')
    #sim_steptime(st)
    # divider dimension
    #divider_dimension(st)
    #
    simtheta(st)

    # plot time vs. various quantities
    #form = os.path.join(pdir, 'singletrack/r_angle_{:04d}.png')
    #form = os.path.join(pdir, 'singletrack/highpvar/angle_{:04d}.png')
    #form = os.path.join(pdir, 'singletrack/highpersist/angle_{:04d}.png')
    #sim_fjr_corr(st, form)


    # comparing R angle distributions for YR and Sim data
    #ks_against_yr()


