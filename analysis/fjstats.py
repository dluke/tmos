#!/usr/bin/env python3 

import os, sys
import copy
import numpy as np
import scipy.stats
import scipy.optimize
import matplotlib.pyplot as plt
import matplotlib as mpl
from collections import OrderedDict

import command 
import stats 
import twutils
import plotutils
import shapeplot
import astrack
import matdef
import _fj
import _analysefj

# dividier dimension + some other junk

#################################################################################
##
#https://stackoverflow.com/questions/47417986/using-scipy-gaussian-kernel-density-estimation-to-calculate-cdf-inverse


##
# distributions over the population but each track contributes one value

partial = '/home/dan/twitching/pypili/src/analysis/fjdata'
plots = '/home/dan/twitching/pypili/src/analysis/fjplots'
metafile = os.path.join(partial, 'meta.txt')


def process(statitem, plotsdir):
    name, stat = statitem
    print('Analysing Stat {} ...'.format(name))
    print(scipy.stats.describe(stat))
    # draw distribution
    kde = scipy.stats.gaussian_kde(stat)
    mspace = np.linspace(stat.min(), stat.max(), 100)
    plt.clf()
    sstat = np.sort(stat)
    pde = list(map(kde.evaluate, mspace))
    # cumulative distribution
    cdf = np.cumsum(pde) 
    cdf /= np.max(cdf)
    plot_pde(name, plotsdir, mspace, pde)
    plt.clf()
    plot_cdf(name, plotsdir, mspace, cdf)

def plot_pde(name, plotsdir, mspace, pde):
    # Add vertical guide lines
    #hlow, hhigh = sstat[int(0.05*stat.size)], sstat[int(0.95*stat.size)]
    # faintkw = {'alpha':0.2, 'linestyle':'--', 'color':'k'}
    #plt.axvline(hlow, **faintkw)
    #plt.axvline(hhigh, **faintkw)
    plt.plot(mspace, pde)
    plt.xlabel(name)
    plt.ylabel('smooth probability density')
    dst_out = os.path.join(plotsdir, name+'_dst.png')
    print('saving distribution to {}'.format(dst_out))
    plt.tight_layout()
    plt.savefig(dst_out)

def plot_cdf(name, plotsdir, mspace, cdf):
    plt.plot(mspace, cdf)
    plt.xlabel(name)
    plt.ylabel('cumulative probability')
    cdf_out = os.path.join(plotsdir, name+'_cdf.png')
    print('saving cdf to {}'.format(cdf_out))
    plt.tight_layout()
    plt.savefig(cdf_out)


# Apply process to all stats
def process_stats(plotsdir):
    stats = regetstats()

    for statitem in list(stats.items()):
        process(statitem, plotsdir)

# aspect crawling/walking

def cutcrawling():
    # by aspect ratio
    stats = regetstats()
    aspratio = stats['aspect_ratio']
    kde = scipy.stats.gaussian_kde(aspratio)
    res = scipy.optimize.minimize(kde, x0=2.)
    crawlingidx = np.argwhere(aspratio > res.x).astype(int)
    print('of {} keep {}'.format(aspratio.size, float(crawlingidx.size)/aspratio.size))
    cba = os.path.join('fjslice/', 'crawling_by_aspect.npy')
    np.savetxt(cba, crawlingidx)


# step time
# define step time as the time taken to move one distance of one step d
step_d = 0.12 # \mu m
arbcut = 98
def calc_timedist(ft, step_d=step_d):
    xy = ft.get_columns(matdef.LEAD) #  list of shape(Time,2)  arrays

    def distance(xy1,xy2):
        return np.linalg.norm(xy2-xy1)

    def compute_step_times(trxy):
        trxy = trxy.T # rows are now points
        interval = 1
        step = trxy[0]
        stimes = [] # result
        for i, cv in enumerate(trxy):
            d = distance(cv, step)
            #print 'distance', d
            if  d > step_d:
                stimes.append(interval)
                step = cv
                interval = 1 # reset
            else:
                # keep going until the distance travelled is sufficient
                interval += 1 
        return stimes
    times = list(map(compute_step_times, xy))

    timedist = np.concatenate(times)
    timedist = timedist.astype(float) * ft.timestep # convert frame no. to time
    return timedist

def steptime(timedist, linekw={}):
    # apply kdeplot
    vcut  = np.percentile(timedist, arbcut)
    timedist = timedist[timedist < vcut] # arbitrary cut
    outd = plotutils.deplot(timedist, linekw, xlims=(None,vcut))
    
    histstyle = {'density':True, 'bins':30, 'alpha':0.4}
    plt.hist(timedist, **histstyle)

    tauspace, pde = outd['space'], outd['pde']
    return tauspace, pde

def _curvefit(tauspace, pde):
    #p0= None
    # fit to Aexp(tau/tau_p_1) + Bexp(tau/tau_p_2)
    #def f(xdata, A, tau_1, B, tau_2):
        #return A*np.exp(-xdata/tau_1) + B*np.exp(-xdata/tau_2)
    def f(xdata, tau_1):
        return (1./tau_1)*np.exp(-xdata/tau_1) 

    p0= None
    print('curve fit, with p0 = {}'.format(p0))
    #bounds = ([0,-np.inf,0,-np.inf],np.inf)
    bounds = (-np.inf, np.inf)
    sigma = 1./pde
    popt, pcov = scipy.optimize.curve_fit(f, tauspace, pde, p0=p0, bounds=bounds, sigma=sigma)
    tau_1 = popt
    print("tau_1 = ", tau_1)
    def func(xdata):
        return f(xdata, *popt)
    fitline = list(map(func, tauspace))
    plt.plot(tauspace, fitline)

#https://stackoverflow.com/questions/3433486/how-to-do-exponential-and-logarithmic-curve-fitting-in-python-i-found-only-poly
def _linfit(tauspace, pde):
    print('linear fit')
    # fit to Aexp(tau/tau_p)
    grad, logA = np.polyfit(tauspace, np.log(pde), 1, w=pde) # linear fit
    tau_p = 1./np.exp(logA)
    g_tau_p = -1./grad
    print('1/A = ', tau_p)
    print('tau_1 = ', g_tau_p)
    plot_expfit(tauspace, g_tau_p, tau_p)
    # return the fitting function
    return {'tau_p':tau_p, 'g_tau_p': g_tau_p} 

def plot_expfit(tauspace, g_tau_p, tau_p):
    def expform(tau):
        return 1/tau_p * np.exp(-tau/g_tau_p)
    y = list(map(expform, tauspace))
    style={'linestyle':'--'}
    plt.plot(tauspace,y, label='fit', **style)

def fj_steptime(ft, step_d, scaling, fit, linekw):
    timedist = calc_timedist(ft, step_d)/step_d * scaling
    tauspace, pde = steptime(timedist, linekw)#this method cuts the tail from the data
    print('No. tau', len(timedist))
    #histstyle = {'density':True, 'bins':10, 'alpha':0.6}
    #plt.hist(timedist, **histstyle)
    vcut  = np.percentile(timedist, arbcut)
    plt.xlim(0,vcut)
    if fit:
        meta = _linfit(tauspace, pde)
    else:
        meta = {}
    plt.xlabel(r"rescaled step time/d ")
    plt.ylabel(r"p(\tau)")
    plt.title(r"step distance d = {} (\mu m)".format(step_d))
    plt.legend()
    plt.tight_layout()

    # stop here
    return {}

    # else fit
    tau_p = meta['tau_p']; g_tau_p = meta['g_tau_p']
    def fitf(tau):
        return 1./tau_p * np.exp(-tau/g_tau_p)

    # ks measure
    fitline = list(map(fitf, tauspace))
    cumfit = scipy.integrate.cumtrapz(fitline, tauspace) 
    stau = (tauspace[1:]+tauspace[:-1])/2.
    cdf = scipy.interpolate.interp1d(stau, cumfit, fill_value='extrapolate')
    cumpde = scipy.integrate.cumtrapz(pde,tauspace)
    # rescale
    cumpde /= cumpde[-1]
    cumfit /= cumfit[-1]
    #plt.clf()
    #plt.plot(stau, cumfit)
    #plt.plot(stau, cumpde)
    #plt.show()
    D, pvalue = scipy.stats.kstest(timedist, cdf)
    print('D = ', D, 'what is stats.kstest doing?')
    ksd = np.max(np.abs(cumpde-cumfit))
    print('sup(abs(A-B)) ', ksd)
    meta.update({'ksd':ksd})
    return meta

def fj_steptime_fit(ft, fit=True):
    timedist = calc_timedist(ft)
    tauspace, pde = steptime(timedist)
    vcut  = np.percentile(timedist, arbcut)
    plt.xlim(0,vcut)

    if fit:
        cdf = _curvefit(tauspace, pde)

def stfits(ft, step_d=0.12, scaling=1., fit=True, linekw={}):
    meta = fj_steptime(ft, step_d, scaling, fit, linekw)
    meta.update({'step_d':step_d})

    plt.savefig(os.path.join(command.pdir, '_stime_distribution.png'))

    form = os.path.join(command.pdir, 'stime/_stime_distribution_{:.04f}.png')
    out = form.format(step_d)
    print('saving to ', out)
    plt.savefig(out)
    command.writemeta(meta, out)

def _r_div_dimension(ft):
    avgvel = np.mean(_analysefj.disps(ft, cols=matdef.LEAD)/ft.timestep)
    dbasis = np.linspace(0.03, 0.30, 11, True)
    def divcount(d_step):
        curlD = len(calc_timedist(ft, d_step))
        return np.log(curlD)
    divider_count = list(map(divcount, dbasis))
    log_dbasis = np.log(dbasis)
    norm_dim = np.mean(np.gradient(divider_count,log_dbasis))
    norm_dim = 1 - norm_dim
    print('dd', norm_dim)
    return norm_dim


def _divider_dimension(ft, dbasis):
    avgvel = np.mean(_analysefj.disps(ft, cols=matdef.LEAD)/ft.timestep)
    print('using dbasis', dbasis)

    def divcount(d_step):
        #curlD = len(calc_timedist(ft, d_step))
        curlD = np.mean(calc_timedist(ft, d_step))
        return np.log(curlD)
    divider_count = list(map(divcount, dbasis))

    dbasis = np.log(dbasis)

    ax = plt.gca()
    def expf(x, pos):
        return '{:0.2f}'.format(np.exp(x))
    fexpf = mpl.ticker.FuncFormatter(expf)
    ax.xaxis.set_major_formatter(fexpf)
    plt.plot(dbasis, divider_count, linestyle='None', marker='o')
    plt.xlabel('d')
    plt.ylabel(r'log(mean(\tau))')
    plt.tight_layout()

    return dbasis, divider_count
        
def finite_grad(basis, a):
    x = np.gradient(a, basis)
    return np.gradient(x, basis)

import scipy.ndimage
def convolve_grad(basis, a):
    dxdx = np.diff(basis)
    #ddf = np.diff(f, 2) / dxdx
    #ccf = np.convolve(f, [1, -2, 1]) / dxdx
    sigma= basis[1]-basis[0]
    ggf = scipy.ndimage.gaussian_filter1d(a[:-1], 
            sigma=1.0, order=2, mode='reflect') / dxdx
    return ggf

@command.defaultsave(pdir=command.pdir)
def divider_dimension(ft):
    #dbasis = np.linspace(0.03, 1.20, 41, True
    dbasis = np.logspace(-3, 1, 41, True)
    dbasis, divider_count = _divider_dimension(ft, dbasis)
    cut = slice(6,16)
    #norm_dim = np.mean(np.gradient(divider_count[cut], dbasis[cut]))
    norm_dim, intercept = np.polyfit(dbasis[cut], divider_count[cut], 1)
    plotutils.abline(plt.gca(), norm_dim, intercept)
    plt.title(r'{} D = {:.04f}'.format(ft.source, norm_dim))
    print('norm divider dimension', norm_dim) 

    dl = _analysefj.disp_length(ft)
    yrpix = 0.042 
    xguide = (yrpix, np.min(dl)/10.)
    for x in xguide:
        plt.axvline(np.log(x), alpha=0.6)
    
    print('end to end distance', np.min(dl), np.max(dl))


if __name__=='__main__':
    """
    """
    ff, thisargs = command.process(locals())
    ff(*thisargs)
