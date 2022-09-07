# utilites for twanalyse

# by adding purely analytical methods in separate file we free up twanalyse.py
# for executable methods which operate directly on the data/filesystem

import os
join = os.path.join
import subprocess
import json
import numpy as np
import scipy.optimize
import scipy.stats


import matplotlib.pyplot as plt

import command

# weighted mean and variance of array
def wmeanvar(arr, ws):
    # weighted mean
    mean = np.average(arr, weights=ws)
    # weighted variance
    variance = np.average((arr-mean)**2, weights=ws) 
    return mean, variance

def sync_directory(directory, force=False):
    cluster_home = "/cluster/rs_lab/dbarton"
    home = os.path.expanduser("~")
    relative = directory.lstrip(home)
    target, _tail = os.path.split(directory) 
    cluster_dir = join(cluster_home, relative)
    if not force and os.path.exists(join(target, "data/")): 
        return 
    command = ["rsync", "-av", "compute:"+cluster_dir, target]
    print("running command", command)
    output = subprocess.check_output(command)
    return output

# important method for reading summary statistics from json
def make_get(localvar):
    keys = localvar.split('.')
    def getter(localdata):
        _dct = localdata
        for k in keys:
            _dct = _dct.get(k, None)
            if _dct is None: 
                return np.nan
        return _dct
    return getter

def pretty_describe(arr):
    d = scipy.stats.describe(arr)
    s = ''
    for attr in ['nobs', 'minmax', 'mean', 'variance']:
        s += '{} = {}\n'.format(attr,getattr(d, attr))
    return s
# 
def trim(arr, alpha=0.05):
    return arr[np.logical_and(arr>np.quantile(arr,alpha), arr<np.quantile(arr,1.0-alpha))]

# asymmetric trim
def asytrim(arr, left=0.05, right=0.05):
    return arr[np.logical_and(arr > np.quantile(arr,left), arr < np.quantile(arr,1.0-right))]

def trim_tail(arr, alpha):
    print('trim using value ', np.nanquantile(arr,1.0-alpha))
    return arr[arr < np.nanquantile(arr,1.0-alpha)]

# write to stdout
def _pretty_dict(sd, cwidth=40):
    string = ''
    form = '{:<%ds} {}\n' % cwidth
    for key, val in sd.items():
        if isinstance(val, dict):
            string += _pretty_dict(val) + '\n'
        else:
            string += form.format(key, val)
    return string

def print_dict(sd, cwidth=40):
    print(json.dumps(sd, indent='\t'))

# fitting exponential functions

def exp_form(a,b,c):
        def exp(t):
            return a + b * np.exp(c*t)
        return exp

# just call plt.show() after this to see check the fit by eye
def fit_decay(time, data, plot=False):
    # a + b * exp(c*t)
    a,b,c = exp_est(time, data)
    exp_r = exp_form(a,b,c)
    # use this result as an initial guess for scipy optimize curve fit
    def sci_exp(t,a,b,c):
            return a + b * np.exp(c*t)
    popt, _ = scipy.optimize.curve_fit(sci_exp, time, data, p0=[a,b,c])
    
    exp_f = exp_form(*popt)

    if plot:
        @command.defaultsave()
        def show_fit_decay():
            fig = plt.figure(figsize=(12,9))
            ax = plt.gca()
            ax.plot(time, data, marker='D', linestyle='None')
            ax.set_xlabel(r'$\tau$')
            ax.set_ylabel(r'Corr(v,v)')
            ax.plot(time, [exp_r(t) for t in time], label='recursive')
            ax.plot(time, [exp_f(t) for t in time], label='scipy fit')
            ax.legend()
        show_fit_decay()
    return popt



"""
https://stackoverflow.com/questions/3938042/fitting-exponential-decay-with-no-initial-guessing

compute an exponential decay fit to two vectors of x and y data
result is in form y = a + b * exp(c*x).
ref. https://gist.github.com/johanvdw/443a820a7f4ffa7e9f8997481d7ca8b3
"""
def exp_est(x,y):
    n = np.size(x)
    # sort the data into ascending x order
    y = y[np.argsort(x)]
    x = x[np.argsort(x)]

    Sk = np.zeros(n)

    for n in range(1,n):
        Sk[n] = Sk[n-1] + (y[n] + y[n-1])*(x[n]-x[n-1])/2
    dx = x - x[0]
    dy = y - y[0]

    m1 = np.matrix([[np.sum(dx**2), np.sum(dx*Sk)],
                    [np.sum(dx*Sk), np.sum(Sk**2)]])
    m2 = np.matrix([np.sum(dx*dy), np.sum(dy*Sk)])

    [d, c] = (m1.I * m2.T).flat

    m3 = np.matrix([[n,                  np.sum(np.exp(  c*x))],
                    [np.sum(np.exp(c*x)),np.sum(np.exp(2*c*x))]])

    m4 = np.matrix([np.sum(y), np.sum(y*np.exp(c*x).T)])

    [a, b] = (m3.I * m4.T).flat

    return [a,b,c]


# Fairly fast for many datapoints, less fast for many costs, somewhat readable
# https://stackoverflow.com/questions/32791911/fast-calculation-of-pareto-front-in-python
def is_pareto_efficient_simple(costs):
    """
    Find the pareto-efficient points
    :param costs: An (n_points, n_costs) array
    :return: A (n_points, ) boolean array, indicating whether each point is Pareto efficient
    """
    is_efficient = np.ones(costs.shape[0], dtype = bool)
    for i, c in enumerate(costs):
        if is_efficient[i]:
        # Keep any point with a lower cost
            is_efficient[is_efficient] = np.any(costs[is_efficient]<c, axis=1)  
            is_efficient[i] = True  # And keep self
    return is_efficient
    