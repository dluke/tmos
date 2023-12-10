
# anneal the mapping

import os, sys
from copy import copy, deepcopy
import random
import pickle
join = os.path.join 
import collections
import numpy as np
import scipy
import scipy.optimize

pi = np.pi
norm = np.linalg.norm

from pili import support
import pili
import mdl

# globals

_debug = False
local_f =  sys.stdout

# ~

# -----------------------------------------------------------------
# wavelet

sigma = 0.03
scikit_config = {"wavelet":'db1', 'method':'VisuShrink', "mode":'soft', "rescale_sigma":False}

def vary_sigma(x, y, sigma, config=scikit_config):
    x_denoise = denoise_wavelet(x, sigma=sigma, **config)
    y_denoise = denoise_wavelet(y, sigma=sigma, **config)
    denoised = np.stack([x_denoise, y_denoise])
    return denoised

def contract(denoised, sigma):
    threshold = 1e-2 * sigma
    N = denoised[0].size
    xyt = denoised.T
    # for some reason if we choose large sigma, the denoised trajectory as repeated values
    # contract repeated values
    diff = xyt[1:] - xyt[:-1]
    eq = norm(diff, axis=1) < threshold
    seg_idx = np.insert(np.argwhere(eq == False).ravel()+1, 0, 0)
    return xyt[seg_idx].T

def model_from_denoised(denoised, sigma=0.04, angle_threshold=15):
    contracted = contract(denoised, sigma)

    # check the distances between nodes
    dt = np.zeros(contracted.shape[1])
    wavemodel = mdl.LPtrack(dt, contracted[0], contracted[1])
    wavemodel = mdl.recursive_coarsen(wavemodel, 2*sigma)

    def contract_by_angle(wavemodel, angle_threshold=5):
        # threshold in degrees
        theta = wavemodel.get_angle() * 180/np.pi
        list_theta = theta.tolist()[::-1]

        keepidx = [0]
        i = 1
        _curr = 0
        while list_theta:
            _curr += list_theta.pop()
            if abs(_curr) > angle_threshold:
                _curr = 0
                keepidx.append(i)
            i += 1
        keepidx.append(len(wavemodel)-1)
        n= len(keepidx)
        dt = np.zeros(n)
        model = mdl.LPtrack(dt, wavemodel.x[keepidx], wavemodel.y[keepidx])
        return model

    bmodel = contract_by_angle(wavemodel, angle_threshold)
    return bmodel


# -----------------------------------------------------------------
# 




