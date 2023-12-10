
# module for analysing pwl trajectories


import os
join = lambda *x: os.path.abspath(os.path.join(*x))
import json
from glob import glob
import numpy as np
norm = np.linalg.norm
import scipy.stats
from scipy.optimize import curve_fit

import pili
from pili import support
import _fj
import mdl
import annealing
import pwlpartition

from skimage.restoration import estimate_sigma, denoise_wavelet

root = join(pili.root, '../sparseml')

# ----------------------------------------------------------------
# utilities

def angle(v1,v2):
    v1_u, v2_u = v1/norm(v1,axis=1)[:,np.newaxis], v2/norm(v2,axis=1)[:,np.newaxis]
    return np.sign(np.cross(v1_u, v2_u)) * np.arccos(np.clip(np.sum(v1_u * v2_u, axis=1), -1.0, 1.0))

# ----------------------------------------------------------------

def load_pwl(path):
    solver = annealing.Solver.load_state(join(path, "solver"))
    config_file = join(path, "config.json")
    with open(config_file, 'r') as f:
        config = json.load(f)
    return config, solver

def load_pwl_model(path):
    solver = annealing.Solver.load_state(join(path, "solver"))
    return solver.anneal.get_current_model()

# load the original data as LPtrack
def load_lptr(trackidx):
    track = _fj.trackload_original([trackidx])[0]
    data = mdl.get_lptrack(track)
    return data

def load_smoothed(trackidx):
    track = _fj.trackload([trackidx])[0]
    return track

def load_sim_pwl(target, refname='candidate'):
    # ! simulation dt is normalised, i.e. time interval is 1
    proc_dir = join(target, 'proc', refname)
    def _load(path):
        data = np.loadtxt(path)
        dt, x, y = data.T
        return mdl.LPtrack(dt, x, y)
    leading = [_load(path) for path in sorted(glob(join(proc_dir, 'lead_pwl_*.npy')))]
    trailing = [_load(path) for path in sorted(glob(join(proc_dir, 'trail_pwl_*.npy')))]
    return leading, trailing

# ----------------------------------------------------------------

whitelist = list(_fj.slicehelper.load("candidates_whitelist"))
whitelist.append(2924)
def is_whitelist(trackidx):
    return trackidx in whitelist


def track_summary(trackidx, wavetrack):
    # wavetrack is the original Fanjin lab wavelet smoothed data
    local = {}
    local['speed'] = np.mean(norm(wavetrack.get_head_v(), axis=1))
    local['whitelist'] = is_whitelist(trackidx)
    return local

def solver_summary(solver, track=None):
    local = {}
    # local['iteration'] = solver.count_iter
    local['dl'] = solver.get_description_length()
    local['loss'] = solver.partition.get_total_loss()
    data = solver.partition.data
    local['estimate_sigma'] = pwlpartition.estimate_error(data.x, data.y)
    local.update(empirical_summary(solver, track))
    # local.update(summary(solver.partition.model))
    return local

def summary(model):
    local = {}
    N = len(model) - 1 
    local['pwlm'] = len(model)
    local['N'] = N
    local['duration'] = model.get_duration()
    local['step_rate'] = N/model.get_duration()

    steps = model.get_step_length()
    local['contour_length'] = np.sum(steps)
    local['median_step_length'] = np.median(steps)
    local['mean_step_length'] = np.mean(steps)

    # * distributions *

    # model lengths
    local["lengths"] = steps
    # compute angles between adjacent segments
    xy = model.get_n2()
    u = xy[1:] - xy[:-1]
    alpha = angle(u[1:], u[:-1])
    local["angles"] = alpha

    return local

def true_summary(model):
    #! TODO: true models should be generated with an accurate dt
    #! this was lazy of me
    if model.dt is None:
        model.reset_dt()
    return summary(model)

def sim_summary(target, refname='candidate'):
    # simulation may generate one or many trajectories
    leading, trailing = load_sim_pwl(target, refname)
    if len(leading) == 0 or len(trailing) == 0:
        return None
    _DT = 0.1
    lookup_reduce = {
        'N' : np.sum,
        'duration' : lambda x: _DT * np.sum(x),
        'contour_length' : np.sum,
        'lengths' : np.concatenate,
        'angles' : np.concatenate
    }
    def _reduce(datalist):
        return {k : lookup_reduce[k]([local[k] for local in datalist]) for k in lookup_reduce}

    def _extra(local, simpwl):

        # !pwl model has very rare zero-length vectors for some reason
        dxdata = [model.get_vectors() for model in simpwl]
        for i, dx in enumerate(dxdata):
            dx = dx[~(norm(dx,axis=1)==0.)] 
            dxdata[i] = dx
        dxu = np.concatenate([dx[:-1] for dx in dxdata])
        dxv = np.concatenate([dx[1:] for dx in dxdata])
        alpha = angle(dxu, dxv)
        local["angle"] = {}
        if len(alpha) > 3:
            coef, pvalue  = scipy.stats.pearsonr(alpha[:-1], alpha[1:])
            local["angle"]["corrcoef"] = (coef, pvalue)
        else:
            local["angle"]["corrcoef"] = (np.nan, np.nan)
        
        return local
    
    lsummary = _extra(_reduce([summary(model) for model in leading]), leading)
    tsummary = _extra(_reduce([summary(model) for model in trailing]), trailing)
    return lsummary, tsummary
    


# ----------------------------------------------------------------
# distribution analysis of individual trajectores


def fit_exp_function(x, B, beta):
    return  np.exp(-(x - B)/beta) 

def fit_exp_curve(right_data):
    # fit exponential decay to histogram
    entries, bins = np.histogram(right_data, bins='auto')
    ydata = entries
    bincenters = bins[:-1] + np.diff(bins)
    # sigma = np.empty_like(ydata)
    # zeroidx = ydata == 0
    # sigma[zeroidx] = np.inf
    # sigma[~zeroidx] = 1/np.sqrt(ydata)[~zeroidx]
    sigma = 1/np.sqrt(ydata)
    popt, pcov = curve_fit(fit_exp_function, xdata=bincenters, ydata=ydata, sigma=sigma)
    return popt

def normal_function(x, A, sigma):
    return  A * np.exp(-x**2/sigma**2) 

def fit_normal_curve(data):
    # fit exponential decay to histogram
    size = data.max()
    bins = np.linspace(-size, size, 41, endpoint=True)

    # 
    xdata = bins[:-1] + np.diff(bins)
    min_data = 0.1
    cut_index = np.searchsorted(xdata, min_data,  side='right')
    xdata = xdata[cut_index:]

    entries, bins = np.histogram(data, bins=bins)

    ydata = entries[cut_index:]
    sigma = 1/np.sqrt(ydata)
    popt, pcov = curve_fit(normal_function, xdata=xdata, ydata=ydata, sigma=sigma)
    # print(bincenters, ydata)
    return popt


def empirical_summary(solver, track=None):
    local = {}
    model = solver.partition.model
    local.update(summary(model))


    if track:
        # compute angle between segment and body axis
        #! should we do some processing on the body axis? compute it from smooth data?
        segmentv = model.get_vectors()
        body = track.get_body_projection()
        time = model.get_time()
        deviation = angle(body[time[:-1]], segmentv)
        local["deviation"] = deviation 

    # 
    solver.partition.get_curve_coord()
    curve_coord = solver.partition.get_curve_coord()
    local["dx"] = np.diff(curve_coord)

    dx = local["dx"]
    max_dx = dx.max()
    right_dx = dx[ np.logical_and(dx > 0, dx < max_dx) ]
    offset, decay = fit_exp_curve(right_dx)
    local["right_exp"] = offset, decay
    local["stall_fraction"] = 1.0 - right_dx.size/dx.size

    return local




# ----------------------------------------------------------------
# analysis comparing solved and true models
# this primarily makes sense for synthetic data


def diff_models(m1, m2, meta):
    m1 = mdl.LPshape.from_lptr(m1)
    m2 = mdl.LPshape.from_lptr(m2)
    sampling_r = meta.get('sampling_r', meta["r"])

    # conting nodes
    d_M = m1.M - m2.M
    l1, l2 = m1.get_contour_length(), m2.get_contour_length()

    # total length
    d_L = l1 - l2

    # chi2 like metric which involves sampling the longer trajectory and then computing the shortest distance betweent them
    m_shorter, m_longer = (m1, m2) if l1 < l2 else (m2, m1)
    max_length = m_longer.get_contour_length()
    sample = np.arange(0, max_length+sampling_r/2, sampling_r)
    # print('sampling at r = ', sampling_r)
    m_pt = np.array([m_longer(v) for v in sample])
    # print('sampling size', m_pt.size)
    distance_matrix = np.empty((sample.size, m_shorter.M-1))

    short_data = m_shorter.get_n2()
    for i in range(m_shorter.M-1):
        a, b = short_data[i].reshape(1,2), short_data[i+1].reshape(1,2)
        s, d = support.line_coord_dist(m_pt, a, b)
        distance_matrix[:,i] = d

    min_distance = np.min(distance_matrix, axis=1)
    mean_separation = np.mean(min_distance)

    # compute a version of the maximum deviation 
    # get the curve coordinates which don't include the end segments
    inner_curve_coords = m_longer.get_cumulative_length()[[1, -2]]
    inner_left_index, inner_right_index = (np.searchsorted(sample, icc) for icc in inner_curve_coords)
    inner_min_distance = min_distance[inner_left_index : inner_right_index]
    max_sep_index = np.argmax(inner_min_distance) + inner_left_index
    max_separation = min_distance[max_sep_index]
    # print("maximum deviation at u({}) {} v = {}".format(sample[max_sep_index], m_pt[max_sep_index], max_separation))

    varname = ['d_M', 'd_L', 'mean_separation', 'max_separation']
    value = [d_M, d_L, mean_separation, max_separation]
    data = {k: v for k, v in zip(varname, value)}
    return data # (d_M, d_L, mean_separation, max_separation) 



# ----------------------------------------------------------------------
# my version of the wavelet transform

scikit_config = {"wavelet":'db1', 'method':'VisuShrink', "mode":'soft', "rescale_sigma":False}
def denoise_track(track, sigma, config=scikit_config):
    track =  track.copy()

    x = denoise_wavelet(track['x'], sigma=sigma, **config)
    y = denoise_wavelet(track['y'], sigma=sigma, **config)
    track['x'] = x
    track['y'] = y
    lmodel = pwlpartition.contract(np.stack([x,y]))
    step_idx = lmodel.get_time()
    step_idx[-1] -= 1 
    _fj.lin(track, step_idx)
    track.step_idx = step_idx

    x = denoise_wavelet(track['trail_x'], sigma=sigma, **config)
    y = denoise_wavelet(track['trail_y'], sigma=sigma, **config)
    track['trail_x'] = x
    track['trail_y'] = y
    tmodel = pwlpartition.contract(np.stack([x,y]))
    trail_step_idx = tmodel.get_time()
    trail_step_idx[-1] -= 1 
    _fj.lin(track, trail_step_idx, xkey='trail_x', ykey='trail_y')
    track.trail_step_idx = trail_step_idx

    return track

# ----------------------------------------------------------------------
# convenience

def load_candidate_sigma_r(target):
    solver = pwlpartition.Solver.load_state(join(target, "solver"))
    return solver.sigma, solver.r

def load_candidate_statistics(target):
    solver = pwlpartition.Solver.load_state(join(target, "solver"))
    # !tmp
    solver.partition.use_probability_loss = False
    solver.partition.inter_term = False
    # !~tmp

    candidate_summary = solver_summary(solver)
    return candidate_summary

def load_model_at(target):
    target_path = join(target, "solver.pkl")
    if not os.path.exists(target_path):
        print(target_path, "not found")
        return None
    solver = pwlpartition.Solver.load_state(join(target, "solver"))
    # !tmp
    solver.partition.use_probability_loss = False
    solver.partition.inter_term = False
    # !~tmp
    return solver.get_model()

