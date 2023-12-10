
import sys, os
import numpy as np
import copy
import collections
import pickle
from tqdm import tqdm


import datapaths
import readmat
import matdef

"""
Loading and saving FJ data for fast retrieval
"""

###################################################################################
# Load FanJin data

def init(first=True, alldirs=False, tracktype=matdef.DENOISED):
    initdirs = ddirs if alldirs else firstdir
    Ftracking = readmat.ExpTracks.usefj(initdirs, first, tracktype=tracktype)
    print('loading total of {} tracks'.format(Ftracking.ntracks))
    return Ftracking

# 5 days of experiments. 10 .mat files each with n cells each a time series

import pili
dataset = os.path.join(pili.root, '../fanjin/working_copy/')
firstdir = [os.path.join(dataset, './0019--2012-09-28/')]
ddirs_loc = [
        '0019--2012-09-28/', 
        '0020--2012-10-07/',
        '0028--2012-10-26/',
        '0031--2012-11-02/',
        '0039--2012-11-03/']
ddirs = [os.path.join(dataset, dn) for dn in ddirs_loc]

# the slow way
def _get_all_data():
    return _fj.init(False, True)

##################################################################################
# 
# Save/Load Pickled data
# Edit. Still Slow
pdir = os.path.join(dataset, 'pickled/')
datapaths.force_mkdir(pdir)
pfile = os.path.join(pdir, 'allft.pkl') 
def pdump():
    allft = init(1, 0)
    with open(pfile, 'w') as pf:
        pickle.dump(allft, pf)

def pload():
    with open(pfile, 'r') as pf:
        return pickle.load(pf)
    
#################################################################################

d = {matdef.ORIGINAL : "originalnpy/", matdef.DENOISED : "trnpy/"}

def trackdump(idx, trs, dumpdir='trnpy/'):
    form = os.path.join(dataset, dumpdir, 'track_{:05d}.pkl')
    for i, tr in zip(idx, trs):
        f_name = form.format(i)
        print('saving track to ', f_name)
        with open(f_name, 'wb') as fb:
            pickle.dump(tr, fb)

def trackload(idx, dumpdir='trnpy/'):
    form = os.path.join(dataset, dumpdir, 'track_{:05d}.pkl')
    trs = []

    seq = tqdm(idx) if len(idx) > 1 else idx
    for i in seq:
        name = form.format(i)
        with open(name, 'rb') as fb:
            trs.append(pickle.load(fb))
    return trs

# convience functions
def trackload_original(idx):
    return trackload(idx, dumpdir='originalnpy/')

def lintrackdump(idx, trs):
    return trackdump(idx, trs, dumpdir='ltrnpy/')

def lintrackload(idx):
    return trackload(idx, dumpdir='ltrnpy/')

def npydump():
    """write FJ tracks as numpy objects"""
    allft = init(0, 1)
    datapaths.force_mkdir(npydir)
    allft.writenpy(npydir)

def _nameify(name):
    return name+'.npy' if not name.endswith('.npy') else name

def npyloadall(idx):
    # idx is an index array
    print("Loading ({}) tracks {}".format(len(idx), idx))
    return readmat.ExpTracks.load_all(npydir, idx)

def npyload(idx):
    # load a single track
    print("loading track {} from .npy file".format(idx))
    return readmat.ExpTracks.load_one(npydir, idx)

N = 3113 # number of tracks in dataset
npydir = os.path.join(dataset, 'npy/')

class Fjslicing(object):
    """Operations relating to extracting a subset of Fj tracks"""

    N = 3113

    def __init__(self):
        here = os.path.dirname(__file__)
        self.fjslicedir = os.path.join(here, 'fjslice/')

    # take the names of files without extension
    def intersect(self, a, b):
        af,bf  = os.path.join(self.fjslicedir, a+'.npy'), os.path.join(self.fjslicedir, b+'.npy')
        aslice, bslice = np.loadtxt(af).astype(int), np.loadtxt(bf).astype(int)
        out= os.path.join(self.fjslicedir, '_'.join([a,b]))  +'.npy'
        isect = np.intersect1d(aslice, bslice)
        print('set 1 has {:d}, set 2 has {:d}'.format(aslice.size, bslice.size))
        print('intersection has {:d}'.format(isect.size))
        print('Writing intersection to ', out)
        np.savetxt(out, isect)

    def load(self, name):
        if name=='all':
            return np.arange(0,self.N,1, dtype=int)
        return np.loadtxt(os.path.join(self.fjslicedir, _nameify(name))).astype(int)

    def load_data(self, name, n=None):
        idx = self.load(name)
        if n != None and n > Fjslicing.N:
            raise ValueError("Asked for {} tracks from dataset with {} tracks".format(n,Fjslicing.N))
        idx = idx[:n]
        return idx, npyloadall(idx)

    def load_original_trs(self, name, n=None):
        idx = np.atleast_1d(self.load(name))
        return idx[:n], trackload(idx[:n], dumpdir='originalnpy/')

    def load_trs(self, name, n=None):
        # idx, trackmanager = self.load_data(name, n)
        # return idx, trackmanager.get_tracklike()
        idx = np.atleast_1d(self.load(name))
        return idx[:n], trackload(idx[:n])

    def load_linearized_trs(self, name, n=None):
        # idx, trackmanager = self.load_data(name, n)
        # return idx, trackmanager.get_tracklike()
        idx = np.atleast_1d(self.load(name))
        return idx[:n], lintrackload(idx[:n])

    def save(self, name, data):
        path = os.path.join(self.fjslicedir, _nameify(name))
        assert(data.ndim == 1)
        print('saving data slice to {}'.format(path))
        np.savetxt(path, data)

# singleton global object
slicehelper = Fjslicing()

def load_subset_idx():
    return collections.OrderedDict([
        ("candidate", np.array([2924])),
        ("top", slicehelper.load("candidates_whitelist")),
        ("half", slicehelper.load("halflist")),
        ("median", slicehelper.load("median")),
        ("walking", slicehelper.load("pure_walking")),
    ])
    
def load_subsets(linearised=True, original=False):
    subset = {}
    for name, idx in load_subset_idx().items():
        if original:
            subset[name] = trackload_original(idx)
        else:
            if linearised:
                subset[name] = lintrackload(idx)
            else:
                subset[name] = trackload(idx)
        
    return subset

# load lvel distributions
def load_subset_speed():
    distrib = {}
    for name, ltrs in load_subsets().items():
        distrib[name] = np.concatenate([ltr.get_step_speed() for ltr in ltrs])
    return distrib

#################################################################################
# FanJin uses linear regression to linearize trajectory for displacements less than 2 pixels (0.12 \mu m)
# he claims a spatial resolution of 0.03 \mu m and velocity resolution of 0.3 \mu m/s

step_d = 0.12
def step_to_distance(tr, xkey='x', ykey='y', step_d=step_d):
    time = tr._track['time']
    x, y = tr._track[xkey], tr._track[ykey]
    step_idx = [0]
    step_displacement = []
    def compute_steps():
        i = 0
        i_f = 1
        while i_f < time.size:
            d = np.sqrt( (x[i]-x[i_f])**2 + (y[i]-y[i_f])**2 )
            if d > step_d:
                step_idx.append(i_f)
                step_displacement.append(d)
                i = i_f
                i_f = i+1
            else:
                i_f += 1
        # the last partial step should be thrown away
        return step_idx, step_displacement
    return compute_steps()
    # step_idx, step_displacement = compute_steps()

def lin(tr, step_idx, xkey='x', ykey='y'):
    x, y = tr._track[xkey], tr._track[ykey]
    step_n = np.diff(step_idx)
    # rewrite the track
    for i, start_i in enumerate(step_idx[:-1]):
        n = step_n[i]
        
        # retrieve track[step_start:step_end] inclusive
        x_arr = x[start_i:start_i+n+1] 
        y_arr = y[start_i:start_i+n+1] 

        # linear regression, keep the endpoints fixed
        line_x = np.linspace(x_arr[0], x_arr[-1], n, endpoint=False)
        line_y = np.linspace(y_arr[0], y_arr[-1], n, endpoint=False)
        # if n > 0:
        # else:
        #     line_x = x[start_i+1]
        #     line_y = y[start_i+1]
        tr[xkey][start_i:start_i+n] = line_x
        tr[ykey][start_i:start_i+n] = line_y

def linearize(tr, step_d=step_d):
    # assert(tr.slice == slice(None))
    # don't modify the original data
    tr = copy.deepcopy(tr)

    step_idx, _ = step_to_distance(tr, step_d=step_d)
    trail_step_idx, _ = step_to_distance(tr, 'trail_x', 'trail_y', step_d=step_d)
    lin(tr, step_idx)
    lin(tr, trail_step_idx, xkey='trail_x', ykey='trail_y')
    # throw away the remaining partial step
    end = max(step_idx[-1], trail_step_idx[-1])
    tr._track = tr._track[:end+1]
    # this may mean that that end of the leading or trailing track is not linearised
    tr.step_idx = step_idx
    tr.trail_step_idx = trail_step_idx
    return tr

def test_linearize():
    idx, trs = slicehelper.load_trs('default_crawling_list', 100)
    ltrs = [linearize(tr) for tr in trs]
    step_n = np.concatenate([np.diff(tr.step_idx) for tr in ltrs])
    plt.hist(step_n)
    plt.show()

def test_linearize_287():
    idx = 288
    print('testing linearization of track', idx)
    tr = trackload([idx])[0]
    print('track size', tr.size)
    linearize(tr)


# ----------------------------------------------------------------
# walking analysis
norm = np.linalg.norm
def pole_travel(tr):
    # negative score implies poles should be flipped 
    adx = tr.get_step_dx()
    bdx = tr.get_step_dx(trail=True)
    atravel = np.sum(norm(adx, axis=1))
    btravel = np.sum(norm(bdx, axis=1))
    return atravel, btravel

def pole_travel_score(tr):
    atravel, btravel = pole_travel(tr)
    denom = (btravel+atravel)
    # hiding a divide by zero possibility
    if denom == 0:
        return 0.
    score = (btravel-atravel)/denom
    return score
    
def body_corr(tr):
    # compute the correlation between body orientation and movement direction
    apole = tr.get_head()[:,:2]
    bpole = tr.get_trail()[:,:2]
    # compute body direction
    body = apole - bpole
    normalize = norm(body,axis=1)[:,np.newaxis]
    with np.errstate(divide='ignore', invalid='ignore') as errstate:
        body = body/normalize 
    body[np.isnan(body)] = 0. # set nan to zero
    # compute movement direction
    center = (apole + bpole)/2
    center_dx = center[1:] - center[:-1]
    normv = norm(center_dx,axis=1)[:,np.newaxis]
    movement = center_dx/normv
    # dot product
    product = np.sum(movement*body[1:],axis=1)
    corr = np.mean(product)
    return corr

def flip_score(tr):
    return body_corr(tr) + pole_travel_score(tr)

# Note: by applying the pole flipping each time we use the data
# I am trying to avoid having to make my own version of the dataset
def redefine_poles(trs):
    # apply appropriate flipping to a whole set of tracks
    # report the indices of the flipped tracks
    flipped_idx = []
    flipped_score = []
    for idx, tr in enumerate(trs):
        score = flip_score(tr)
        if score < 0:
            tr.flip_poles()
            flipped_idx.append(idx)
            flipped_score.append(score)
            
    a, b = len(flipped_idx), len(trs)
    print("flipped {}/{} tracks ({:.1f}%)".format(a, b, 100*float(a)/b))
    return flipped_idx, flipped_score


if __name__=='__main__':
    pass
    # npydump()
    
    def save_trs():
        debug=None
        idx, trackmanager = slicehelper.load_data('all', debug)
        trs = trackmanager.get_tracklike()
        trackdump(idx, trs)

    def save_linearised():
        debug = None
        idx, trs = slicehelper.load_trs('all', debug)
        print('linearising tracks')
        ltrs = [linearize(tr) for tr in tqdm(trs)]
        lintrackdump(idx, ltrs)

    # save_trs()
    save_linearised()

    # test_linearize_287()
    


