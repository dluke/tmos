

"""
sobol analysis
"""

from SALib.sample import saltelli
from SALib.analyze import sobol
import numpy as np
import pandas as pd
import json
import os, sys
join = os.path.join
import copy
import subprocess
import multiprocessing
import time
import resource
import collections

# tmos
import stats
import range_tfp
import readtrack
import command
import twutils
import rtw
import twanalyse
import _fj
from parameters import ParameterList

# ----------------------------------------------------------------
# parameter data and limits 

paramdata = [
    ['k_ext_off', 1.6],
    ['k_ret_off', 9.1],
    ['k_ext_on', 2.4],
    ['k_ret_on',  0.4],
    ['dwell_time', 1.0],
    ['pilivar', None],
    ['anchor_angle_smoothing_fraction', 1.0],
    ['k_spawn', 1.0],
    ['k_resample', 1.0]
]
 
_names = [p[0] for p in paramdata]

problem = {
    'num_vars': 9,
    'names': _names,
    'bounds': [
        [1.0,2.0],
        [5.0,15.0],
        [1.0,10.0],
        [0.0,1.0],
        #
        [0.5,3.0],
        [1.0,20.0],
        [0.125,1.0],
        [0.5,5.0],
        [1.0,10.0]
    ]
}

# ----------------------------------------------------------------
# utilities

def write_sampled(values):
    sample_data_name = 'problem_samples.npy'
    np.savetxt(sample_data_name, values)

def write_problem(problem):
    with open('problem.json', 'w') as f:
        json.dump(problem, f, indent=1)
    
def read_sampled(targetdir='.'):
    sample_data_name = join(targetdir, 'problem_samples.npy')
    values = np.loadtxt(sample_data_name)
    return values

def read_problem(targetdir='.'):
    with open(join(targetdir,'problem.json'), 'r') as f:
        problem = json.load(f)
    return problem

# generator constructor function
def new_uuid():
    N = 8
    chk = list("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789")
    check = set()
    while True:
        while True:
            new = ''.join(np.random.choice(chk, N))
            if new not in check:
                check.add(new)
                break
        yield new

# ----------------------------------------------------------------
#

prefix = '_u_'

class SimulationPool(object):
    """
    handle tmos setup
    """

    def __init__(self, root_sim_dir='.'):
        self.root_sim_dir = root_sim_dir
        # problem as defined for SALib
        self.problem = None
        # tmos configuration parameters
        self.args = None
        # take control of args.system.repeats parameter 
        self.nrepeats = 1
        # max cores for running repeats in parallel 
        self.maxcores = 1
        # 
        self.uuid = new_uuid()
        #
        self.sampled_params = None
        # 
        self.step_target = None
        # 
        self.max_repeats = 20


    def is_ready():
        if (args == None):
            raise RuntimeError('Pypili configuration is not set')
        return True

    def make_jobs(self):
        jobs = []
        lookup = [list(), dict()]
        for params in self.sampled_params:
            uid = prefix + next(self.uuid)
            jobs.append([uid, params]) 
            lookup[0].append(uid)
            lookup[1][uid] = params.tolist()
        with open('lookup.json', 'w') as f:
            print('writing lookup table to ', 'lookup.json')
            json.dump(lookup, f, indent=1)
        return jobs

    def _run(self, uid, x):

        args = copy.deepcopy(self.args)

        param_names = self.problem['names']
        for param, value in zip(param_names, x):
            args.pset(param, value)
        args.system.repeats = self.nrepeats

        new_dir = os.path.join(self.root_sim_dir, uid)
        if self.step_target:
            total_steps = 0
            for idx in range(self.max_repeats):
                job  = range_tfp.Job(idx, idx, new_dir) 
                job.setup(args)
                job.assign_process()
                job.p.join() # wait
                tr = readtrack.read_last_track(new_dir)
                _n = _fj.linearize(tr).get_nsteps()
                total_steps += _n
                if total_steps > self.step_target:
                    break 
        else:
            singlepool = range_tfp.buildsinglepool(args, exit_condition=None, dirname=new_dir)
            # waits until all jobs are finished
            range_tfp.runpool(singlepool, args, maxcores=1)

        # compute summary in parallel as well
        with command.chdir(new_dir):
            trs = readtrack.trackset()
            twanalyse._summary(trs)
            # help speed up memory usage tracking later
            os.system('du -s . > du.txt')

    def parallel_run(self, jobs):
        parallel_run(self._run, jobs, maxcores=self.maxcores)
        print('ran all simulations')

def read_lookup(targetdir):
    lfile = join(targetdir,'lookup.json')
    with open(lfile, 'r') as f:
        lookup = json.load(f)
    return lookup

# -------------------------------

def compute_stats():
    lookup = read_lookup()
    l_order, l_uid = lookup
    for uid in l_order:
        with command.chdir(uid):
            trs = readtrack.trackset()
            sd = twanalyse._summary(trs)

def parallel_compute_stats(maxcores=64):
    if not os.path.exists('lookup.json'):
        print("ERROR: lookup.json does not exist")
        sys.exit()
    with open("lookup.json", 'r') as f:
        lookup = json.load(f)
    def runner(i, udir):
        with command.chdir(udir):
            print("job {} for directory {}".format(i, udir))
            trs = readtrack.trackset()
            sd = twanalyse._summary(trs)
    jobs = list(enumerate(lookup[0]))
    parallel_run(runner, jobs, maxcores=maxcores)

def parallel_compute(analyze, maxcores=64):
    if not os.path.exists('lookup.json'):
        print("ERROR: lookup.json does not exist")
        sys.exit()
    with open("lookup.json", 'r') as f:
        lookup = json.load(f)
    def runner(i, udir):
        with command.chdir(udir):
            print("job {} for directory {}".format(i, udir))
            analyze()
    jobs = list(enumerate(lookup[0]))
    parallel_run(runner, jobs, maxcores=maxcores)

def parallel_compute_accepted(analyze, maxcores=64):
    def runner(i, udir):
        with command.chdir(udir):
            print("job {} for directory {}".format(i, udir))
            analyze()
    jobs = list(enumerate(np.loadtxt("standard_accepted_uid.npy", dtype=object).tolist()))
    parallel_run(runner, jobs, maxcores=maxcores)



# -------------------------------

def collect(objectives, targetdir='.', alldata=False, tdct={}, local_stats_path=None):
    l_order, l_uid = read_lookup(targetdir)
    sampled_values = read_sampled(targetdir)
    if local_stats_path is None:
        local_stats_path = stats.localjson 

    getters = []
    Y = {}
    for name in objectives:
        getters.append( rtw.make_get(name) )
        Y[name] = np.empty(sampled_values.shape[0], dtype=tdct.get(name, float))
    lduid = collections.OrderedDict()
    for i, uid in enumerate(l_order):
        with command.chdir(join(targetdir, uid)):
            ld = stats.load(json_path=local_stats_path)
            nodata = not bool(len(ld))
            if nodata:
                # print("warning: no data at ", uid)
                ld["failed"] = True
                ld["failed_condition"] = "no data"
            lduid[uid] = ld
            for k, name in enumerate(objectives):
                try:
                    value = getters[k](ld)
                except AttributeError:
                    value = np.nan
                Y[name][i] = value
    if alldata:
        return Y, lduid
    else:
        return Y

# turns out often easier to just collect all the summary data 
# and the pull out what we want later
def collect_lduid(targetdir='.', local_stats_path='local.json'):
    _, lduid = collect([], targetdir, alldata=True, local_stats_path=local_stats_path)
    return lduid

# -------------------------------
# utilities

def parallel_run(runner, jobs, maxcores=64):
    running = []

    def check_finished(running):
        still_running = []
        for p in running:
            if p.is_alive():
                still_running.append(p)
        running = still_running
        return running
    
    while(running or jobs):
        while(len(running) >= maxcores):
            running = check_finished(running)
            time.sleep(0.01)
        if jobs:
            job = jobs.pop(0)
            print('starting job with params {}'.format(job))
            p = multiprocessing.Process(target=runner, args=job)
            p.start()
            running.append(p)
            if not jobs:
                print('waiting for the remaining {} jobs to finish...'.format(len(running)))
        else:
            running = check_finished(running)
            time.sleep(0.01)


# -------------------------------
# analysis routines 
# call from notebooks
import SALib.analyze as sal
import pandas as pd

# convenience function to check for nan data in objectives
def check_missing(lookup, Y):
    missing_table = []
    has_nan = False
    for name in Y.keys():
        nan_idx = np.isnan(Y[name])
        nans = np.argwhere(nan_idx).ravel()
        # nan_uid = [uid for i, uid in enumerate(lookup[0]) if nan_idx[i] == True]
        missing_table.append([name, np.count_nonzero(nan_idx), nans])
        if np.any(nan_idx):
            has_nan = True
            # find first 'failed' analysis
            # nanid = np.argwhere(nan_idx)[0][0]
            # nanuid = lookup[0][nanid]
    return missing_table if has_nan else None
 
def collect_obs(lookup, lduid, subnames, scores):
    Yf = {}
    for i, subset in enumerate(subnames):
        getter = rtw.make_get(scores[i])
        Yf[subset] =  np.array([getter(lduid[uid]) for uid in lookup[0]])
    return Yf

def filter_missing(problem, lookup, nan_idx, calc_second_order=False):
    # nan_idx is the index of the simulation to remove from lookup
    # we also need to remove simulations in those affected blocks
    D = problem['num_vars'] 
    step = 2 * D + 2 if calc_second_order else D + 2
    N = len(lookup[0])//step
    block = np.zeros(N, dtype=bool)
    for idx in nan_idx:
        block_idx = idx // step
        block[block_idx] = 1
    mask = np.repeat(block, step)
    return ~mask
    # _uid_list = [lookup[0][i] for i in range(mask.size) if mask[i] == False]
    # return _uid_list
    # return lookup

class MCData(object):

    def __init__(self, problem, lookup, lduid):
        self.problem = problem
        self.lookup = lookup
        self.lduid = lduid
        self.valid = lookup[0]

    def get(self, objective):
        getter = twutils.make_get(objective)
        return np.array([getter(self.lduid[uid]) for uid in self.valid])

    def paramsdf(self, objectives):
        _lookup = self.lookup
        _parlist = list(zip(self.problem["names"],  zip(*[_lookup[1][_u] for _u in _lookup[0]])))
        _col = {}
        _col["uid"] = [_u  for _u in _lookup[0]]
        _col.update({k:v for k, v in _parlist})
        _data =  collect_obs(_lookup, self.lduid, objectives, objectives)
        _col.update({name:_data[name] for name in objectives})
        return pd.DataFrame(_col)

class SobolData(MCData):

    def __init__(self, problem, lookup, lduid, second_order):
        super().__init__(problem, lookup, lduid)
        self.second_order = second_order
        self.is_sobol = True

    def compute(self, Y):
        return compute_sobol(self.problem, Y, second_order=self.second_order)

    def format(self, S):
        return format_sobol(self.problem, S)


def compute_sobol(problem, Yf, second_order=True):
    Sf = {}
    for name in Yf.keys():
        Sf[name] = sal.sobol.analyze(problem, Yf[name], 
            calc_second_order=second_order)
    return Sf

def format_sobol(problem, Sf):
    _icol = {"parameter" : problem["names"]}
    table1 = _icol.copy()
    table2 = _icol.copy()
    table1.update({subset : S['S1'] for subset, S in Sf.items()})
    table2.update({subset : S['ST'] for subset, S in Sf.items()})
    dftable1 = pd.DataFrame(table1)
    dftableT = pd.DataFrame(table2)
    return dftable1, dftableT

def sortscore(problem, lookup, Yf, scores):
    # need to sort each column seperately, can't do this in one dataframe
    # construct a dataframe with cols [i, simulation_index, score] for each subset
    paramlist = problem["names"]
    sortdf = {}
    for subset, data in Yf.items():
        sortidx = np.argsort(data)
        udir = [lookup[0][idx] for idx in sortidx]
        _cols = {"index": sortidx, "dir": udir, "score": data[sortidx]}   
        _parlist = zip(problem["names"],  zip(*[lookup[1][_u] for _u in udir]))
        _cols.update({k:v for k, v in _parlist})
        _df = pd.DataFrame(_cols)
        sortdf[subset] = _df
    return sortdf

# -------------------------------
# (MOVED TO twutils)
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

# -------------------------------
# (MOVED TO abcimplement)
def mcsampling(problem, N):
    # sample a flat prior in the space defined by problem
    d = problem["num_vars"]
    rv = np.random.rand(N,d)
    for i, bounds in enumerate(problem["bounds"]):
        xn, xm = bounds
        rv[:,i] = rv[:,i] * (xm-xn) + xn
    return rv



if __name__=='__main__':
    pass



