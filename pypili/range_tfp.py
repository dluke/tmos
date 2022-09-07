import sys, os
from os import path
from collections import OrderedDict

import random
import copy
from multiprocessing import Process

# we also want twanalyse routines to run on sets of folders with the same name
# and some intelligent ways of instantly comparing all data as subplots

import numpy as np

import parameters
import tfputils
import ctfp3d as tfp

import pili.database.sqldb as sqldb

############################################################################

class Batch():

    def __init__(self, jobs):
        self.jobs = jobs
        self.finished = []

    def set_job_as_finished(job):
        assert(job in self.jobs)
        self.finished.append(job)

    def check_finished(self):
        return [job in self.finished for job in self.jobs]


class Job():

    def __init__(self, idx, cellid, directory, parameters=None, exit_condition=None, batch=None, algorithm=None):
        # the id of the job
        self.idx = idx
        # the cell id enumerates jobs with the same parameters
        self.cellid = cellid
        self.directory = directory
        # parameters is a dictionary of key value pairs
        self.parameters = parameters
        # reference to a group of jobs which are duplicates of this one
        self.batch = batch
        # reference to the algorithm which created this job
        self.algorithm = algorithm


        if exit_condition is None:
            self.exit_condition = tfputils.no_exit_condition 
        else:
            self.exit_condition = exit_condition

        self.args = None
        self.started = False
        self.finished = False

    def is_alive(self):
        return p.is_alive()

    def set_as_finished():
        if self.batch:
            self.batch.set_job_as_finished(self)
        self.finished = True

    def get_batch_list():
        return self.batch.jobs if self.batch else [self]

    def is_batch_started(self):
        if self.batch:
            return any(job.started for job in self.batch.jobs)
        else:
            return self.started

    def is_batch_complete(self):
        if self.batch:
            return all(self.batch.check_finished())
        else: # If job is not part of batch then it is it's own batch
            return True 

    # use psetter object to set the C++ paramters
    def setup(self, args):
        self.args = copy.deepcopy(args)
        if self.parameters:
            for name, value in list(self.parameters.items()):
                self.args.pset(name, value)

        if not os.path.exists(self.directory):
            os.mkdir(self.directory)
        # obtain random seed
        if self.args.system.randominit:
            self.args.system.seed = random.randrange(sys.maxsize)

    def new_process(self):
        if self.args.system.main == 'production':
            return tfp.run(self.args, self.cellid, self.exit_condition)
        elif self.args.system.main == 'attachment':
            args = self.args
            p = Process(target=tfp.run_attachment, args=(args, None, args.system.n_attach_stop))
            return p

    def assign_process(self):
        assert(self.args is not None)
        #
        curdir = os.getcwd()
        os.chdir(self.directory)
        self.p = self.new_process()
        self.p.start()
        self.started = True
        os.chdir(curdir)


# naming rule for the simulation directories
# looks like parameter1_value1_parameter2_value2_....
def naming_rule(parameters, args):
    format_vals = [args.params[name].form.format(val) for name, val in list(parameters.items())]
    iters = [iter(list(parameters.keys())), iter(format_vals)]
    parts = itertools.chain.from_iterable(itertools.repeat(iters, len(parameters)))
    dirname = '_'.join([next(part) for part in parts])
    return dirname


# allparams should be a dictionary {parameter_name : [values]}
# Construct all the Job objects for a given parameter set and add them to the list
import itertools
def buildpool(args, allparams, exit_condition=None, overwrite=False):
    # create a number of identical jobs for each parameter combination
    pool = []

    # this order of loops should give jobs of different types first
    # before repititions of the same type
    cwd = os.getcwd()
    pnames = list(allparams.keys())
    idx = 0

    for cellid in range(args.system.repeats):
        for paramset in itertools.product(*list(allparams.values())):
            parameters = OrderedDict(list(zip(pnames, paramset)))
            dirname = os.path.join(cwd, naming_rule(parameters, args))
            if not os.path.exists(dirname) or overwrite:
                pool.append( Job(idx, cellid, dirname, parameters, exit_condition) )
                idx += 1
            else:
                print('Will not overwrite data in ', dirname)
                pass 


    # later we pop() from the pool, so it's reversed here
    return list(reversed(pool))




# We can use this Job system for all simulation that use the tfputils.simulate_cell mainloop
# But debugging runs need their own mainloop
def buildsinglepool(args, exit_condition=None, dirname='./'):
    pool = [Job(i, i, dirname, None, exit_condition) for i in range(args.system.repeats)]
    return list(reversed(pool))

# Pretty self-explanatory
# Take a list of Jobs and run them all using <maxcores> cores
def runpool(pool, args, maxcores=12):
    # make a copy incase args gets modified by the simulation setup
    args = copy.deepcopy(args)

    # total 
    total = len(pool)
    ncells= args.system.repeats
    oneround = total//ncells

    print('There are {} jobs over {} repeats'.format(total, ncells))
    # print('Thats {} parameter combinations'.format(oneround))
    print('...')

    import time
    start_t = time.time()

    using = 0
    active = {}
    need_submit = 1
    while True:

        if len(pool) == 0:
            if need_submit == 1:
                print('Finished submitting all the jobs.')
            need_submit = 0
        if len(pool) == 0 and len(active) == 0:
            print('All processes finished.')
            break
        if using <= maxcores and need_submit:
            job = pool.pop()
            job.setup(args)
            print('starting job No. {}/{}, round number is {}/{}'.format(
                    job.idx+1, total, (job.idx//oneround)+1, ncells))
            job.assign_process()
            active[job.idx] = job
            using += 1

        # now check all the active jobs to see if they finished
        for job in list(active.values()):
            # might be more useful to use exitcode() than is_alive()
            #job.p.exitcode 
            if job.p.is_alive() == False:
                # this job finished
                print('job No. {}/{} finished'.format(job.idx+1, total))
                job = active.pop(job.idx)
                # can all job resources be freed now?
                #job.terminate()
                using -= 1

        time.sleep(0.01)

    end_t = time.time()
    exec_t = end_t - start_t
    print('Exited normally after {} seconds'.format(exec_t))
    with open("timing.txt", 'w') as ft:
        ft.write("full simtime {:11.4f} cores_used {}\n".format(exec_t, maxcores))




#########
# for optimisation
def_var = ['tau', 'pilivar', 'free_tau_eovers']

def read_record(rfile, variables=None):
    record = []
    with open(rfile, 'r') as f:
        firstline = next(f)
        if (firstline[0] == '#'):
            variables = firstline[1:].split()
        else:
            lvals= list(map(float, firstline.split()))
            record.append( (lvals[:-1], lvals[-1]) )
            variables = def_var
        #
        for line in f:
            lvals= list(map(float, line.split()))
            record.append( (lvals[:-1], lvals[-1]) )

    assert(variables is not None)
    return variables, record


#################################################################################

# output utilities

falp = "parameters.meta"
def pretty_allparams(allparams):
    """
    Write out the parameter sets used in this parameter search
    """
    with open(falp, 'w') as fp:
        fp.write("Parameters used in this parameter search.\n")
        fp.write("#\n")
        for name, values in list(allparams.items()):
            fp.write('{: <25s} {}\n'.format(name, list(values)))

#################################################################################

# Database interface
def getdirs(path):
    return sorted( next(os.walk(path))[1] )

def index_of(path):
    ind = path.split('_')[-1]
    assert(ind.isdigit())
    return int(ind)

class Dynamicpool(object):
    """
    construct jobs on the fly. Record their configuration and later their results into the database.
    """

    rundir = "/data/dan/run/db"

    def __init__(self, args, test=True, maxthreads=12, notify=False):
        self.test = test
        self.maxthreads = maxthreads
        self.notify = notify
        self.args = args
        #
        self.parentdir = None
        self.pending = []
        self.active = {}
        self.finished = []
        #
        self.localidx = 0

    def connect(self):
        return sqldb.connect(self.test)

    def newparent(self):
        parentdirs = getdirs(self.rundir)
        last_index = index_of(parentdirs[-1])
        parent_form = 'parent_{:012d}'
        new_parent = parent_form.format(last_index+1)
        os.mkdir(new_parent)
        self.parentdir = new_parent

    def get_lastrowid(self):
        conn = self.connect()
        sql = "SELECT MAX(ID) FROM simulation"
        cur = conn.execute(sql)
        max_id = cur.fetchone()[0]
        if max_id is None:
            max_id = 0
        return max_id

    def newsim(self):
        assert(self.parentdir != None)
        form = '{}/simulation_{:012d}'
        simdir = form.format(self.parentdir, get_lastrowid()+1)
        if job.batch:
            for job in job.batch.jobs:
                job.directory = simdir
        return simdir

    def start_next_job(self):
        if not self.pending:
            raise RuntimeError("No job to start.")
        if self.parentdir is None:
            self.newparent()
        job = self.pending.pop(0)

        if not job.is_batch_started():
            newsimdir = self.newsim() 
            job.directory = newsimdir

        active[job.idx] = job
        self.insert_sim_config(job)
        job.setup(self.args)
        job.assign_process()

    def insert_sim_config(job):
        sqldb.insert_simulation(job)

    def add_pending(self, job):
        job.idx = self.localidx
        self.localidx += 1
        self.pending.append(job)

    def update():
        for job in list(active.values()):
            if job.p.is_alive() == False:
                if self.notify:
                    print(('job {} at directory {} finished.'.format(job.idx, job.directory)))
                # print notification
                job = active.pop(job.idx)
                finished.append(job)
                # 
                self._on_finished(job) 
                #  
                if job.is_batch_complete():
                    self._on_batch_finished(job.get_batch_list())

                if len(self.active) < self.maxthreads and self.pending:
                    self.start_next_job()

    def _on_finished(self, job):
        pass

    def _on_batch_finished(self, jobs):
        """read in the data and compute metrics, etc. """
        pass


# null classes
class Metric(object):
    pass

class NoMetric(Metric):

    def __init__(self):
        self.name = NoMetric.__name__

class Algorithm(object):
    pass

class NoAlgorithm(Algorithm):
    def __init__(self):
        self.version = 0.0
        self.hexsha = None
        self.metric = NoMetric()


def test_dyanmic_pool():
    pool = Dynamicpool(None, test=True)
    rowid = pool.get_lastrowid()


if __name__=='__main__':

    test_dyanmic_pool()

