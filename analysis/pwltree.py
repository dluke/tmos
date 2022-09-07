

import sys, os
import time
import itertools
join = os.path.join

from copy import copy, deepcopy
import numpy as np
norm = np.linalg.norm

import scipy
import scipy.optimize

import pickle

from pili import support

import matplotlib.pyplot as plt
import matplotlib as mpl


import mdl
import pwlpartition


local_f = sys.stdout

def local_print(*args):
    print(*args, file=local_f)


# ----------------------------------------------------------------
# use the ideas from decision tree regression to construct a tree based approach to classifying our data into linear features

class Node(object):

    def __init__(self, idx, findex, tindex, parent=None):
        self.idx = idx
        self.findex = findex
        self.tindex = tindex
        self.parent = parent
        self.children = [] # could use first,next,prev,last (linked list tree)
        self.loss = None

        # the fitness of joining with the next node
        self.future_loss = None

    def length(self):
        return self.tindex - self.findex

    def __len__(self):
        return len(self.children)

    def  __lt__(self, other):
        return self.idx < other.idx

    def  __str__(self):
        return f"Node({self.idx})"

    def  __repr__(self):
        return str(self)


class Loss(object):
    # LeastSquares loss function

    def __init__(self, data, jac=False):
        self.data = data
        self.jac = jac
        self.sqresidual = np.array([0.0])

    def __call__(self, x):
        theta, cx = x

        if self.jac:
            px, py = self.data.T
            ux, uy = np.cos(theta), np.sin(theta)
            s = ux * (px - cx) + uy * py
            dx = cx + s*ux - px
            dy = s*uy - py
            self.sqresidual = dx*dx + dy*dy
            value = np.sum(self.sqresidual)
            _pdv_theta = np.sum( s * (dy  * ux -  dx * uy) )
            _pdv_cx = np.sum( dx * (1 - ux) - dy*uy )
            grad = 2 * np.array([_pdv_theta, _pdv_cx])
            return value, grad
        else:
            u = np.array([[np.cos(theta), np.sin(theta)]])
            C = np.array([[cx, 0]]) 
            self.sqresidual = support.line_distance(self.data, u, C)
            value = np.sum(self.sqresidual)
            return value

    def get_sqresidual(self):
        return self.sqresidual

class MaxLoss(Loss):

    def __call__(self, x):
        super().__call__(x)
        return np.sqrt(self.get_sqresidual().max())


def stop_at(r):
    def stop_condition(value):
        return value > r
    return stop_condition

class TreeSolver(object):

    def __init__(self, data, rng=None, params={}, optconf={}, 
            overlap=True, Loss=Loss, jac=True):
        self.rng = rng if rng!=None else np.random.RandomState(0)
        self.point_data = data.get_n2()
        self.N = len(self.point_data)
        self.jac = jac

        # config
        self.params = params
        self.optconf = optconf

        self.overlap = overlap
        self.Loss = Loss

        # 
        self.root = None
        # lookup the node by index
        # self.lookup = {}

        self._current_history_item = {}
        self.history = []


    def update_history(self, **kw):
        self._current_history_item.update(kw)

    def push_history(self):
        self.history.append(self._current_history_item)
        self._current_history_item = {}

    def get_history(self):
        return self.history

    def solve(self, stop_condition):

        while True:

            entry = self.pqueue.pop_entry()
            if entry is None:
                print("WARNING: priority queue is empty")
                break # warning queue is empty

            # local_print(entry)
            loss, (fnode, tnode) = entry
            self.update_history(loss=loss)

            if stop_condition(loss):
                break

            self.join(fnode, tnode)

            self.push_history()

    def build_priority(self):

        # for each adjacent pair
        nodes = self.root.children
        pairs = []
        for i in range(len(nodes)-1):
            fnode = nodes[i]
            tnode = nodes[i+1]
            loss = self.pair_loss(fnode, tnode)
            fnode.future_loss = loss
            pairs.append([loss, (fnode, tnode)])
        q = sorted(pairs, key=lambda t: t[0])
        self.pqueue = support.Modheapq.from_sorted(q)

    def priority_join(self):
        entry = self.pqueue.pop_entry()
        if entry is None:
            return None # empty queue
        loss, (fnode, tnode) = entry
        self.join(fnode, tnode)
        return loss

    def join(self, a, b):
        # * a, b are data ordered
        # local_print('join {} {}'.format(a,b), 'M = ', len(self.root))
        assert(a.parent == b.parent)
        parent = a.parent

        new = Node(next(self.index), a.findex, b.tindex, parent)
        a.parent = new
        b.parent = new
        new.children = [a, b]

        # remove these two nodes from parent and replace with the new node
        index = parent.children.index(a)
        M = len(self.root.children)

        if index > 0:
            self.pqueue.remove_task((parent.children[index-1], a))

        if index < M-2:
            self.pqueue.remove_task((b, parent.children[index+2]))

        parent.children.pop(index+1)
        M = M-1
        parent.children[index] = new

        if index > 0:
            fnode = parent.children[index-1]
            loss = self.pair_loss(fnode, new)
            self.pqueue.add_entry([loss, (fnode, new)])

        if index < len(self.root.children)-1:
            tnode = parent.children[index+1]
            loss = self.pair_loss(new, tnode)
            self.pqueue.add_entry([loss, (new, tnode)])
     
     
    def pair_loss(self, a, b):
        fi, ti = a.findex, b.tindex+1 if self.overlap else b.tindex
        data = self.point_data[fi:ti]
        loss = self.Loss(data, jac=self.jac)
        result = self.linear_fit(loss)
        
        # return result.fun
        value = np.sqrt(loss.get_sqresidual().max())
        return value

    def build_initial_tree(self, wavemodel):
        N = len(self.point_data)
        time = wavemodel.get_time()
        root = Node(None, 0, N)
        index = itertools.count()
        for i in range(wavemodel.M-1):
            node = Node(next(index), time[i], time[i+1], parent=root)
            root.children.append(node) 
            # self.lookup[node.idx] = node

        self.index = index
        self.root = root

    def compute_loss(self, node):
        if node.length() == 1:
            node.loss = 0.0
        else:
            data = self.point_data[node.findex:node.tindex]
            result = self.linear_fit(self.Loss(data, jac=self.jac))
            node.loss = result.fun


    def linear_fit(self, loss):
        data = loss.data
        if len(data) <= 2:
            return scipy.optimize.OptimizeResult({'fun': 0.0})
        v = data[-1] - data[0]
        guess_theta = np.arctan2(v[1],v[0])
        guess_cx = data[0][0]
        x0 = [guess_theta, guess_cx]
        
        res = scipy.optimize.minimize(loss, x0, jac=self.jac)
        return res
        
    def get_disconnected_model(self):
        lines = []
        for node in self.root.children:
            fi, ti = node.findex, node.tindex-1
            a = self.point_data[fi]
            b = self.point_data[ti]
            lines.append([a,b])
        return lines

    def get_model(self):
        # note we snap to the data here but so far we only solver for infinite linear features
        dt = [0]
        pt = [self.point_data[0]]
        for node in self.root.children:
            dt.append(node.length())
            ti = node.tindex if self.overlap else node.tindex-1
            ti = min(ti, self.N-1)
            pt.append(self.point_data[ti])
        x, y = np.stack(pt,axis=1)
        return mdl.LPtrack(dt, x, y)

# ---------------------------------------------------------------- 
# top level methods

def solve_lead_trail(track):
    dt = np.insert(np.diff(track['time']), 0, 0)
    lead_data = mdl.LPtrack(dt, track['x'], track['y'])
    trail_data = mdl.LPtrack(dt, track['trail_x'], track['trail_y'])
    #  ...


def solve_simdata(data, sigma, r):
    wavemodel = pwlpartition.sim_wavemodel(data, sigma)
    tsolver = TreeSolver(data, overlap=True)
    tsolver.build_initial_tree(wavemodel)
    tsolver.build_priority()
    tsolver.solve(stop_at(r))
    return tsolver

import readtrack
import multiprocessing

def pwl_proc(target, sigma, r, refname, parallel=False):
    """
    # target is the simulation directory
    # r is the threshold to use for pwl solver
    # refname is the unique name for this r-value, i.e. 'candidate'
    """

    # create a new directory for processed data
    proc_dir = join(target, 'proc', refname)
    if not os.path.exists(proc_dir):
        os.makedirs(proc_dir)
        
    # load the simulation tracks
    tracks = readtrack.trackset(ddir=join(target,'data/'))
    form1 = join(proc_dir, 'lead_pwl_{:05d}.npy')
    form2 = join(proc_dir, 'trail_pwl_{:05d}.npy')

    def solve(i, track):
        dt = np.insert(np.diff(track['time']), 0, 0)
        lead_data = mdl.LPtrack(dt, track['x'], track['y'])
        trail_data = mdl.LPtrack(dt, track['trail_x'], track['trail_y'])
        
        print(f'solving leading pole for track {i}')
        with support.PerfTimer() as t:
            solver = solve_simdata(lead_data, sigma, r)
        lead_model = solver.get_model()
        lead_model.save_array(form1.format(i))
        print(f'solved in {t.get_time()}s, M = {len(lead_model)}')

        print(f'solving trailing pole for track {i}')
        with support.PerfTimer() as t:
            solver = solve_simdata(trail_data, sigma, r)
        trail_model = solver.get_model()
        trail_model.save_array(form2.format(i))
        print(f'solved in {t.get_time()}s, M = {len(trail_model)}')

    if not parallel:
        for i, track in enumerate(tracks):
            solve(i, track)
    else:
        for i, track in  enumerate(tracks):
            p = multiprocessing.Process(target=solve, args=(i,track))
            p.start()

if __name__=='__main__':
    import pwlstats
    path =  join(pwlstats.root, "run/partition/candidate/_candidate_pwl/")
    sigma, r = pwlstats.load_candidate_sigma_r(path)
    local_pwl_proc = lambda : pwl_proc('./', sigma, r, 'candidate', parallel=True)
    local_pwl_proc()

