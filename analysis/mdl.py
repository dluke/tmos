
import pickle
import time
import itertools
import numpy as np
import scipy
norm = np.linalg.norm

import llist
from pili import support
import pandas as pd

import matplotlib.pyplot as plt


class Tracknode(object):
    def __init__(self, uidx, dt, x , y):
        self.uidx = uidx #  not needed?
        self.dt = dt
        self.x = x
        self.y = y
    
    def distance(self, other):
        return np.sqrt((self.x - other.x)**2 + (self.y - other.y)**2)

    def __hash__(self):
        return self.uidx

    def __eq__(self,other):
        return self.uidx == other.uidx
    
    def __repr__(self):
        return "Tracknode({},{},{},{})".format(self.uidx,self.dt,self.x,self.y)

class LPtrack(object):
    def __init__(self, dt, x, y, sidx=None):
        # linear piecewise trajectory
        self.M = len(x)
        self.dt =  dt 
        self.x =  x
        self.y =  y
        if sidx is not None:
            self.sidx = sidx
            self.max_index = np.max(sidx)
        else:
            self.max_index = self.M-1
            self.sidx = self._default_sidx()
        

    def reset_sidx(self):
        self.sidx = self._default_sidx()

    def reset_dt(self):
        self.dt = np.array(list(range(self.M)))

    def _default_sidx(self):
        sidx = np.array(list(range(self.M)))
        sidx[-1] = -1 # invalid value for index
        return sidx

    def _new_index(self):
        self.max_index += 1
        return self.max_index

    def expand_index(self, si, sf, local=1):
        return max(0, si-local), min(self.M, sf+local)

    def lengths_at(self, si, sf):
        # include the length of the segment to the endpoint
        _x = self.x[si:sf]
        _y = self.y[si:sf]
        L = np.sqrt( np.diff(_x)**2 + np.diff(_y)**2 )
        return L

    def get_vectors(self):
        disp = np.stack([self.x,self.y],axis=1)
        return disp[1:] - disp[:-1]

    def get_distance(self):
        disp = np.stack([self.x,self.y],axis=1)
        v = disp[1:] - disp[:-1]
        return norm(v, axis=1)

    def get_step_length(self):
        return self.get_distance()

    def get_break_coord(self):
        return np.insert(np.cumsum(self.get_step_length()), 0, 0)
    
    def get_ab_at(self, si, sf):
        # get model segment start and end positions
        _x = self.x[si:sf+1]
        _y = self.y[si:sf+1]
        _ab = np.stack([_x, _y]).T
        return _ab

    def get_point(self, index):
        return np.array([self.x[index], self.y[index]])

    def delete(self, si):
        self.x = np.delete(self.x, si)
        self.y = np.delete(self.y, si)
        _next = si+1 if si+1 < self.M else si-1
        self.dt[_next] += self.dt[si]
        self.dt = np.delete(self.dt, si)
        self.M -= 1

        self.sidx = np.delete(self.sidx, si)

    def create(self, si, x, y, dt1):
        assert(si < self.M-1) 
        sdt = self.dt[si+1]
        dt2 = sdt - dt1   #  split the count into two parts
        # print('dt1, dt2', dt1, dt2)
        self.x = np.insert(self.x, si+1, x)
        self.y = np.insert(self.y, si+1, y)
        self.dt[si+1] = dt2
        self.dt = np.insert(self.dt, si+1, dt1)
        self.M += 1

        self.sidx = np.insert(self.sidx, si+1, self._new_index())

    def move_data_at(self, index, dn):
        if dn > 0:
            self.move_data_left_at(index, dn)
        elif dn < 0:
            self.move_data_right_at(index, -dn)

    def move_data_right_at(self, index, dn):
        self.dt[index] -= dn
        self.dt[index+1] += dn

    def move_data_left_at(self, index, dn):
        self.dt[index+1] -= dn
        self.dt[index] += dn
         
    def insert(self, si, xydata):
        if len(xydata) == 0:
            return
        print('insert', si, xydata)
        newx, newy = xydata.T
        newdt = np.zeros(newx.size, dtype=int)
        # insert new points at index si
        # todo can use np.insert 
        self.x = np.concatenate([self.x[:si+1], newx, self.x[si+1:]])
        self.y = np.concatenate([self.y[:si+1], newy, self.y[si+1:]])
        self.dt = np.concatenate([self.dt[:si+1], newdt, self.dt[si+1:]])
        self.M = len(self.dt)
        

    def get_limits(self, buffer=(0,0)):
        buffx, buffy = buffer
        lx = [self.x.min()-buffx, self.x.max()+buffx]
        ly = [self.y.min()-buffy, self.y.max()+buffy]
        return (lx, ly)

    def get_bounding_size(self, buffer=(0,0)):
        lx, ly = self.get_limits(buffer)
        return np.array([lx[1]-lx[0], ly[1]-ly[0]])

    def __len__(self):
        return self.M

    def __eq__(self, other):
        if other is None:
            return False
        return all([(self.dt == other.dt).all(), (self.x == other.x).all(), 
            (self.y == other.y).all(), (self.sidx == other.sidx).all()])

    def __iter__(self):
        return iter(zip(self.dt, self.x, self.y))

    def cut(self, a, b):
        # a,b in time units
        time = np.cumsum(self.dt)
        idx1, idx2 = np.searchsorted(time,a), np.searchsorted(time,b,side='right')+1
        new = LPtrack(self.dt[idx1:idx2].copy(), self.x[idx1:idx2].copy(), self.y[idx1:idx2].copy(), self.sidx[idx1:idx2].copy())
        return new

    def cut_index(self, idx1, idx2):
        idx2 = idx2+1
        if self.dt is None:
            new = LPtrack(None, self.x[idx1:idx2], self.y[idx1:idx2], self.sidx[idx1:idx2])
        else:
            new = LPtrack(self.dt[idx1:idx2], self.x[idx1:idx2], self.y[idx1:idx2], self.sidx[idx1:idx2])
        return new


    def get_n2(self):
        # format x,y as (n,2) shape
        return np.stack([self.x, self.y]).T

    def get_time(self, DT=0.1):
        # return np.insert(np.cumsum(self.dt[:-1]), 0, 0)
        return np.cumsum(self.dt)

    def get_duration(self, DT=0.1):
        return np.sum(self.dt)

    def get_angle(self):
        xyt = self.get_n2()
        u = xyt[1:] - xyt[:-1]
        def angle(v1,v2):
            v1_u, v2_u = v1/norm(v1,axis=1)[:,np.newaxis], v2/norm(v2,axis=1)[:,np.newaxis]
            return np.sign(np.cross(v1_u, v2_u)) * np.arccos(np.clip(np.sum(v1_u * v2_u, axis=1), -1.0, 1.0))
        theta = angle(u[1:], u[:-1]) 
        return theta


    def get_llist(self):
        idx = range(len(self.dt))
        _seq = [Tracknode(_i, _dt, _x, _y) for _i, _dt, _x, _y in zip(idx, self.dt, self.x, self.y)]
        return llist.dllist(_seq)

    def to_dataframe(self):
        return pd.DataFrame({"dt": self.dt, "x": self.x, "y": self.y, "sidx": self.sidx})

    @classmethod
    def from_track(cls, track):
        dt = np.insert(np.diff(track['time']), 0, 0)
        return cls(dt, track['x'], track['y'])

    def save_array(self, path):
        data = np.column_stack([self.dt, self.x, self.y])
        np.savetxt(path, data)

    @classmethod
    def load_array(cls, path):
        dt, x, y = np.loadtxt(path).T
        return cls(dt, x, y)
    

    # def __str__(self):
    #     return "LPtrack with {} points".format(len(self.x))

# todo this inheritance relationship should be the other way around (?)

class LPshape(LPtrack):

    def __init__(self, x, y):
        super().__init__(None, x, y)
        self._lengths = None
        self._cumulative_length = None
        self.update_lengths()

    def update_lengths(self):
        self._lengths = self.get_distance()
        self._cumulative_length = np.insert(np.cumsum(self._lengths), 0, 0)

    def get_cumulative_length(self):
        return self._cumulative_length

    # def __str__(self):
    #     return "LPtrack with {} points".format(len(self.x))

# todo this inheritance relationship should be the other way around (?)

class LPshape(LPtrack):

    def __init__(self, x, y):
        super().__init__(None, x, y)
        self._lengths = None
        self._cumulative_length = None
        self.update_lengths()

    def update_lengths(self):
        self._lengths = self.get_distance()
        self._cumulative_length = np.insert(np.cumsum(self._lengths), 0, 0)

    def get_cumulative_length(self):
        return self._cumulative_length

    def get_contour_length(self):
        return self._cumulative_length[-1]

    def get_index(self, s):
        return np.searchsorted(self._cumulative_length, s, 'right') -1


    @classmethod
    def from_lptr(cls, lptr):
        new = cls(lptr.x, lptr.y)
        return new

    def __call__(self, v):
        # todo vectorise
        cumd = self._cumulative_length
        x, y = self.x, self.y
        max_v = cumd[-1]
        if v > max_v:
            return np.array([x[-1],y[-1]])
        s_i = np.searchsorted(cumd[1:], v, side='left')
        ex_v = v - cumd[s_i]
        # line equation
        _x, _y = x[s_i], y[s_i]
        _xn, _yn = x[s_i+1], y[s_i+1]
        m = ex_v/self._lengths[s_i] 
        return np.array([_x  +  m * (_xn - _x), _y  +  m * (_yn - _y)])


        

# lptrack friends

def get_lptrack(tr):
    timestep = 0.1
    dt = np.full(len(tr['x']), timestep)
    dt[0] = 0
    return LPtrack(dt, tr['x'], tr['y'])

def _regularise_extension(filename, ext):
    _name = filename
    if not filename.endswith(ext):
        _name = filename + ext
    return _name

def save_model(model, filename):
    _name = _regularise_extension(filename, '.pkl')
    print("writing to ", _name)
    with open(_name, 'wb') as f:
        pickle.dump(model, f)

def load_model(filename):
    _name = _regularise_extension(filename, '.pkl')
    print("reading from ", _name)
    with open(_name, 'rb') as f:
        return pickle.load(f)


class OperableTrack(object):
    def __init__(self):
        self.N = None
        self.lplist = None
        self.node_lookup = None
        self.dlist = None

    def update_distance(self, prev, node):
        prev_entry = self.dlist.entry_finder[prev.value.uidx]
        prev_node = self.node_lookup[prev_entry[1]]
        _distance = prev_node.value.distance(node.value)
        self.dlist.update_priority(prev_entry, _distance)

def link_to_lptr(lplist):
    dt = []
    x = []
    y = []
    for item in lplist:
        dt.append(item.dt)
        x.append(item.x)
        y.append(item.y)
    return LPtrack(np.array(dt), np.array(x), np.array(y))
    

def nodesum(op, node_a, node_b):
    lplist = op.lplist
    a = node_a.value
    b = node_b.value
    a.x = (a.x + b.x)/2
    a.y = (a.y + b.y)/2
    a.dt = a.dt + b.dt
    # what if node b is the end node? 
    lplist.remove(node_b)
    # update distance from previous node
    if node_a.prev:
        op.update_distance(node_a.prev, node_a)
    # remove distance associated with node b
    if node_a.next != None:
        op.dlist.remove_task(node_b.value.uidx)
        # update distance to next node
        _distance = a.distance(node_a.next.value)
        entry = [_distance, a.uidx]
        op.dlist.add_entry(entry)
    else:
        pass
        # print('do not remove task at', node_b.value.uidx)
    
def recursive_coarsen(_lptr, par_value, parameter="l"):
    assert(parameter in ['l', 'M'])

    # repeatedly connect points in the trajectory that are less than r apart
    op = OperableTrack()
    N = len(_lptr)
    op.lplist = _lptr.get_llist()
    op.N = N
    dist = _lptr.get_distance()
    noderange = list(range(len(op.lplist)))
    op.node_lookup = dict([(_i, op.lplist.nodeat(_i)) for _i in noderange])
    pairs = [[d, node] for d, node in zip(dist, noderange[:-1])]
    dtup = sorted(pairs, key=lambda t: t[0])
    op.dlist = support.Modheapq.from_sorted(dtup)

    #
    def l_condition(op, par_value):
        return op.dlist.pq[0][0] < par_value
    def M_condition(op, par_value):
        return len(op.lplist) > par_value
    select_condition = {
        'l' : l_condition,
        'M' : M_condition
    }
    condition = select_condition[parameter]

    #
    while condition(op, par_value):
        entry = op.dlist.pop_entry()
        node = op.node_lookup[entry[1]]
        nodesum(op, node, node.next)

    # print("new track with {}={} has {} nodes".format(parameter, par_value, len(op.lplist)))
    model = link_to_lptr(op.lplist)
    model.dt[1] += model.dt[0]
    model.dt[0] = 0
    return model


# TODO dep this for faster methods
# need a more robust compute_separation function
def robust_compute_separation(data_lptr, _x, intermediate=None):
    # _x contains the model points
    n = _x.size//2
    ab = _x.reshape((n,2))
    a = ab[:-1]
    b = ab[1:]
    point_tracker = [[] for _ in range(len(a))]
    # record_distances = [{} for _ in range(len(a))]
    distance_matrix = np.zeros((len(a),len(data_lptr)))
    distance = []
    pt_i = 0
    for _dt, _x, _y in data_lptr:
        p = np.array([_x, _y])
        dls = support.lineseg_dists(p, a, b)
        distance_matrix[:,pt_i] = dls
        segidx = np.argmin(dls)
        point_tracker[segidx].append(pt_i)
        dist = dls[segidx]
        distance.append(dist)
        pt_i += 1
    if intermediate != None:
        intermediate["point_tracker"] = point_tracker
        intermediate["distance_matrix"] = distance_matrix
    return np.array(distance)

# define the MDL cost function 
def description_length(data_lptr, coarse_lptr, r, intermediate=None):
    # data_lptr contains the original data
    # coarse_lptr is a piecewise linear model 
    if intermediate != None  and "seglist" in intermediate:
        seglist = intermediate["seglist"]
        separation = seglist.get_closest_segment()
    else:
        _x = coarse_lptr.get_n2()
        separation = robust_compute_separation(data_lptr, _x.flatten())
    noise = separation > r
    if intermediate != None:
        intermediate["is_outlier"] = noise
    return len(coarse_lptr) + np.sum(noise)

def new_intermediate():
    return {}

# --------------------------------------------------------------
# --------------------------------------------------------------
# spatio-temporal indexing

# neighbour list ideas
class SegmentList(object):
    
    # model comes with approximate times preloaded
    def __init__(self, data, model):
        self.model = model # LPtrack
        self.data = data # LPtrack
        N = len(self.data)
        M = len(self.model)
        #
        self.data_time = np.cumsum(data.dt)
        self.segment_time = np.cumsum(self.model.dt)
        self.point_data = data.get_n2()
        self.M = M
        self.N = N
        #
        self._reset()
        self.closest_segment = None


    def _reset(self):
        M, N = self.M, self.N
        # self.distance_matrix = np.full((M-1,N), np.inf)
        # can almost use scipy.sparse for this ... Almost ...
        self.distance_matrix = [{} for _ in range(N)]

    def _update_time(self):
        _dt = 0.1
        _time = np.zeros(self.M-1)
        for segment in self.closest_segment:
            _time[segment] += _dt
        self.segment_time = np.zeros(self.M)
        self.segment_time[1:] = np.cumsum(_time)

    @staticmethod
    def minargmin(distance_matrix):
        # needed because scipy.sparse not quite sufficient
        min_d = []
        min_k = []
        for dct in distance_matrix:
            d = np.inf
            k = 0
            for _k, _v in dct.items():
                if _v < d:
                    d = _v
                    k = _k
            min_d.append(d)
            min_k.append(k)
        return np.array(min_k), np.array(min_d)

        
    def update_distances(self, x):
        self.x = x
        self._reset()
        def _to_index(dt):
            return int(round(10*dt))
        def _get_adjacent_segments(): 
            segidx = np.arange(self.M-1, dtype=int)
            pre = np.roll(segidx, 1)
            aft = np.roll(segidx, -1)
            adjacent = np.stack([pre,segidx,aft],axis=1)
            return adjacent
        adjacent = _get_adjacent_segments()
        m_pt = x.reshape((self.M,2))
        current_segment = 0
        def _local_distance(adj_seg_idx):
            # compute the distances from all nearby points to this segment update distance matrix
            a = m_pt[current_segment].reshape((1,2))
            b = m_pt[current_segment+1].reshape((1,2))
            # select points by time-proximity to adjacent segments
            first, last = adj_seg_idx[0], adj_seg_idx[-1]+1
            ft, lt = self.segment_time[first], self.segment_time[last] 
            fidx, lidx = _to_index(ft), _to_index(lt)+1
            p = self.point_data[fidx:lidx]
            dist = support.lineseg_dists(p, a, b)
            # update distance matrix
            for pt_i, l_d in zip(range(fidx, lidx), dist):
                d = self.distance_matrix[pt_i].get(current_segment, np.inf)
                if l_d < d:
                    self.distance_matrix[pt_i][current_segment] = l_d
            return dist
        # 
        adj_seg_idx = adjacent[0][1:]
        _local_distance(adj_seg_idx)
        current_segment += 1
        for seg_i in range(1,self.M-2):
            adj_seg_idx = adjacent[current_segment]
            _local_distance(adj_seg_idx)
            current_segment += 1
        adj_seg_idx = adjacent[-1][:-1]
        _local_distance(adj_seg_idx)
        current_segment += 1
        # closest segment to each point
        # self.closest_segment = np.argmin(self.distance_matrix, axis=0)
        self.closest_segment, self.d =  self.minargmin(self.distance_matrix)
        self._update_time()
    
    def get_closest_segment(self):
        # 
        # indices = zip(self.closest_segment, range(self.N))
        # d = np.array([self.distance_matrix[a,b] for a,b in indices])
        return self.d

    def get_closest_data(self):
        # closest data points for each segment
        point_tracker = [[] for _ in range(self.M-1)]
        for pt_i, segment in enumerate(self.closest_segment):
            point_tracker[segment].append(pt_i)
        m_pt = self.x.reshape((self.M, 2))
        m_distance = []
        for i in range(self.M):
            pt_list = []
            if i > 0:
                pt_list.extend( point_tracker[i-1] )
            if i < self.M-1:
                pt_list.extend( point_tracker[i] )
            d = 0
            if len(pt_list) > 0:
                p = np.array([self.point_data[pt_i] for pt_i in pt_list])
                dist = norm(p - m_pt[i], axis=1)
                d = np.min(dist)
            m_distance.append(d)
        return np.array(m_distance)


# --------------------------------------------------------------
# --------------------------------------------------------------
# scipy optimisation

# construct optimisation function
def model_lsq(data_lptr, model_lptr, conf={}, intermediate=None):
    coarse_pt = model_lptr.get_n2()
    x0 = coarse_pt.flatten()
    quiet = conf.get("quiet", True)

    seglist = SegmentList(data_lptr, model_lptr)
    intermediate["seglist"] = seglist

    def end_terms(x, data, inter):
        point_tracker = inter["point_tracker"]
        distance_matrix = inter["distance_matrix"]
        extra_terms = []
        _M = x.size//2
        for segment_i in range(_M-1):
            close_points = point_tracker[segment_i]
            if not close_points: continue
            if segment_i > 0:
                pt_i = point_tracker[segment_i][0]
                d = distance_matrix[segment_i-1][pt_i]
                extra_terms.append(d)
            if segment_i < distance_matrix.shape[0]-1:
                pt_i = point_tracker[segment_i][-1]
                d = distance_matrix[segment_i+1][pt_i]
                extra_terms.append(d)
        return np.array(extra_terms)

    # be careful! This may not work for overlapping trajectories
    def constrain_model_to_data(x, data, inter):
        point_tracker = inter["point_tracker"]
        data_pt = data.get_n2()
        _M = x.size//2
        ab = x.reshape((_M,2))
        # slow!
        extra_terms = []
        for model_pt in ab:
            d = np.min(norm(data_pt - model_pt, axis=1))
            extra_terms.append(d)
        return np.array(extra_terms)

    add_end_terms = False # doesn't work
    add_model_terms = True 
    def loss(x, data, inter):
        # distance = robust_compute_separation(data, x, intermediate=inter)
        Seglist = inter["seglist"]
        Seglist.update_distances(x)
        distance = Seglist.get_closest_segment()
        #
        model_distance = np.array([])
        if add_model_terms:
            model_distance = Seglist.get_closest_data()
        #
        s1, s2 = np.sum(distance**2), np.sum(model_distance**2)
        # print('terms', s1, s2)
        return s1 + s2 

    count = itertools.count()
    def observe(xk):
        n = next(count)
        if n % 10 == 0:
            if not quiet:
                print("{} iterations".format(n))
    
    args = (data_lptr, intermediate)
    _conf = {"tol":1e-2, "method":"BFGS", "callback":observe}
    _conf.update(conf)
    res = scipy.optimize.minimize(loss, x0, args=args, **_conf) 
    return res


def get_result_lptr(res, _model):
    x = res.x
    n = (x.size//2)
    x = x.reshape((n,2))
    _x, _y = x.T
    dt = _model.dt
    return LPtrack(dt, _x, _y)


def local_optimal_lsqmdl(_data, par_value, r, conf={}):
    # print("conf", conf)
    inter = new_intermediate()
    coarsen_parameter = conf.get("coarsen_parameter", 'l')
    _coarse = recursive_coarsen(_data, par_value, parameter=coarsen_parameter)
    print("solving model with M = {}, N = {}".format(len(_coarse),len(_data)))
    with support.PerfTimer() as t:
        res = model_lsq(_data, _coarse, intermediate=inter, conf=conf.get("scipy", {}))
    print("result.success", res.success)
    _model = get_result_lptr(res, _coarse)
    dl = description_length(_data, _model, r, intermediate=inter)
    inter["coarse"] = _coarse
    inter["res"] = res
    inter["exect"] = t.time
    return dl, _model, inter

class IntegerFunction(object):

    def __init__(self, N, objective, args=tuple()):
        # construct
        intr = range(0,N)
        self.llist = llist.dllist([[i] for  i in intr])
        self.args = args
        self.objective = objective
        #
        self.cache_value = {}
        self.cache_model = {}
        self.cache_inter = {}
        # execution time
        self.cache_exect = {}
        # self.cache_res = {}
        # optimisation 
        self.minx = None 


    def __call__(self, x):
        if x in self.cache_value:
            return self.cache_value[x]
        _start = time.perf_counter()
        value, model, inter = self.objective(x, *self.args)
        _end = time.perf_counter()
        self.cache_value[x] = value
        self.cache_model[x] = model
        self.cache_inter[x] = inter
        self.cache_exect[x] = _end-_start
        return value
    
    def pre(self, x):
        x_m = self.llist.nodeat(x).prev.value[0]
        return x_m

    def aft(self, x):
        x_p = self.llist.nodeat(x).next.value[0]
        return x_p

    def contract(self, x, x_m):
        r_node = self.llist.nodeat(x_m)
        node = self.llist.nodeat(x)
        node.value.append(x_m)
        self.llist.remove(r_node)


gr = (np.sqrt(5) + 1) / 2
def gss(f, bracket):
    # https://en.wikipedia.org/wiki/Golden-section_search
    """Golden-section search
    to find the minimum of f on [a,b]
    f: a strictly unimodal function on [a,b]
    """
    a, b = bracket
    _c = b - (b - a) / gr
    _d = a + (b - a) / gr
    c = int(np.floor(_c))
    d = int(np.ceil(_d))
    while (b-a) > 2:
        print('a, b', a,b)
        if f(c) < f(d):
            b = d
        else:
            a = c
        # We recompute both c and d here to avoid loss of precision which may lead to incorrect results or infinite loop
        c = b - (b - a) / gr
        d = a + (b - a) / gr
        c = int(np.floor(c))
        d = int(np.ceil(d))
    # at the end evaluate the triple of adjecent integers and return argmin
    triple = np.arange(a, b+1, 1, dtype=int)
    print("evaluate triple", triple)
    ev = np.array([f(_x) for _x  in triple])
    return triple[np.argmin(ev)]


# test on x^2
def test_gss():
    f = lambda x: (x-3.4)**2
    gss(f, [-10, 10])

def convex_minimise(Intf, bracket):
    a, b = bracket
    # pointless evaluations just so we know the evaulations at the bounds
    Mm = Intf(a)
    Mp = Intf(b)
    return gss(Intf, bracket)


# golden section search
def mdlgss(_data, r, bracket):
    args = (_data, r)
    N = len(_data)
    _conf = {"scipy": {"method":"Nelder-Mead", "tol":1e-2}, "coarsen_parameter": 'M'}
    def loss(_M, _data, r):
        dl, _model, inter = local_optimal_lsqmdl(_data, _M, r, conf=_conf)
        return dl, _model, inter

    Intf = IntegerFunction(N, loss,  args=args)

    xM = convex_minimise(Intf, bracket)
    Intf.minx = xM
    return Intf

#
def best(Intf):
    # select model with the local minimum value of DL but also check for evaluations with lower M
    xM = Intf.minx
    xdl = Intf.cache_value[xM]
    cand = [(m, dl)  for (m, dl) in Intf.cache_value.items() if dl <= xdl]
    return list(sorted(cand))[0]

def new_Mbracket(_data, l0):
    Mmax = len(recursive_coarsen(_data, 0.3*l0))
    Mmin = len(recursive_coarsen(_data, 2*l0))
    bracket = [Mmin, Mmax]
    return bracket

# ----------------------------------------------------------------
# plotting

def new_plot_model_on_data(*args, **kw):
    fig, ax = plt.subplots(figsize=(20,5))
    return plot_model_on_data(ax, *args, **kw)

# note that  we must call plt.draw if we want to use the data transform
# https://stackoverflow.com/questions/29090836/set-aspect-and-coordinate-transforms-in-matplotlib
def plot_model_on_data(ax, model, data, intermediate=None, config={}, lw=8):
    r = config.get("r", 0.03)
    print("plotting using r = ", r)

    def _draw(lw):
        lkw = {"alpha":0.1, "linewidth":lw, "solid_capstyle":'round'}
        highlight_segments = config.get("h_seg", False)
        highlight_outliers =  config.get("h_outlier", False)
        highlight_points =  config.get("h_points", [])
        highlight_nodes =  config.get("h_nodes", [])
        shaded =  config.get("shaded", True)
        labels = []
        handles = []
        ptlkw = {"linestyle":'--', "marker":"o", "alpha":0.4}
        if isinstance(data, LPtrack):
            x = data.x
            y = data.y
        else:
            x = data['x']
            y = data['y']
        # ~~
        l1, = ax.plot(x, y, **ptlkw)
        if highlight_outliers and intermediate is not None and "is_outlier" in intermediate:
            is_outlier = intermediate["is_outlier"]
            outpt = np.stack([x,y],axis=1)[is_outlier]
            for pt in outpt:
                _x, _y = pt
                ax.scatter(_x, _y, color='r')
        #
        for pt_i in highlight_points:
            _x, _y = data.x[pt_i], data.y[pt_i]
            ax.scatter(_x, _y, color="k")
        #
        handles.append(l1)
        labels.append("data")
        # ~~
        if model:
            if highlight_segments:
                color = itertools.cycle(["#FFEC36", "#FB3C2B", "#FF6F00"])
            else:
                color = itertools.cycle(["#FEC216"])
            for i in range(model.x.size-1):
                _x = [model.x[i], model.x[i+1]]
                _y = [model.y[i], model.y[i+1]]
                _col = next(color)
                if hasattr(model, "sidx") and model.sidx[i] in highlight_nodes:
                    green = "#00B41B"
                    _col  = green
                # shaded region
                if shaded:
                    l3, = ax.plot(_x, _y, c=_col, **lkw)
                # center line
                l3, = ax.plot(_x, _y, c=_col, lw=1.5, alpha=0.6)
            handles.append(l3)
            labels.append("model")
        ax.set_aspect('equal')
        return handles, labels

    handles, labels = _draw(lw)
    plt.draw()

    #! why does scaling lw still seem to be bugged?
    #! perhaps we should give up on this idea and re-write this code with patches?

    def scale_lw(ax, r):
        transform = ax.transAxes.transform
        pxsize = transform((1,1)) - transform((0,0))
        lx, ly = ax.get_xlim(), ax.get_ylim()
        datasize = (lx[1]-lx[0], ly[1]-ly[0])
        scale = pxsize[0]/datasize[0]
        # print('pxsize', 'datasize', pxsize, datasize)
        lw = 2 * r * scale
        return lw

    # compute the linewidth in data coordinates
    lw = scale_lw(ax, r)
    print('lw = ', lw)
    ax.clear()
    handles, labels = _draw(lw)
    #
    ax.tick_params(axis='both', labelsize=24)
    ax.legend(handles, labels, fontsize=26)
    # making use of python trick that we can assign arbitrary values to objects
    # ax._handles = handles
    # ax._labels = labels
    
    
def plot_mapped(ax, mapped, _data, _model):
    x, y = mapped.T
    ax.plot(x, y, marker='D', linestyle='None', markersize=4)
    ax.plot(_data.x, _data.y, linestyle='None', marker='o', alpha=0.5)
    ax.plot(_model.x, _model.y, linewidth='4', alpha=0.3, c='k')
    pt = _data.get_n2()
    for i in range(len(_data)):
        pair = np.array([pt[i], mapped[i]])
        _x, _y = pair.T
        ax.plot(_x, _y, linestyle='--', c='orange')
    ax.set_aspect('equal')

# ----------------------------------------------------------------
# text output

def describe_result(res):
    sl = []
    if res.success:
        sl.append( "optimisation finished successfully" )
    else:
        sl.append( "optimisation failed with message {}".format(res.message) )
    sl.append( f"after {res.nit:d} iterations, {res.nfev:d} evaluations" ) 
    return ' '.join(sl)

def describe_optimize_result(inter):
    res = inter["res"]
    exect = inter["exect"]
    sl = []
    if res.success:
        sl.append( "optimisation finished successfully" )
    else:
        sl.append( "optimisation failed" )
    sl.append( f"after {res.nit:d} iterations" ) 
    sl.append( f"({exect:.3f}s)" )
    return ' '.join(sl)