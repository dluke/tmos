
# annealing module

import os, sys
from copy import copy, deepcopy
import random
import pickle
join = os.path.join 
import itertools
import collections
import numpy as np
import scipy
import scipy.optimize

pi = np.pi
norm = np.linalg.norm

import llist
from pili import support
import pili
import mdl

import logging

# important globals

_DT = 0.1
_debug = False

local_f =  sys.stdout

# ~


class AnnealNode(object):

    def __init__(self, dt, x , y, sidx=None):
        self.dt = dt
        self.x = x
        self.y = y
        self.sidx = sidx
        # annealing
        self.temp = 1.0
        # derived member variables
        self.t = None # time coordinate
    
    def get_time(self):
        return _DT * self.t

    def get_priority_task(self):
        return (self.temp, self.sidx)

    def distance(self, other):
        return np.sqrt((self.x - other.x)**2 + (self.y - other.y)**2)

    def __repr__(self):
        return "({}, x {:.3f},y {:.3f}, {})".format(self.dt,self.x,self.y,self.sidx)

    # hashable

    def __hash__(self):
        return self.sidx

    def __eq__(self,other):
        return self.sidx == other.sidx


# ----------------------------------------------------------------
# helpers

#
def _expand_nodes(llpart, include_end=True):
    new = copy(llpart)
    prv = new[0].prev
    nxt = new[-1].next
    if prv:
        new.insert(0, prv)
    if nxt:
        if include_end or nxt.next:
            new.append(nxt)
    return new

def _expand_right(llpart):
    new = copy(llpart)
    nxt = new[-1].next
    if nxt:
        new.append(nxt)
    return new 

def _expand_left(llpart):
    new = copy(llpart)
    prv = new[0].prev
    if prv:
        new.append(prv)
    return new 


# get list of nodes from linked list
def _unroll(llmodel):
    nodes = []
    node = llmodel.first
    while node.next != None:
        nodes.append(node)
        node = node.next
    return nodes

def get_logger_path(name):
    home = join(pili.root, "notebook")
    return join(home, "annealing/{}.log".format(name))



class FLog(object):

    DEBUG = 10
    INFO = 20
    WARNING = 30

    names = {DEBUG: "DEBUG", INFO: "INFO", WARNING: "WARNING"}

    def __init__(self, path):
        self.path = path
        self.level = FLog.DEBUG
        self.form = "{level}: {message}"
        self.f = None

    def open(self):
        self.f = open(self.path, 'w')

    def __str__(self):
        return "Flog(level = {})".format(self.level)

    def set_path(self, path):
        if self.f:
            self.close()
        self.path = path

    def set_level(self, level):
        self.level = level

    def write(self, line):
        print("line", line)
        # ignore logging levels and forms
        if self.f:
            self.f.write(line + "\n")

    def write_array(self, arr, fmt="{:.8f}"):
        form = ' '.join([fmt]*len(arr)) + "\n"
        if self.f:
            self.f.write(form.format(*arr))

    def close(self):
        if self.f:
            self.f.close()

    def clear(self):
        print("clear")
        self.close()
        self.f = open(self.path, 'w')
        self.close()
    
    def log(self, message, level):
        levelname = FLog.names[level]
        message = self.form.format(**{"level": levelname, "message": message})
        if self.level  <= level:
            if self.f:
                self.f.write(message + '\n')

    def debug(self, message):
        self.log(message, FLog.DEBUG)

    def info(self, message):
        self.log(message, FLog.INFO)

# this all doesn't work in global scope

lognames = ["AnnealSegmentList", "Anneal", "Scipy", "State", "Solve"]

ASLLog = FLog(get_logger_path("AnnealSegmentList"))
ASLLog.info("New Logger %s" % ASLLog)

ALog = FLog(get_logger_path("Anneal"))
ALog.info("New Logger %s" % ALog)

ScipyLog = FLog(get_logger_path("Scipy"))
ScipyLog.info("New Logger %s" % ScipyLog)

StateLog = FLog(get_logger_path("State"))
StateLog.info("New Logger %s" % StateLog)

SolveLog = FLog(get_logger_path("Solve"))
SolveLog.info("New Logger %s" % SolveLog)

logs = [ASLLog, ALog, ScipyLog, SolveLog, StateLog]

def set_logging_directory(targetdir):
    for i, log in enumerate(logs):
        name = lognames[i]
        log.set_path(join(targetdir, name + '.log'))

def set_logging_level(level):
    for log in logs:
        log.set_level(level)

def clear_logs():
    for log in logs:
        log.clear()

def close_logs():
    for log in logs:
        log.close()

# util

def nodeiter(llmodel):
    node = llmodel.first
    while node != None:
        yield node
        node = node.next

# ----------------------------------------------------------------
# similar to mdl.SegmentList class but I think it's easier to re-write than to generalise
class AnnealSegmentList(object):


    def __init__(self, data):
        N = len(data)
        self.N = N
        self.point_data = data.get_n2()
        self.distance_matrix = [{} for _ in range(N)] 
        #  ^ (data pt, segments)
        self.s_matrix = [{} for _ in range(N)] 
        # ^ s is the local "coordinate" of the point mapped onto the line
        self.closest_segment = np.empty(N, dtype=int)
        self.segment_distance = np.empty(N)
        self.local_coord = np.empty(N)

        # flags
        self.allow_remapping = True
        self.mapping_constraint = True
        
        # experiment
        self.remap = {}

    def clear_distance_matrix(self):
        self.distance_matrix = [{} for _ in range(self.N)] 
        self.s_matrix = [{} for _ in range(self.N)] 

    def initialise(self, model):
        # self.model = model
        self.M = len(model)
        # ^ model comes with an initial time distribution
        # ^ needs to be updated ...
        self._init(model)
        # self.update_at(list(nodeiter(model)))

    def _init(self, model):
        # ! wont this fail for self intersecting trajectory?
        # todo we need to use update_at() instead but there are some sublte bugs to work through
        # M x N
        lptr = array_model(model)
        M = len(lptr)
        xy = np.column_stack([lptr.x, lptr.y])
        a = xy[:-1]
        b = xy[1:]
        time = np.zeros(M, dtype=int)

        for pt_i in range(self.N):
            p = self.point_data[pt_i]
            p = p[:,np.newaxis].T
            dist = np.zeros(M-1)
            pwls = np.zeros(M-1)
            for i in range(M-1):
                _a = a[i].reshape(1,2)
                _b = b[i].reshape(1,2)
                s, d = support.line_coord_dist(p, _a, _b)
                dist[i] = d[0]
                pwls[i] = s[0]
            #
            min_k = np.argmin(dist)
            time[min_k] += 1
            min_d = dist[min_k]
            sidx = lptr.sidx[min_k]
            self.closest_segment[pt_i] = sidx
            self.segment_distance[pt_i] = min_d
            self.local_coord[pt_i] = pwls[min_k]

        # update dt (segment counts)
        node = model.first
        node.value.dt = 0
        i = 0
        while node.next:
            node = node.next
            node.value.dt = time[i]
            i += 1
        self.update_node_time_at(node)

    
    @staticmethod
    def _current_time(node):
        # O(n)
        dt = 0
        _node = node
        while _node:
            dt += _node.value.dt
            _node = _node.prev
        return dt
        
    def update_node_time_at(self, node):
        model = node.owner()
        t = 0
        _node =  model.first
        _node.value.t = t 
        while _node.next:
            _node = _node.next
            t += _node.value.dt
            _node.value.t = t
        assert(_node.value.t == self.N)

    def minargmin(self, exllpart):
        """
        # distance
        # closest segement idx
        # PWL local coordinate
        """
        
        m = len(exllpart)
        leftnode = exllpart[0]
        rightnode = exllpart[-1]
        record_sidx = []

        # extract the moddable chunk 
        part_ftidx = leftnode.value.t
        part_ltidx = rightnode.value.t
        part_sidx = [_node.value.sidx for _node in exllpart[:-1]]
        ijmap = np.copy(self.closest_segment[part_ftidx:part_ltidx])
        newmap = np.empty_like(ijmap)
        mindist = np.empty(ijmap.size)
        scoord = np.empty(ijmap.size)
        
        remapping = []
        for i in range(m-1):
            node = exllpart[i]

            # expand nodes
            adj_nodes = []
            if node != leftnode and node.prev:
                adj_nodes.append(node.prev)
            adj_nodes.append(node)
            if  node.next != rightnode:
                adj_nodes.append(node.next)

            slidx = [_node.value.sidx for _node in adj_nodes]
            ftidx, ltidx = node.value.t, node.next.value.t
            # points in the locality of chosen segment
            for pt_i in range(ftidx, ltidx):
                # As an experiment: remap a specific point to a specific segment
                if self.allow_remapping:
                    if pt_i in self.remap:
                        # !
                        # print(slidx, pt_i, self.remap)
                        p = self.point_data[pt_i].reshape(1,2)
                        node = self.remap[pt_i]
                        sidx = node.value.sidx
                        _a = np.array([node.value.x, node.value.y]).reshape((1,2))
                        _b = np.array([node.next.value.x, node.next.value.y]).reshape((1,2))
                        k = sidx
                        s, d = support.line_coord_dist(p, _a, _b)
                    else:
                        enum_slidx = enumerate(slidx)
                        d = np.inf
                        k = None
                        s = None
                        # chosen and adjacent segments
                        for s_i, sidx in enum_slidx:
                            #!  all the distances should be here # or not?
                            # d = self.distance_matrix[pt_i][sidx]
                            # s = self.s_matrix[pt_i][sidx]
                            # todo : efficiency
                            sidx_node = self.sidxmap[sidx]
                            _s, _d = self.distance_matrix_element(pt_i, sidx_node)

                            if _d < d:
                                k = sidx
                                d = _d
                                s = _s
                    if k == None:
                        raise ValueError("sidx is None")
                    local_index = pt_i - part_ftidx
                    newmap[local_index] = k
                    if self.mapping_constraint:
                        if newmap[local_index] != ijmap[local_index]:
                            remapping.append((local_index, ijmap[local_index], newmap[local_index]))
                    mindist[local_index] = d
                    scoord[local_index] = s
                    record_sidx.append(k)
                else:
                    #! prevent remapping : for debugging
                    sidx = self.closest_segment[pt_i]
                    mindist[local_index] = self.distance_matrix[pt_i][sidx]
                    scoord[local_index] = self.s_matrix[pt_i][sidx]
                    record_sidx.append(sidx)

        # counts = {node.value.sidx : node.next.value.dt for node in _expand_nodes(exllpart[:-1], include_end=False)}
        # counts = {node.value.sidx : node.next.value.dt for node in exllpart[:-1]}
        # !tmp O(n)
        if self.mapping_constraint:
            counts = collections.OrderedDict([(node.value.sidx , node.next.value.dt) for node in nodeiter(exllpart[0].owner()) if node.next])

        _stop = False
        if _debug:
            copycounts = counts.copy()
            copymap = ijmap.copy()

        # constrain segments to have at least two associated points
        if self.mapping_constraint:
            for remap in remapping:
                local_index, sj, sjprime = remap
                pt_i = local_index + part_ftidx
                if counts[sj] > 2:
                    # allow remap
                    # print('map point {}, {} -> {}'.format(pt_i, sj, sjprime))
                    ijmap[local_index] = sjprime

                    if counts[sj] < 1:
                        print('count problem', counts[sj])
                        print('map point {}, {} -> {}'.format(pt_i, sj, sjprime))
                        print()
                        _stop = True

                    counts[sj] -= 1
                    counts[sjprime] += 1
                else:
                    # don't allow remap
                    # in this case we need to recompute the d/s for the old mapping
                    node = self.sidxmap[sj]
                    s, d = self.distance_matrix_element(pt_i, node)
                    mindist[local_index] = d
                    scoord[local_index] = s
        else:
            ijmap = newmap


        if _debug and _stop:
            print('ft, lt', part_ftidx, part_ltidx)
            print('counts before', copycounts, file=f)
            print('part_sidx', part_sidx, file=f)
            print('sidx', [v.sidx for v in exllpart[0].owner()], file=f)
            print('oldmap', self.closest_segment, file=f)
            print('modmap', ijmap, file=f)
            print('counts after', counts, file=f)
            print(file=f)

        if _debug and _stop:
            sys.exit()
        # ! put this assert back in
        # assert(min([counts[sidx] for sidx in part_sidx]) > 1)

        # replace the moddable chunk
        self.closest_segment[part_ftidx:part_ltidx] = ijmap
        self.segment_distance[part_ftidx:part_ltidx] = mindist
        self.local_coord[part_ftidx:part_ltidx] = scoord


        # update 'time' (counting points adjacent to segments)
        # !tmp O(n)
        if self.mapping_constraint:
            for node in nodeiter(exllpart[0].owner()):
                if node.next:
                    assert(node.value.sidx in counts)
                    node.next.value.dt = counts.get(node.value.sidx)
        else:
            counts = collections.Counter(record_sidx)
            for node in exllpart[:-1]:
                node.next.value.dt = counts.get(node.value.sidx, 0)

        # debugging
        if _debug:
            N = sum([v.dt for v in exllpart[0].owner()])
            assert(N == self.N)

    def distance_matrix_element(self, pt_i, node):
        sidx = node.value.sidx
        s, d = self.s_matrix[pt_i].get(sidx), self.distance_matrix[pt_i].get(sidx)
        if s == None or d == None:
            if _debug:
                print('warning: distance not found in matrix for ', pt_i, sidx)
            v, vn = node.value, node.next.value
            a = np.array([[v.x, v.y]])
            b = np.array([[vn.x, vn.y]])
            p = self.point_data[pt_i].reshape((1,2))
            s, d = support.line_coord_dist(p, a, b)
        return s, d

    def get_remapping(self, ijmap, newmap):
        index = np.argwhere(ijmap != newmap).ravel()
        dct = {}
        for i in index:
            tup = (i, ijmap[i], newmap[i])
            if ijmap[i] in dct:
                dct[ijmap[i]].append( tup )
            else:
                dct[ijmap[i]] = [tup]
        return dct


    def get_model_distance(self):
        llpart = self._llpart[1:-1]
        tidx = self._tidx
        m_distance = []
        for i, node in  enumerate(llpart):
            _d = node.value
            m_pt = np.array([_d.x, _d.y])
            ftidx, ltidx = tidx[i], tidx[i+2]
            p = self.point_data[ftidx:ltidx]
            dist = norm(p - m_pt, axis=1)
            d = np.min(dist)
            m_distance.append(d)
        return np.array(m_distance)


    # todo : optimise : split linked list and array representations (speed up when not adding/removing points)
    def update_at(self, llpart):
        self.update_node_time_at(llpart[-1])
        self.update_distance_at(llpart)

    def update_distance_at(self, llpart): 
        """
        argument: moveable nodes
        first expand: distances to update
        second expand: distances to include in objective function
        """
        # first expand includes all the moving segments
        llpart = _expand_nodes(llpart.copy())
        # second expand includes adjacent affected segments
        exllpart = _expand_nodes(llpart.copy())
        # todo : need to save these?
        self._llpart = llpart
        self._exllpart = exllpart
        
        # get segments
        ab = np.array([[node.value.x,node.value.y] for node in exllpart])
        a = ab[:-1]
        b = ab[1:]

        # compute sparse matrices
        self.clear_distance_matrix()
        m = len(exllpart)
        for i in range(m-1):
            node = exllpart[i]
            sidx = node.value.sidx
            adj_nodes = _expand_nodes([node])
            rightnode = adj_nodes[-1].next if adj_nodes[-1].next != None else adj_nodes[-1]
            ftidx, ltidx = adj_nodes[0].value.t, rightnode.value.t
            p = self.point_data[ftidx:ltidx]
            _a = a[i].reshape(1,2)
            _b = b[i].reshape(1,2)
            _s, _d = support.line_coord_dist(p, _a, _b)
            if (ftidx < 0 or ltidx > self.N):
                print(list(nodeiter(exllpart[0].owner())))
                print("fill matrix for (ftidx, ltidx, sidx) ", ftidx, ltidx, sidx, a.shape, b.shape, p.shape, _s.shape, _d.shape, file=f)
                sys.exit()

            for i, pt_i in enumerate(range(ftidx, ltidx)):
                self.distance_matrix[pt_i][sidx] = _d[i]
                self.s_matrix[pt_i][sidx] = _s[i]
        
        #
        self.minargmin(exllpart)

    def get_segment_distance(self):
        return self.segment_distance

    # O(n)
    # * compute the pwl curve coordinate 't'
    # * using this definition modifying the length of partial segment affects all the coordinates ahead of it
    # * therefor the computation is always O(n)
    def get_local_coord(self):
        return self.local_coord


    def __getstate__(self):
        state = self.__dict__.copy()
        # don't deepcopy temporaries
        if "_llpart" in state: del state["_llpart"]
        if "_exllpart" in state: del state["_exllpart"]
        del state['sidxmap']
        return state



def linked_model(model):
    model = deepcopy(model)
    model.dt = np.array([int(round(_t/_DT)) for  _t in model.dt])
    # ^ convert time to an index ( model DT = 1 )
    _seq = [AnnealNode(*_tuple) for _tuple in zip(model.dt, model.x, model.y, model.sidx)]
    return llist.dllist(_seq)

def partial_array_model(nodes):
    mdata = np.empty((len(nodes),4))
    for i, node in enumerate(nodes):
        v = node.value
        mdata[i] = np.array([v.dt, v.x, v.y, v.sidx])
    dt, x, y, sidx = mdata.T
    dt = _DT * dt.astype(float) # return to the original time units
    sidx = sidx.astype(int)
    return mdl.LPtrack(dt, x, y, sidx)

def array_model(llmodel):
    M = len(llmodel)
    mdata = np.empty((M,4))
    for i, nodedata in enumerate(llmodel):
        mdata[i] = np.array([nodedata.dt,nodedata.x,nodedata.y,nodedata.sidx])
    dt, x, y, sidx = mdata.T
    dt = _DT * dt.astype(float) # return to the original time units
    sidx = sidx.astype(int)
    return mdl.LPtrack(dt, x, y, sidx)

def new_sidx_map(llmodel):
    # construct mapping {sidx : node} from linked list model
    sidxmap = {}
    node = llmodel.first
    while node:
        sidxmap[node.value.sidx] = node
        node = node.next
    return sidxmap


# Loss function

class Loss(object): 

    def __init__(self, anneal, seglist, llpart, 
            use_ordering=True, 
            snap_model=False,
            contour_term=0,
            debug = False
            ):
        self.anneal = anneal
        self.seglist = seglist
        self.llpart = llpart 
        self.use_ordering = use_ordering
        self.snap_model = snap_model
        self.contour_term  = contour_term
        self.debug = debug
        self.terms = []

        # clear
        path = "annealing/Loss.log"
        if os.path.exists(path):
            os.remove(path)
            

        # tmp
    def log(self, line):
        with open("annealing/Loss.log", 'a+') as f:
            f.write(line)

    def __getstate__(self):
        # Remove the unpicklable entries.
        state = self.__dict__.copy()
        if 'llpart' in state: del state['llpart']
        return state

    def _update_llmodel(self, xy):
        llpart = self.llpart
        for i in range(len(xy)):
            nodedata = llpart[i].value
            nodedata.x, nodedata.y = xy[i]
            # logging.debug("update {}, {} -> {}".format(i, ab[i], xy[i]))

    def ordering_loss(self):
        # todo: make this local (implemented in seglist)
        pwlt = self.anneal.get_pwl_coord()
        pwlpair = list(enumerate(pwlt))
        sortpwl = sorted(pwlpair, key=lambda t:t[1])
        sort_s = np.array([t[1] for t in sortpwl])
        residuals = (sort_s -  pwlt)**2
        return residuals

    def model_loss(self):
        # snap the model to the data using the local coordinate
        # print("model_loss")
        # print("length, np.min(s_adjacent), np.max(s_adjacent), sqs1, sqs2")
        
        # optional function gives zero contribution for deviations less than r
        # kr = 2 * self.anneal.r
        kr = 0
        def offset_func(a): 
            # assume a > 0
            if (a < kr):
                return 0
            return (a - kr)**2
            
        slocal = self.seglist.get_local_coord()
        part = _expand_left(self.llpart)
        # iterate the segments
        term = 0
        for node in part:
            v = node.value
            vn = node.next.value
            ftidx, ltidx = node.value.t, node.next.value.t
            s_adjacent = slocal[ftidx:ltidx]
            if len(s_adjacent) == 0:
                continue
            length = np.sqrt((v.x - vn.x)**2 + (v.y - vn.y)**2)

            # special case for one adjacent point
            if len(s_adjacent) == 1:
                x = slocal[0]
                if x > length/2:
                    term += offset_func(length - x)
                else: 
                    term += offset_func(x)
                continue

            # 
            s1 = np.min(s_adjacent)
            s2 = length - np.max(s_adjacent)
            sqs1 = offset_func(s1)
            sqs2 = offset_func(s2)
            term += sqs1 + sqs2
            # print("({:.3f},{:.3f}) -> ({:.3f}, {:.3f})".format(v.x, v.y, vn.x, vn.y))
            # print(length, np.min(s_adjacent), np.max(s_adjacent), sqs1, sqs2)
        # print()
        return term

    def contour_length(self):
        part = _expand_left(self.llpart)
        # convert to array
        xy = np.empty((len(part),2))
        for i, node in enumerate(part):
            v = node.value
            xy[i][0] = v.x
            xy[i][1] = v.y
        return np.sum(norm(xy[1:] - xy[:-1], axis=1))

    def __call__(self, x):
        seglist = self.seglist
        llpart = self.llpart
        # update in callback?
        m = len(llpart)
        xy = x.reshape((m,2))
        self._update_llmodel(xy)
        seglist.update_at(llpart) #! node.value.t coordinate is updated
        s1 = np.sum(seglist.segment_distance**2)
        
        s2 = 0.
        if self.use_ordering:
            s2 = np.sum(self.ordering_loss())
        
        s3 = 0.
        if self.snap_model:
            s3 = self.model_loss()
            
        s4 = 0.
        if self.contour_term > 0:
            s4 = self.contour_term * self.contour_length()

        if self.debug:
            line = "{:.8f} {:.8f} {:.8f} {:.8f}\n".format(s1, s2, s3, s4)
            self.log(line)

        self.terms = [s1, s2, s3, s4]
        return sum(self.terms)



class Anneal(object):

    def __init__(self, r=0.03, k=1./160, use_description_length=True):
        self.r = r
        self.k = k 
        self.loss = None
        self.default_loss_conf = {}
        self.use_description_length = use_description_length

    def set_mapping_constraint(self, value):
        self.seglist.mapping_constraint = value

    def get_sidxmap(self):
        return self.seglist.sidxmap
    def set_sidxmap(self, _sidxmap):
        self.seglist.sidxmap = _sidxmap
    sidxmap = property(get_sidxmap, set_sidxmap)

    def initialise(self, model, data):
        """
        model is LPtrack
        data is LPtrack
        """
        self.M = len(model)
        self.N = len(data)
        # ALog.info("Anneal (M, N) = ({},{})".format(self.M, self.N))

        self.llmodel = linked_model(model)
        # ^ construct a mapping {sidx : node}
        # this unables O(1) node retrieval, the cost is we need to copy the sidx_map along with the linked list
        # ! the state of these objects can't be pickled
         
        _start_temp = 1
        self.seglist = AnnealSegmentList(data)
        self.seglist.initialise(self.llmodel)
        # ^ initialise distance matrix
        #
        self.count_segments = itertools.count(self.M - 1)
        # ^ setup segment counter

        self.sidxmap = new_sidx_map(self.llmodel)

        x, y =  data.x, data.y
        pad = 2 * self.r
        rxl = (x.min()-pad, x.max()+pad)
        ryl = (y.min()-pad, y.max()+pad)
        self.point_data_lims = (rxl, ryl)
        # ^ setup model limits

    # do not modify the linked list model except by these methods
    def insertafter(self, nodedata, node):
        new = self.llmodel.insertafter(nodedata, node)
        self.sidxmap[nodedata.sidx] = new
        self.M += 1
        return new

    def remove(self, node):
        self.llmodel.remove(node)
        del self.sidxmap[node.value.sidx]
        self.M -= 1

    def modify_sidx(self, node, sidx):
        _sidx = node.value.sidx
        node.value.sidx = sidx
        self.sidxmap[sidx] = node
        del self.sidxmap[_sidx]

    def _next_sidx(self):
        return next(self.count_segments)

    # ~

    def get_mapping_counts(self):
        counts = collections.OrderedDict([(node.value.sidx , node.next.value.dt) for node in nodeiter(self.llmodel) if node.next])
        return counts

    def accept_probability(self):
        pass

    def sample(self):
        # randomly choose a new state
        pass

    def destroy(self, node):
        # * destroy two segments and create one new 
        # !tmp
        # print('destroy', node.value.sidx)
        # print('mapping before', self.seglist.closest_segment)
        
        llmodel = self.llmodel
        next_node = node.next
        prev_node = node.prev
        if node.prev is None:
            node.next.next.value.dt += node.next.value.dt
            node.next.value.dt = 0
            affected = [next_node]
            scs = self.seglist.closest_segment
            scs[scs == node.value.sidx] = node.next.value.sidx

        elif node.next is None:
            raise ValueError("cannot destroy end node")

        else:
            new_sidx = self._next_sidx()
            affected = [prev_node, next_node]
            scs = self.seglist.closest_segment
            scs[scs == node.value.sidx] = new_sidx
            scs[scs == node.prev.value.sidx] = new_sidx

            self.modify_sidx(node.prev, new_sidx)  
            node.next.value.dt += node.value.dt

        self.remove(node)

        # self.seglist.update_at(affected)
        self.seglist.update_node_time_at(next_node)

        # # !tmp
        # a = [v.dt for v in next_node.owner()]
        # t = [v.t for v in next_node.owner()]
        # print('after destroy', sum(a), a)
        # print(t)
        # print('mapping after', self.seglist.closest_segment)
        N = sum([v.dt for v in next_node.owner()])
        assert(N == self.N)

        # self.seglist.initialise(self.llmodel) # !tmp
        return affected

    def create(self, node):
        # * destroy one segment and create two new ones, equally spaced
        llmodel = self.llmodel
        if node.next is None:
            raise ValueError("cannot create end node")
        next_node = node.next 
        
        #  new position
        _data = node.value
        x = (_data.x + node.next.value.x)/2
        y = (_data.y + node.next.value.y)/2

        # update adjacent point counts
        _hdt = node.next.value.dt//2
        # assert(_hdt >= 2)
        node.next.value.dt -= _hdt
        
        new_left_sidx = self._next_sidx()
        new_right_sidx = self._next_sidx()

        # update mapping
        ft = node.value.t
        scs = self.seglist.closest_segment
        _left = scs[:ft+_hdt]
        _left[_left == node.value.sidx] = new_left_sidx
        _right = scs[ft+_hdt:]
        _right[_right == node.value.sidx] = new_right_sidx
        scs[:ft+_hdt] = _left
        scs[ft+_hdt:] = _right
    
        # create the new node & update sidx
        nodedata = AnnealNode(_hdt, x, y)
        nodedata.sidx = new_right_sidx
        new_node = self.insertafter(nodedata, node)
        self.modify_sidx(node, new_left_sidx)

        assert(next_node.value.dt >= 0)
        assert(new_node.value.dt >= 0)

        N = sum([v.dt for v in node.owner()])
        assert(N == self.N)

        # ^ node also has a new segment
        affected = [node, new_node, next_node]
        # self.seglist.initialise(self.llmodel) # !tmp
        self.seglist.update_at(affected)
        return affected

    def loss_function(self, nodes, n_local, use_ordering=True, loss_conf=None):
        if not loss_conf:
            loss_conf = self.default_loss_conf
        for _ in range(n_local):
            nodes = _expand_nodes(nodes.copy())
        llpart = nodes
        args = (self, self.seglist, llpart)
        loss = Loss(*args, use_ordering=use_ordering, **loss_conf)
        ab = np.array([[node.value.x,node.value.y] for node in llpart])
        x0 = ab.flatten()
        loss(x0)
        return loss

    def local_lsq_optimise(self, nodes, 
            n_local=1, use_bounds=True, use_ordering=True, 
            options={}, optconf={}, loss_conf=None, _verbose=False):
        # _verbose = True # !tmp
        if not loss_conf:
            loss_conf = self.default_loss_conf
        """
        optconf: arguments to  scipy.optimise.minimize
        options: scipy.optimize.minimise options argument (could be nested in optconf)
        config: configure this method
        loss_conf: the configuration for the loss function
        """
        # optimise locally around nodes
        # 
        optconf["method"] = "Nelder-Mead"
        # if "maxiter" not in options: options["maxiter"] = 100
        if _verbose:
            print("local lsq optimise using options: ", options)
        #
        for _ in range(n_local):
            nodes = _expand_nodes(nodes.copy())
        llpart = nodes
        if _verbose:
            print("free {} nodes to optimise".format(len(llpart)) )
        ab = np.array([[node.value.x,node.value.y] for node in llpart])
        x0 = ab.flatten()

        # bounds
        bounds = None
        if use_bounds:
            rxl, ryl = self.point_data_lims
            n = len(x0)
            bounds = [None for _ in range(n)]
            bounds[0::2] = [rxl for _ in range(n//2)]
            bounds[1::2] = [ryl for _ in range(n//2)]
        optconf["bounds"] = bounds
        #

        args = (self, self.seglist, llpart)
        loss = Loss(*args, use_ordering=use_ordering, **loss_conf)
        res = scipy.optimize.minimize(loss, x0, options=options, **optconf)
        print(mdl.describe_result(res), file=local_f)
        if _verbose:
            print(res.message)
            print("nit, nfev",  res.nit, res.nfev)
        self.loss = loss
        return res.success
    
    def get_loss(self):
        if self.loss:
            return self.loss
        else:
            raise RuntimeError("No loss function available")

    def get_residuals(self):
        # return the residual array, size N
        return self.seglist.get_segment_distance()**2

    def get_outliers(self):
        distance = self.seglist.get_segment_distance()
        outlier = distance > self.r
        return outlier

    def get_score(self):
        # generic alternative to description length
        # todo
        if self.use_description_length:
            return self.get_description_length()
        else:
            #! this idea appears to be a failure
            lsq = sum(self.loss.terms)
            return lsq + self.k * self.r**2 * self.M


    def get_description_length(self):
        # todo : optimise by local updating
        distance = self.seglist.get_segment_distance()
        outlier = distance > self.r
        self._outlier = outlier
        n = int(np.sum(outlier))
        return self.M + n

    def get_pwl_coord(self, get_mapped=False):
        # TODO: make this calculation local (?)
        """
        get_mapped: 
            additionally return the mapped xy positions
        return (N,) mapped 'time' coordinates 
        """
        model = self.llmodel
        scoord = self.seglist.get_local_coord()
        txy = self._get_mdata()
        # ^ slow
        xy = txy[:,1:3]
        model_disp = xy[1:] - xy[:-1]
        segment_length = norm(model_disp, axis=1)
        i_length  = np.insert(np.cumsum(segment_length)[:-1], 0, 0)
        # we mapped points to segments already, now we need to map segment ids on to their order
        order_map = {}
        for i, nodedata in enumerate(model):
            order_map[nodedata.sidx] = i
        pwl = []
        mapped = []
        for pt_i in range(self.N):
            sidx = self.seglist.closest_segment[pt_i]
            _si = order_map[sidx]
            il = i_length[_si]
            s = scoord[pt_i]
            coord = il + s
            pwl.append(coord)
            if get_mapped:
                d = segment_length[_si]
                mpt = xy[_si] + model_disp[_si] * s/d
                mapped.append(mpt) 
        # return
        if get_mapped:
            return np.array(pwl), np.array(mapped)
        else:
            return np.array(pwl)

    def _get_mdata(self):
        mdata = np.empty((self.M,4))
        for i, nodedata in enumerate(self.llmodel):
            mdata[i] = np.array([nodedata.dt,nodedata.x,nodedata.y,nodedata.sidx])
        return mdata

    def _update_llmodel(self, mdata):
        for i, nodedata in enumerate(self.llmodel):
            nodedata.dt, nodedata.x, nodedata.y = mdata[i][:3]
            nodedata.sidx = int(mdata[3])

    def get_segment_node(self, sidx):
        try:
            return self.sidxmap[sidx]
        except KeyError:
            raise KeyError("segment {} not in model".format(sidx))


    def expand_nodes(self, nodes, n_local=1):
        for _i in range(n_local):
            nodes = _expand_nodes(nodes)
        return nodes

    def get_nodes_at(self, sidx, n_local=1):
        # get a number of nodes before and after the target node
        node = self.get_segment_node(sidx)
        prev = []
        aft = []
        _node = node
        for _i in range(n_local):
            _p = _node.prev
            if _p == None:
                break
            prev.append(_p)
            _node = _p
        _node = node
        for _i in range(n_local):
            _p = _node.next
            if _p == None:
                break
            aft.append(_p)
            _node = _p
        return [*prev, node, *aft]
        
    def get_current_model(self):
        return array_model(self.llmodel) 

    # https://docs.python.org/3/library/pickle.html#what-can-be-pickled-and-unpickled
    # used by deepcopy
    def __getstate__(self):
            # Copy the object's state from self.__dict__ which contains
            # all our instance attributes. Always use the dict.copy()
            # method to avoid modifying the original state.
            state = self.__dict__.copy()
            # Remove the unpicklable entries.
            del state['llmodel']
            return state


    def dump_state(self, path="annealing/current"):
        annealpath = '_'.join([path, 'anneal.pkl'])
        modelpath = '_'.join([path, 'model.pkl'])
        with open(annealpath, 'wb') as f:
            print("dump state to", annealpath)
            pickle.dump(self, f)
        with open(modelpath,'wb') as f:
            print("dump model to", annealpath)
            # ! this is a trick we use to pickle the model
            model = self.get_current_model()
            pickle.dump(model, f)

    @classmethod
    def load_state(cls, path="annealing/current"):
        annealpath = '_'.join([path, 'anneal.pkl'])
        modelpath = '_'.join([path, 'model.pkl'])
        with open(annealpath, 'rb') as f:
            new = pickle.load(f)
        with open(modelpath,'rb') as f:
            model = pickle.load(f)
        new.set_model(model)
        return new
    
    def set_model(self, model):
        self.llmodel = linked_model(model) 
        self.sidxmap = new_sidx_map(self.llmodel)

    def clone(self):
        new = deepcopy(self)
        # explicitely create a new linked list because it is not copiable
        new.llmodel = llist.dllist([deepcopy(_n) for _n in self.llmodel])
        new.sidxmap = new_sidx_map(new.llmodel)
        return new
        

# ----------------------------------------------------------------
# least squares solve by solving small chunks of trajectory


class BaseSolver(object):

    def __init__(self, anneal, rng=None, config={}):
        self.anneal = anneal
        self.rng = rng if rng!=None else np.random.RandomState(0)
        self.config = config

    def dump_state(self, path="annealing/current"):
        # todo record the state of the solver, not just some of its components
        rngpath = '_'.join([path, "rng.pkl"])
        self.anneal.dump_state(path)
        with open(rngpath, 'wb') as f:
            print("dump random state  to", rngpath)
            pickle.dump(self.rng.get_state(), f)
        solvepath = '_'.join([path, "solver.pkl"])
        with open(solvepath, 'wb') as f:
            print("dump solver state  to", solvepath)
            pickle.dump(self,f )

    def __getstate__(self):
        state = self.__dict__.copy()
        if "anneal" in state: del state["anneal"]
        if "rng" in state: del state["rng"]
        return state

    @classmethod
    def load_state(cls, path="annealing/current"):
        anneal = Anneal.load_state(path)
        rngpath = '_'.join([path, "rng.pkl"])
        with open(rngpath, 'rb') as f:
            rng = pickle.load(f)
            rngstate = np.random.RandomState()
            rngstate.set_state(rng)
        solverpath = '_'.join([path, "solver.pkl"])
        with open(solverpath, 'rb') as f:
            new = pickle.load(f)
        new.rng = rngstate
        new.anneal = anneal
        # new = cls(anneal, rng=rngstate)
        return new

    def get_current_model(self):
        return self.anneal.get_current_model()

    def default_local_lsq(self, anneal, nodes, n_local=1):
        use_ordering = self.config.get('use_ordering', True)
        nm_options = {"fatol":0.01, "xatol":0.01}
        return anneal.local_lsq_optimise(nodes, n_local=n_local, options=nm_options, use_ordering=use_ordering)

class ChunkLsq(BaseSolver):
    
    def linear_solve(self, **opt):
        n_local = opt.get("n_local", 2)
        nm_options = {"fatol":0.01, "xatol":0.01}
        node = self.anneal.llmodel.first.next
        success = False
        while node.next !=  None:
            success = self.anneal.local_lsq_optimise([node], n_local=n_local, options=nm_options, use_bounds=True)
            # success = self.default_local_lsq(self.anneal, [node])
            if not success:
                pass
            node = node.next
        return success

    def multiple_linear_solve(self, n=1, **opt):
        for _i in range(n):
            self.linear_solve(**opt)

# --------------------------------
# solver with create and destroy function

class Solver(ChunkLsq):

    def __init__(self, anneal, rng=None, config={}):
        super().__init__(anneal, rng, config=config)
        self.reset() # todo breaks dump/load solver states
        self.greedy = config.get("greedy", True)
    
    def reset(self):
        segmentidx = [_n.sidx for _n in self.anneal.llmodel][:-1] # last node has no segment
        # segmentidx = segmentidx[1:-1] # discard end node
        self.deleted = set()
        self.rng.shuffle(segmentidx) # shuffle in place
        self.queue = collections.deque(segmentidx)
        self.use_priority = False
        # observer
        self._record = []

    def new_priority_task(self, node, priority):
        temp = 1.0 - priority
        node.value.temp = temp
        if node.value.sidx not in self.pqueue:
            # print("new task", (node.value, priority))
            self.pqueue.add_entry([priority, node.value.sidx])
        else:
            # print("update task", (node.value, priority))
            self.pqueue.update_priority([None, node.value.sidx], priority)
    
    def random_accept(self, score_prime, score, T):
        if self.greedy:
            return score_prime < score
        else:
            s = np.exp(-(score_prime - score)/T) 
            u = self.rng.uniform()
            # print('s, u, T', s, u ,T)

            return s >= u
        return 

    def create(self, sidx):
        # initial state
        node = self.anneal.get_segment_node(sidx)
        self.anneal.loss_function([node], n_local=1) #! must be the same as default_local_lsq
        score = self.anneal.get_score()

        # compute final state
        _anneal = self.anneal.clone() #! O(n)
        _node = _anneal.get_segment_node(sidx)
        if _node.next == None: 
            return
        if self.anneal.seglist.mapping_constraint and node.next.value.dt < 4:
            return

        affected = _anneal.create(_node)
        self.default_local_lsq(_anneal, affected)
        score_prime = _anneal.get_score()

        # bookkeeping
        if self.random_accept(score_prime, score, _node.value.temp):
            accepted = True
            print(f"accept create segment at {sidx}", file=local_f)
            # accept
            self.anneal = _anneal
            for _node in affected:
                _temp = 1.0
                self.new_priority_task(_node, 1.0 - _temp)
            self.deleted.add(sidx)
        else:
            accepted = False
            # reject change
            _temp = self.fc * node.value.temp
            self.new_priority_task(node, 1.0 - _temp)
        return accepted
            
    def destroy(self, sidx):
        # initial state
        node = self.anneal.get_segment_node(sidx)
        self.anneal.loss_function([node], n_local=1) #! must be the same as default_local_lsq
        score = self.anneal.get_score()

        # compute final state
        _anneal = self.anneal.clone() #! O(n)
        _node = _anneal.get_segment_node(sidx)
        affected = _anneal.destroy(_node)
        self.default_local_lsq(_anneal, affected)
        score_prime = _anneal.get_score()

        # bookkeeping
        if self.random_accept(score_prime, score, _node.value.temp):
            # accept
            accepted = True
            print(f"accept destroy segment {sidx}", file=local_f)
            self.anneal = _anneal
            for _node in affected:
                _temp = 1.0
                self.new_priority_task(_node, 1.0 - _temp)
            self.deleted.add(sidx)
            if node.prev != None: 
                self.deleted.add(node.prev.value.sidx)
        else:
            accepted = False
            # reject change
            _temp = self.fd * node.value.temp
            self.new_priority_task(node, 1.0 - _temp)
        return accepted

    def random_transition(self, sidx):
        # choose create or destroy with equal probability
        options = [('create', self.create), ('destroy', self.destroy)]
        i = self.rng.choice([0,1])
        name, transition_method = options[i]

        print(f'random {name} transition at {sidx}', file=local_f)
        accepted = transition_method(sidx)
        self._record.append((sidx, name, accepted))

    def create_and_destroy(self, sidx):
        # * try creating, then try destroying, then return the anneal object with the best DL, M 
        _anneal = self.anneal
        # SolveLog.debug("init model {}".format(self.anneal.llmodel))
        _node = _anneal.get_segment_node(sidx)
        _anneal.loss_function([_node], n_local=1) #! must be the same as default_local_lsq
        init_dl = _anneal.get_score()
        store = [_anneal]
        store_affected = [[self.anneal.get_segment_node(sidx)]]

        # ! create
        _anneal = self.anneal.clone()
        _node = _anneal.get_segment_node(sidx)
        if _node.next == None: # end node test
            return
        prev_sidx = _node.prev.value.sidx if _node.prev else None
        # SolveLog.debug('create at {}'.format(_node))
        affected = _anneal.create(_node)
        self.default_local_lsq(_anneal, affected)
        create_dl = _anneal.get_score()
        store.append(_anneal)
        store_affected.append(affected)

        # ! destroy
        _anneal = self.anneal.clone()
        _node = _anneal.get_segment_node(sidx)
        # SolveLog.debug('destroy at {}'.format(_node))
        affected = _anneal.destroy(_node)
        self.default_local_lsq(_anneal, affected)
        destroy_dl = _anneal.get_score()
        store.append(_anneal)
        store_affected.append(affected)

        # ! choose
        dls = [init_dl, create_dl, destroy_dl]
        Ms = [len(_a.llmodel) for _a in store]
        scores = list(zip(range(3), dls, Ms))
        best = sorted(scores, key=lambda t: (t[1],t[2]))[0]
        bestidx = best[0]

        # ! bookkeeping
        result_string = ["no change", "create node", "delete node"]
        self._record.append(best)

        if bestidx == 0:
            pass # do nothing
            affected = store_affected[bestidx]
            _node = affected[0]
            _temp = self.fc * _node.value.temp
            self.new_priority_task(_node, 1. - _temp)
        elif bestidx == 1:
            # create
            self.deleted.add(sidx)
            affected = store_affected[bestidx]
            if self.use_priority:
                # update the prority of the first and last node
                for _node in affected:
                    _temp = 1.0
                    self.new_priority_task(_node, 1. - _temp)
            else:
                self.queue.append(affected[0].value.sidx)
                self.queue.append(affected[1].value.sidx)
        elif bestidx == 2:
            # destroy
            self.deleted.add(sidx)
            if prev_sidx != None: 
                self.deleted.add(prev_sidx)
            affected = store_affected[bestidx]
            if self.use_priority:
                for _node in affected:
                    _temp = 1.0
                    self.new_priority_task(_node, 1. - _temp)
            else:
                self.queue.append(affected[0].value.sidx)
        print("scores {} best (idx, DL, M) = {}".format(scores, best))
        print('result: {}'.format( result_string[bestidx]) )
        # SolveLog.info("scores {} best (idx, DL, M) = {}".format(scores, best))
        # SolveLog.info('result: {}'.format( result_string[bestidx]) )

        # ! return
        self.anneal = store[bestidx]
        # SolveLog.debug("new  model {} {}".format(id(self.anneal.llmodel), self.anneal.llmodel))
        return affected

    def random_queue_solve(self):
        def _next():
            if len(self.queue) == 0:
                return
            while True:
                sidx = self.queue.popleft()
                if sidx not in self.deleted:
                    return sidx
            return 
        self.sidx = _next()
        def _step():
            print("anneal at sidx = {}".format(self.sidx))
            try:
                self.create_and_destroy(self.sidx)
            except:
                print("failed at sidx =  {}".format(self.sidx))
                self.dump_state()
                raise
            print()
            self.sidx = _next()
        while self.sidx is not None:
            _step()


    def priority_solve(self, control={}):
        """
        # solve using a priority queue and "temperature" variable
        # our "tasks" will be segments about which to try and improve
    
        control parameters:
        Tc : stop temperature
        fc : fraction to decrease temperature after creation
        fd : fraction to decrease temperature after destroy
        fr : fraction to decrease temperature after remapping
        """
        self.greedy = control.get('greedy', self.greedy)
        self.Tc = control.get('Tc', 0.01)
        self.fc = control.get('fc', 0.3)
        self.fd = self.fc

        # setup priority queue
        self.use_priority = True
        for value in self.anneal.llmodel:
            value.temp = 1.0
        segment_node_data = [value for value in self.anneal.llmodel][:-1] 
        self.rng.shuffle(segment_node_data)
        q = [[1.0 - value.temp, value.sidx] for value in segment_node_data]
        pqueue = support.Modheapq.from_sorted(q)
        self.pqueue = pqueue
        count_iter = 0 
        while True:
            entry = self.pqueue.pop_entry()
            if entry is None:
                break # empty queue
            priority, sidx = entry
            T = 1.0 - priority
            # report every 100 iterations
            if count_iter % 100 == 0:
                print(f"iter {count_iter} T = {T} Tc = {self.Tc}", file=local_f)

            # exit condition
            if T <= self.Tc:
                break
            if sidx in self.deleted:
                continue
            node = self.anneal.get_segment_node(sidx)
            if node.next is None:
                continue

            try:
                self.random_transition(sidx)
            except:
                self.dump_state("partial")
                raise
            count_iter += 1
        print("End Priority Solve", file=local_f)

        
    def prep_queue(self, control={}):
        self.greedy = control.get('greedy', self.greedy)
        self.Tc = control.get('Tc', 0.01)
        self.fc = control.get('fc', 0.3)
        self.fd = self.fc
        segment_node_data = [value for value in self.anneal.llmodel][:-1] 
        self.rng.shuffle(segment_node_data)
        q = [[1.0 - value.temp, value.sidx] for value in segment_node_data]
        pqueue = support.Modheapq.from_sorted(q)
        self.pqueue = pqueue

    def cleanup(self, control={}):
        self.prep_queue()
        # a 'quick' final pass on the trajectory to see if any nodes can be deleted
        for node in nodeiter(self.anneal.llmodel):
            if node.next is None:
                continue
            self.destroy(node.value.sidx)


    # remap a point to a segment and recompute the local minimum
    def remap(self, pt_i, sidx):
        _anneal = self.anneal.clone()
        node = _anneal.sidxmap[sidx] 
        loss = _anneal.loss_function([node], n_local=0)
        pre_ssr= loss.terms[0]
        # ^ must compare with score prior to remapping
        _anneal.seglist.remap = {pt_i : node}

        # optimise single node position
        nm_options = {"fatol":0.001, "xatol":0.001}
        _anneal.local_lsq_optimise([node], n_local=1, options=nm_options, use_ordering=True)
        # self.default_local_lsq(_anneal, [node], n_local=0)

        aft_ssr = _anneal.get_loss().terms[0]
        if aft_ssr < pre_ssr:
            print("accepted remapping point {} segment {}".format(pt_i, sidx), file=local_f)
            self._record.append(('remap', pt_i, sidx, pre_ssr, aft_ssr))
            _anneal.seglist.remap = {}
            self.anneal = _anneal
        else:
            print('noremap', pt_i, sidx, pre_ssr, aft_ssr, file=local_f)
            # no change
        return self.anneal

    def remapping(self, control={}):
        """
        We can apply the remapping transition in various ways.
        Lets start small and just remap the nodes which are close to adjacent segments
        # i.e. have local coordiant (s < r, r > l-r)
        """
        lptr = self.anneal.get_current_model()
        lptr.sidx

        # split the points by closest segment, points will be time sorted (alternatively sort by local coord)
        M = self.anneal.M
        closest_segment = self.anneal.seglist.closest_segment
        # local_coord = self.anneal.seglist.get_local_coord()
        bysegment = {lptr.sidx[i] : [] for i in range(M-1)}
        for i, sidx in enumerate(closest_segment):
            bysegment[sidx].append(i)

        # it start with, create a single queue and then run down it
        queue = []
        for i in range(M-1):
            sidx = lptr.sidx[i]
            if i > 0:
                queue.append((bysegment[sidx][0], lptr.sidx[i-1]))
            if i < M-1 and lptr.sidx[i+1] != -1:
                queue.append((bysegment[sidx][-1], lptr.sidx[i+1]))

        #! tmp
        # queue = queue[:8]

        for pair in queue:
            newstate = self.remap(*pair)
            # !tmp
            # if newstate != self.anneal:
            #     self.anneal = newstate
            #     return


# annealing version of mdlgss from mdl module
# golden section search


def mdlgss(_data, r, bracket, rngstate=None):
    args = (_data, r)
    N = len(_data)
    if not rngstate:
        rngstate = np.random.RandomState(0)
    def loss(_M, _data, r):
        # setup and solve
        _guess = mdl.recursive_coarsen(_data, _M, parameter='M')
        SolveLog.info("new guess with M = {}".format(_M))
        anneal = Anneal()
        anneal.initialise(_guess, _daa)
        solver = Solver(anneal, rng=rngstate)
        # solver.linear_solve()
        SolveLog.info("start linear solve")
        solver.multiple_linear_solve(n=1)
        SolveLog.info("start annealing solve")
        solver.random_queue_solve()
        
        # 
        dl = solver.anneal.get_description_length()
        _model = solver.anneal.get_current_model()
        inter = {}
        inter["is_outlier"] = solver.anneal.get_outliers()
        return dl, _model, inter

    Intf = mdl.IntegerFunction(N, loss,  args=args)

    xM = mdl.convex_minimise(Intf, bracket)
    Intf.minx = xM
    return Intf

# ----------------------------------------------------------------

# convenience method for plotting from solver
def model_plot(ax, anneal, data, _config={"h_outlier":True}):
    # note specific plotting function
    is_outlier = anneal.get_outliers()
    mdl.plot_model_on_data(ax, anneal.get_current_model(), data, 
        intermediate={'is_outlier':is_outlier}, config=_config)
    


# convenince method for solving actual systems
def solve_exp(_data, M, r, default_loss_conf={}):
    def _print_break(message):
        print("----------------------------------------------------------", file=local_f)
        print(message, file=local_f)

    config = {"use_ordering": True}
    rngstate = np.random.RandomState(0)
    _guess = mdl.recursive_coarsen(_data, M, parameter='M')
    anneal = Anneal(r)
    anneal.default_loss_conf = default_loss_conf
    anneal.initialise(_guess, _data)
    solver = Solver(anneal, rng=rngstate, config=config)

    _print_break("initial solve")
    solver.multiple_linear_solve(n=1, n_local=2)

    _print_break("priority solve")
    pqsolve = {
        'greedy' : False, 
        'Tc' : 0.01, 
        'fc' : 0.5
        }
    solver.priority_solve(control=pqsolve)

    _print_break("greedy solve")
    greedy = {
        'greedy' : True, 
        'Tc' : 0.01, 
        'fc' : 0.2
    }
    solver.priority_solve(control=greedy)

    _print_break("remapping")
    solver.remapping()

    _print_break("cleanup")
    solver.cleanup()
    return solver
    

