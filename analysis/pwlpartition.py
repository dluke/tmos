
import sys, os
import time
join = os.path.join
import time
from copy import copy, deepcopy
import itertools
import numpy as np
import pickle
norm = np.linalg.norm
import scipy
import scipy.optimize
import scipy.stats

from pili import support
import matplotlib.pyplot as plt



import mdl

from skimage.restoration import estimate_sigma, denoise_wavelet

# array/mapping implementation for a piecewise linear solver
# use mdl.LPtrack

hide_output = False
local_f = sys.stdout

def local_print(*args):
    if hide_output:
        return
    print(*args, file=local_f)

# ----------------------------------------------------------------
# start with the wavelet initial condition

def estimate_error(x, y):
    sigma_est_x = estimate_sigma(x)
    sigma_est_y = estimate_sigma(y)
    sigma_est = np.mean([sigma_est_x, sigma_est_y])
    return sigma_est

scikit_config = {"wavelet":'db1', 'method':'VisuShrink', "mode":'soft', "rescale_sigma":False}
def vary_sigma(x, y, sigma, config=scikit_config):
    x_denoise = denoise_wavelet(x, sigma=sigma, **config)
    y_denoise = denoise_wavelet(y, sigma=sigma, **config)
    denoised = np.stack([x_denoise, y_denoise])
    return denoised


# some simple thresholds might improve or at least simplify the wavelet model
def contract(denoised, threshold=1e-4):
    # print('contraction threshold', threshold)
    N = denoised[0].size
    xyt = denoised.T
    #  the denoised trajectory has repeated values (VisuShrink) method
    diff = xyt[1:] - xyt[:-1]
    eq = norm(diff, axis=1) < threshold
    # seg_idx = np.concatenate([[0], np.argwhere(eq == False).ravel()+1, [N-1]])
    seg_idx = np.concatenate([[0], np.argwhere(eq == False).ravel()+1])
    if seg_idx[-1] != N-1:
        np.append(seg_idx, N-1)


    x, y = xyt[seg_idx].T
    dt = np.zeros(x.size, dtype=int)
    dt[1:] = np.diff(seg_idx)
    dt[1] += 1 # the 0th point represents a single timepoint
    return mdl.LPtrack(dt, x, y)

def model_from_denoised(denoised, sigma, angle_threshold=5):
    contracted = contract(denoised)

    # check the distances between nodes
    wavemodel = mdl.recursive_coarsen(contracted, sigma)

    def contract_by_angle(wavemodel, angle_threshold=5):
        # threshold in degrees
        theta = wavemodel.get_angle() * 180/np.pi
        list_theta = theta.tolist()[::-1]

        keepidx = [0]
        i = 1
        _curr = 0
        dtlist = [0]
        while list_theta:
            _curr += list_theta.pop()
            if abs(_curr) > angle_threshold:
                _curr = 0
                keepidx.append(i)
            i += 1
        keepidx.append(len(wavemodel)-1)
        n= len(keepidx)
        dt = np.zeros(n, dtype=int)
        dt[1:] = np.diff(np.cumsum(wavemodel.dt)[keepidx])
        model = mdl.LPtrack(dt, wavemodel.x[keepidx], wavemodel.y[keepidx])
        return model

    bmodel = contract_by_angle(wavemodel, angle_threshold)
    return bmodel

# alt transform is polar coordinate transform
# x.shape is (n, 2)
def alt_transform(x):
    # convert  [x0, y0, x1, y1, ...] -> [x0, y0, angle_1, length_1, ... ]
    alt = np.zeros(x.size)
    vector = x[1:] - x[:-1]
    length = norm(vector, axis=1)
    u = vector/length[:,np.newaxis]
    _x, _y = u.T
    theta = np.arctan2(_y, _x)
    alt[0], alt[1] = x[0]
    alt[2::2] = theta
    alt[3::2] = length
    return alt

def inv_transform(x):
    # convert [x0, y0, angle_1, length_1, ... ] -> [x0, y0, x1, y1, ...]
    n = x.size//2
    inv = np.zeros((2,n))
    theta, length = x[2:].reshape((n-1,2)).T
    _x = np.empty(n)
    _y = np.empty(n)
    _x[0] = x[0]
    _y[0] = x[1]
    _x[1:] = length * np.cos(theta)
    _y[1:] = length * np.sin(theta)
    vector = np.stack([_x, _y])
    inv = np.cumsum(vector, axis=1)
    return inv
    

# ----------------------------------------------------------------
# 
class PartAnneal(object):
    # partition anneal

    def __init__(self, model, data, r=0.03, 
            use_bounds=False, use_alternative=True, inter_term=False, 
            use_probability_loss=False, loss_conf={}):
        self.data = data
        self.point_data = data.get_n2()
        self.model = deepcopy(model)
        self.model_time = self.model.get_time()
        self.N = len(self.data)

        # transforms used local optimiser
        self.set_alternative(use_alternative)

        # cache
        self.cache_residuals = np.empty(self.N)
        self.cache_local_coord = np.empty(self.N)

        self._cache_inter_residual = None
        self._cache_terms = None
        self._cache_pvalue = None

        self._cache_local_optimise_result = []
        
        # setup model limits
        x, y =  data.x, data.y
        pad = 0.01 * max(abs(x[-1]-x[0]), abs(y[-1]-y[0]))
        rxl = (x.min()-pad, x.max()+pad)
        ryl = (y.min()-pad, y.max()+pad)
        self.point_data_limits = (rxl, ryl)

        # config
        self.r = r
        self.sigma = None
        self.use_bounds = use_bounds
        self.loss_conf= loss_conf
        self.inter_term = inter_term
        self.use_probability_loss = use_probability_loss

        self.update_residuals()
        
    def set_alternative(self, use_alternative):
        self._use_alternative = use_alternative
        # self.x_transform = alt_transform if use_alternative else lambda x: x.ravel()
        # self.x_inverse = inv_transform if use_alternative else lambda x: x.reshape((-1,2)).T

    def update_model_time(self):
        self.model_time = self.model.get_time()

    def update_residuals(self):
        self.residual_at(0, self.model.M)
        self.inter_residual_at(0, self.model.M-1)

    def get_residuals(self):
        return self.cache_residuals

    def get_local_coord(self):
        return self.cache_local_coord

    def get_curve_coord(self):
        local_coord = self.get_local_coord()
        mapping = self.get_mapping()
        lengths = self.model.get_step_length()
        L = np.insert(np.cumsum(lengths[:-1]), 0, 0)
        return L[mapping] + local_coord

    def compute_continuity_at(self, si, sf):
        exsi, exsf = self.model.expand_index(si, sf)
        boundary_left_local_idx = 1 if si == 0 else 0
        boundary_right_local_idx = -1 if exsf == self.model.M else None

        time = self.model_time[si:exsf]
        ind = time[boundary_left_local_idx : boundary_right_local_idx]
        
        lengths = self.model.lengths_at(exsi, exsf)
        model_curve_coord = np.insert(np.cumsum(lengths[:boundary_right_local_idx]), 0, 0)

        local_coord = self.get_local_coord()
        curve_coord_plus = model_curve_coord[1:] + local_coord[ind]
        curve_coord_minus = model_curve_coord[:-1] + local_coord[ind-1]
        adjacent_curve_distance = curve_coord_plus - curve_coord_minus

        data = self.point_data
        adjacent_data_distance = norm(data[ind] - data[ind-1], axis=1)
        continuity = adjacent_curve_distance - adjacent_data_distance
        return continuity

    def continuity_term_at(self, si, sf):
        # compute continuity
        continuity = self.compute_continuity_at(si, sf)
        # ! warning: the 2*self.r offset is needed ~ and may not be robust ~
        # offset = 2*self.r
        offset = self.r
        continuity = np.clip(continuity - offset, 0, None)
        return np.sum(continuity**2)

    def get_mapping(self):
        time = self.model.get_time()
        mapping = np.empty(self.N,dtype=int)
        for i in range(self.model.M-1):
            mapping[time[i]:time[i+1]] = i
        return mapping

    # compute end point distances here or in inter_residual
    def residual_at(self, si, sf, inter_term=None):
        if inter_term is None: inter_term = self.inter_term
        # compute residuals  
        # todo remove this O(n) operation?
        self.update_model_time()
        tidx = self.model_time[si:sf+1]
        ab = self.model.get_ab_at(si,sf)
        distance_list = []
        coord_list = []
        for i in range(tidx.size-1):
            ti, tf = tidx[i], tidx[i+1]
            p = self.point_data[ti:tf]
            a, b = ab[i], ab[i+1]
            a, b = a[np.newaxis], b[np.newaxis]
            s, d = support.line_coord_dist(p, a, b)
            # * adjustment for endpoints
            # * be careful, description length uses residuals too!
            if inter_term and len(d) > 1:
                s[0] = 0
                d[0] = np.sqrt(np.sum((p[0] - a)**2)) / np.sqrt(2.)
            distance_list.append(d)
            coord_list.append(s)
        residual = np.concatenate(distance_list)
        # print('residual', residual)
        self.cache_residuals[tidx[0]:tidx[-1]] = residual
        scoord = np.concatenate(coord_list)
        self.cache_local_coord[tidx[0]:tidx[-1]] = scoord
        return residual

    # def inter_residual_at(self, si, sf):
    #     # return sf - si - 1 residuals
    #     tidx = self.model_time[si:sf+1]
    #     ab = self.model.get_ab_at(si,sf)
    #     distance_list = []
    #     for i in range(sf-si-1):
    #         tf = tidx[i+1]
    #         p = self.point_data[tf-1]
    #         a, b = ab[i+1], ab[i+2]
    #         a, b = a[np.newaxis], b[np.newaxis]
    #         s, d = support.line_coord_dist(p, a, b)
    #         distance_list.append(d)
    #     arr = np.array(distance_list)
    #     self._cache_inter_residual = arr
    #     return arr

    def inter_residual_at(self, si, sf):
        # return sf - si - 1 residuals
        tidx = self.model_time[si:sf+1]
        ab = self.model.get_ab_at(si,sf)
        distance_list = []
        for i in range(sf-si-1):
            tf = tidx[i+1]
            p = self.point_data[tf-1]
            b = ab[i+1]
            # s = 0
            d = np.sqrt(np.sum((p - b)**2)) / np.sqrt(2.)
            # * could divide by sqrt(2.) since error distribution is 2d here
            distance_list.append(d)
        arr = np.array(distance_list)
        self._cache_inter_residual = arr
        return arr
         

    def get_loss_function(self, si, sf):
        if self.use_probability_loss:
            return self.get_probability_loss(si, sf)
        else:
            return self.get_loss(si, sf)

    def eval_loss_function(self, si, sf):
        return self.get_loss_function(si, sf)(self.get_x0(si, sf))

    #todo: should we use the standard residuals or residuals without clipping?
    def get_probability_loss(self, si, sf):
        exsi, exsf = self.model.expand_index(si, sf)
        
        x_inverse = inv_transform if self._use_alternative else lambda x: x.reshape((-1,2)).T

        def loss(x):
            m = sf-si

            _x, _y = x_inverse(x)
            self.model.x[si:sf] = _x
            self.model.y[si:sf] = _y

            pvalues = self.estimate_segment_probability_at(exsi,exsf)
            if (np.sum(np.isnan(pvalues)) != 0):
                print(pvalues)
                sys.exit()

            log_likelihood = np.sum(np.log(pvalues))

            # print('pvalues', pvalues, log_likelihood)

            if log_likelihood == np.inf:
                print('ll', log_likelihood)
                print('pvalues', pvalues)
                sys.exit()

            s0 = -1 * log_likelihood
            terms = [s0]

            # print('s0', s0, 'pvalues', pvalues)

            # todo boundary term

            self._cache_terms = terms
            
            value = np.sum(terms)
            # print('log likelihood', log_likelihood, pvalues)
            return value

        return loss

    def get_loss(self, si, sf):
        # configuration
        config = self.loss_conf
        contour = config.get('contour_term', None)
        continuity = config.get('continuity_term', None)

        x_inverse = inv_transform if self._use_alternative else lambda x: x.reshape((-1,2)).T

        # 
        exsi, exsf = self.model.expand_index(si, sf)
        def loss(x):
            # todo expand index range for LSQ computation
            m = sf-si

            _x, _y = x_inverse(x)
            self.model.x[si:sf] = _x
            self.model.y[si:sf] = _y

            residual = self.residual_at(exsi, exsf)
            s0 = np.sum(residual**2)
            terms = [s0]
            
            if self.inter_term:
                endsf = min(exsf, self.model.M-1)
                inter_residual = self.inter_residual_at(exsi, endsf)
                terms.append(np.sum(inter_residual**2))

            if contour != None and contour > 0:
                lengths = self.model.lengths_at(exsi, exsf)
                s1 = contour * sum(lengths) / len(lengths)
                terms.append(s1)

            if continuity != None and continuity > 0:
                s2 = continuity * self.continuity_term_at(exsi, exsf)
                terms.append(s2)


            # fix endpoints
            boundary = 0
            if exsi == 0:
                boundary += (self.model.x[0] - self.data.x[0])**2 + (self.model.y[0] - self.data.y[0])**2
            if exsf == self.model.M:
                boundary += (self.model.x[-1] - self.data.x[-1])**2 + (self.model.y[-1] - self.data.y[-1])**2

            # local_print(sum(terms), 'terms', terms, boundary)

            self._cache_terms = terms
            return sum(terms) + boundary
        return loss

    def get_x0(self, si, sf):
        x0 = self.model.get_ab_at(si,sf-1)

        x_transform = alt_transform if self._use_alternative else lambda x: x.ravel()
        return x_transform(x0)

    def get_loss_terms(self, si=0, sf=None):
        if sf is None: sf = self.model.M
        self.get_loss_function(si, sf)(self.get_x0(si, sf))
        return self._cache_terms

    def get_total_loss(self):
        si, sf = 0, self.model.M
        return self.eval_loss_function(si, sf)


    def local_optimise_at(self, si, sf, options={}, optconf={}):
        if "method" not in optconf:
            optconf["method"] = "Nelder-Mead"

        #! use LSQ regardless use_probability_loss?
        # 
        # if self.use_probability_loss:
        if False:
            loss = self.get_loss_function(si, sf)
        else:
            # be careful of this default
            tol = self.sigma/5. if self.sigma else 0.01
            options.update({"fatol": tol, "xatol": tol})
            loss = self.get_loss(si, sf)


        x0 = self.get_x0(si, sf)

        # bounds
        bounds = None
        if self.use_bounds:
            # if self._use_alternative:
            #     rxl, ryl = self.point_data_limits
            #     lim = max(rxl[1]-rxl[0], ryl[1]-ryl[0])
            #     # bounds on length in polar coords
            #     n = len(x0)
            #     bounds = [None for _ in range(n)]
            #     bounds[0] = rxl
            #     bounds[1] = ryl
            #     bounds[2::2] = [(0, np.inf) for _ in range(n//2 -1)]
            #     bounds[3::2] = [(0, lim) for _ in range(n//2 -1)]
            # else:

            # (x,y) coords so use square bounds
            rxl, ryl = self.point_data_limits
            n = len(x0)
            bounds = [None for _ in range(n)]
            bounds[0::2] = [rxl for _ in range(n//2)]
            bounds[1::2] = [ryl for _ in range(n//2)]
            optconf["bounds"] = bounds
        #

        start = time.perf_counter()
        res = scipy.optimize.minimize(loss, x0, options=options, **optconf)
        end = time.perf_counter()
        
        # hijacking scipy.optimize.OptimizeResult object to record the execution time
        res.exec_time = end - start 

        self._cache_local_optimise_result.append(res)
        local_print(mdl.describe_result(res))

        return loss

    def delete_at(self, index):
        self.model.delete(index)
        self.update_residuals()

    def create_at(self, index):
        # compute new x, y positions
        # several ways to consider implementing this
        # lets try bisecting on the middle data point
        if index >= self.model.M-1:
            return False
        sdt = self.model.dt[index+1]

        # todo. bug: relaxing this constraint allows segments with 0 data to be generated at the start of the model
        if sdt < 2:
            return False

        dt1 = sdt//2

        # get the position of the middle data point
        time = self.model.get_time() # todo O(n)
        mid_data_index = time[index] + dt1
        x = self.data.x[mid_data_index]
        y = self.data.y[mid_data_index]

        self.model.create(index, x, y, dt1)

        self.update_residuals()
        return True


    def clone(self):
        return deepcopy(self)

    def get_outliers(self, r):
        # return self.get_residuals() > r
        residuals = self.residual_at(0, self.model.M, inter_term=False)
        return residuals > r

    def estimate_segment_probability(self):
        return self.estimate_segment_probability_at(0, self.model.M-1)

    def estimate_segment_probability_at(self, si, sf):
        residuals = self.get_residuals()
        exsf = min(sf+1, self.model.M-1)
        inter_residual = np.append(self.inter_residual_at(si, exsf), 0)

        self.update_model_time()
        tidx = self.model_time[si:sf+1]

        # chunk size
        n = tidx.size-1
        pvalue = np.empty(n)
        cdfv = np.empty(n)

        for i in range(n):
            ti, tf = tidx[i], tidx[i+1]
            n = tf - ti
            if n == 0:
                local_print('warning: segment with no data = 0')
                pvalue[i] = 1
                cdfv[i] = 0
            else:
                if self.inter_term:
                    res = np.append(residuals[ti:tf], inter_residual[i])
                else:
                    res = residuals[ti:tf]


                sample_residual = np.sum(res**2/self.sigma**2)
                chi = scipy.stats.chi2(res.size)
                cdf = chi.cdf(sample_residual)

                # print('residual i', i, 'n', res.size, 'ti tf', ti, tf, res, cdf)

                # avoid numerical issues with log(0)
                cdf = np.clip(cdf, 0, 1.0-1e-10)

                cdfv[i] = cdf
                pv = 1-cdf if cdf > 0.5 else cdf
                pvalue[i] = pv
        

        # print('cdfv', cdfv)

        self._cache_cdfv = cdfv

        return pvalue

    def get_segment_cdfv(self):
        # call immediately after estimate_segment_probability
        return self._cache_cdfv

    # * analysis helpers

    def count_segment_data(self, r):
        # for each segment in the model, count not only the number of associated data but a
        # all data within r of either end by looking at the adjacent segments
        time = self.model.get_time()
        M = self.model.M
        cts = []
        for i in range(M-1):
            # adjacent segments
            # adjacent_si = [_i in [i-1, i, i+1] if _i > 0 and _i < M-1]
            # primary count
            si, sf = i, i+1
            count = self.model.dt[sf]

            if i-1 > 0:
                si, sf = i-1, i
                m_pt = self.model.get_point(sf)[np.newaxis]
                p = self.point_data[time[si]:time[sf]]
                count += int(np.sum(norm(p - m_pt, axis=1) < r))

            if i+1 < M-1:
                si, sf = i+1, i+2
                m_pt = self.model.get_point(sf)[np.newaxis]
                p = self.point_data[time[si]:time[sf]]
                count += int(np.sum(norm(p - m_pt, axis=1) < r))
            
            cts.append(count)
        return np.array(cts)



class Output(object):

    def __init__(self, outputdir, index=None, freq=100, allow_write=True):
        self.outputdir = outputdir
        self.freq = freq
        self.form = "iteration_{:08d}"
        self.allow_write = allow_write
        self.index = index
        
        # 
        if not os.path.exists(outputdir):
            print('creating directory', outputdir)
            os.makedirs(outputdir)

    def ready(self, count_iter):
        return count_iter % self.freq == 0

    def report(self, count_iter, solver):
        if self.index is None:
            print(f"solver iteration {count_iter:08d}")
        else:
            print(f"solver {self.index} iteration {count_iter:08d}")

    def dump(self, count_iter, solver):
        if self.allow_write:
            path = join(self.outputdir, self.form.format(count_iter))
            solver.dump_state(path)
        

# direction to move points on the partition
FORWARD = 0
REVERSE = 1

class Solver(object):

    DESTROY = 0
    CREATE = 1
    PERCOLATE = 2

    mv_name = ['destroy', 'create', 'percolate']

    def __init__(self, partition, rng=None, r=0.03, sigma=None, 
            use_modified_dl=False, use_binary_percolate=False, 
            min_constraint=0, use_description_length=True, optconf={}):
        self.rng = rng if rng!=None else np.random.RandomState(0)
        self.partition = partition

        # config
        self.sigma = sigma
        self.partition.sigma = sigma
        self.r = r
        self.partition.r = r
        self.use_description_length = use_description_length
        self.use_modified_dl = use_modified_dl

        # history
        self._score_history = []
        self._loss_history = []
        self._mv_history = [] # record the monte carlo move selected at this step, i.e.

        # priority queue solve
        self.optconf = optconf
        self.greedy = None
        self.T = 0
        self.min_constraint = min_constraint
        self.use_binary_percolate = use_binary_percolate


    def update_history(self):
        self._score_history.append(self.get_score())
        self._loss_history.append(self.partition.get_loss_terms())

    def get_history(self):
        return np.array(self._score_history), np.array(self._loss_history)

    def get_outliers(self):
        return self.partition.get_outliers(self.r)

    def get_score(self, partition=None):
        if not partition: partition = self.partition

        if self.use_description_length:
           return self.get_description_length(partition)
        else:
           return self.get_chi2_score(partition)

    def get_chi2_score(self, partition=None):
        if not partition: partition = self.partition
        # todo remove O(M) dependence

        pvalue = partition.estimate_segment_probability()
        cdfv = partition.get_segment_cdfv()

        outlier = cdfv > self.p_threshold
        k = int(np.sum(outlier))

        loss = partition.get_total_loss()

        M = partition.model.M
        score = (k, M, loss)

        # local_print("score", score)
        return score


    def get_description_length(self, partition=None):
        if not partition: partition = self.partition
        outliers = partition.get_outliers(self.r) 
        n = int(np.sum(outliers))
        M = partition.model.M
        if self.use_modified_dl:
            mo = int(np.sum(self.get_model_residual() > 2 * self.r))
            return (n + mo + M, mo, M)
        else:
            return (n + M, M)

    def get_model_residual(self):
        # compute the min distance from the model points to the data
        time = self.partition.model_time
        model = self.partition.model
        residuals = []
        for index in range(model.M):
            si, sf = max(index-1, 0), min(index+1, model.M-1)
            ti, tf = time[si], time[sf]
            data = self.partition.point_data[ti:tf]
            pt = model.get_point(index)

            # (n, 2) - (1,2)
            vector = data - pt[np.newaxis,:]
            distance = norm(vector, axis=1).min()
            residuals.append(distance)

        return np.array(residuals)

    def estimate_segment_probability(self):
        return self.partition.estimate_segment_probability()



    # ----------------------------------------------------------------
    # probability solve


    def chi2_optimise_chunk(self, si, sf, maxiter=10):
        # if we minimise loss then we want to maximise gain, right?
        exsi, exsf = self.partition.model.expand_index(si, sf)

        loss = self.partition.get_probability_loss(si, sf)

        # 
        for i in range(maxiter):
            print('iter', i)

            # choose a random node
            index = self.rng.randint(si, sf)

            self.percolate_at(index)

            # # choose a random direction
            # dirc = self.rng.randint(2)
            # pre = loss(self.partition)
            # print("random index", index, "direction", dirc)

            # print("before", pre, self.partition.estimate_segment_probability_at(exsi,exsf)

            # clone = self.percolate_one(index, direction=dirc, dn=1)
            # aft = loss(clone) 
            
            # print("after", aft, clone.estimate_segment_probability_at(exsi:exsf))

            # if aft < pre:
            #     print('accept', aft, pre)
            #     self.partition = clone

            # print()



    def chi2_solve(self, control, output=None):
        # self.use_description_length = False

        self.p_threshold = control.get('p_threshold', 0.98)
        maxiter = control.get('maxiter', np.inf)

        tolerance = control.get('tolerance', 1e-6)
        window = control.get('window', 100)

        p_percolate = 0.5
        # p_create = 0.5

        self.count_iter = 0
        if output and output.ready(self.count_iter):
            output.report(self.count_iter, self)
            output.dump(self.count_iter, self)

        while True:
            M = self.partition.model.M
            
            # 
            index = self.rng.randint(0,M-1)

            pvalue = self.estimate_segment_probability()
            cdfv = self.partition.get_segment_cdfv()
            k = int(np.sum(pvalue > self.p_threshold))

            local_print(self.count_iter, 'index', index, 'pvalue', pvalue[index])

            # the monte carlo move
            mv = None

            if self.rng.uniform() > p_percolate:
                mv = Solver.PERCOLATE
                self.percolate_at(index)
            
            else:
                if cdfv[index] < self.p_threshold:
                    mv = Solver.DESTROY
                    self.destroy_at(index, percolate=True)
    
                else:
                    mv = Solver.CREATE
                    self.chi2_create_at(index)
                    # 


            self._mv_history.append(mv)
            self.count_iter += 1
            self.update_history()

            if output and output.ready(self.count_iter):
                output.report(self.count_iter, self)
                output.dump(self.count_iter, self)

            if tolerance != None and k == 0 and self.check_converged(window=window, tol=tolerance):
                local_print("converged at tolerance after {} iterations".format(self.count_iter))
                break

            if self.count_iter >= maxiter:
                local_print("exit at max iterations ", self.count_iter)
                break
        
    def chi2_create_at(self, index, local=2, recurr=False):
        clone = self.partition.clone()
        created = clone.create_at(index)
        if not created: 
            local_print("did not create at", index)
            return False

        si, sf = max(0, index-local), min(clone.model.M, index+local)
        clone.local_optimise_at(si, sf, optconf=self.optconf)

        # todo: anneal
        # clone = self.percolate_at(index, partition=clone)
        clone = self.binary_percolate_at(index, partition=clone)
            
        # in chi2 solve we always accept because we require k = 0
        accepted = True
        self.partition = clone
        local_print('chi2 create at ', index)

        if recurr:
            # if there are still segments with large quantiles, split them as well
            clone.estimate_segment_probability_at(index, index+2)
            cdfv = clone.get_segment_cdfv()
            print('after create cdfv', cdfv)
            where = cdfv > self.p_threshold
            if any(where):
                recurr_index = index + where.tolist().index(True)
                print('create recursively at index', recurr_index)
                return self.chi2_create_at(recurr_index)

        return accepted


    class Thermostat(object):

        def __init__(self, solver):
            # destroy / create percolate
            self.solver = solver
            self.tlist = [[1., 1., 1.] for _ in range(self.model.M)]
            self.tenum = [Solver.DESTROY, Solver.CREATE, Solver.PERCOLATE]
            q = [[0, sidx] for sidx in self.model.sidx]
            self.pqueue = support.Modheapq.from_sorted(q)

        @property
        def model(self):
            return self.solver.partition.model
        
        def pop_entry(self):
            entry = self.pqueue.pop_entry()
            if entry is None:
                return None
            priority, task = entry
            index = np.argwhere(self.model.sidx == task)[0][0] # O(n)
            # print('get index', index, 'task', task, self.model.sidx)
            return priority, task, index
    
        def get_temperature_at(self, index):
            return sum(self.tlist[index])/3

        def get_next_move_at(self, index):
            node = self.tlist[index]
            t = max(node)
            choose = [i for i in range(len(node)) if node[i] == t]
            mvid = self.solver.rng.choice(choose)
            return t, self.tenum[mvid]

        def diffuse_temperature_at(self, index):
            local = [index, min(index+1, len(self.tlist)-1)]
            temps = [self.tlist[index] for index in local]
            print(local, temps)
            def mean(lst): return sum(lst)/len(lst)
            return list(map(mean, zip(*temps)))

        def create_at(self, index, excite=True):
            if excite:
                node = [1, 1.0, 1]
            else:
                node = self.diffuse_temperature_at(index)
            self.tlist.insert(index+1, node)
            priority = 1 - self.get_temperature_at(index)
            self.pqueue.add_entry([priority, self.model.sidx[index+1]])

            # add back the index that was popped
            if excite:
                self.tlist[index] = [1, 1, 1]
            priority = 1 - self.get_temperature_at(index)
            self.pqueue.add_entry([priority, self.model.sidx[index]])

            if excite:
                # update the right side index
                if index + 2 < self.model.M-1:
                    self.tlist[index+2] = [1, 1, 1]
                    priority = 1 - self.get_temperature_at(index+2)
                    self.pqueue.update_task(self.model.sidx[index+2], priority)


        def destroy_at(self, index):
            del self.tlist[index]

            adjacent = [_i for _i in [index-1, index] if _i > 0 and _i < self.model.M-1]
            for adj_index in adjacent:
                self.tlist[adj_index] = [1, 1, 1]
                priority = 1 - self.get_temperature_at(adj_index)
                self.pqueue.update_task(self.model.sidx[adj_index], priority)
            
        def update(self, index, mvid, k):
            self.tlist[index][mvid] *= k
            priority = 1 - self.get_temperature_at(index)
            # because we used pop_entry, so at this point we need to add the task back the queue
            self.pqueue.add_entry([priority, self.model.sidx[index]])

    def get_thermostat(self):
        return Solver.Thermostat(self)

    # ----------------------------------------------------------------
    # priority queue implementation
    def priority_solve(self, control, output=None):
        """ 
        This time will be slightly different than annealing.py
        (i) track a 'temperature' for each operation and when a node is popped pick the transition with largest t
        not random enough?

        """ 
                
        self.update_history()
        thermostat = Solver.Thermostat(self)


        # exit conditions
        maxiter = control.get('maxiter', np.inf)
        tolerance = control.get('tolerance', None)
        window = control.get('window', 100)
        t_end = control.get('t_end', 1/8. + 1/16.)

        # 
        self.greedy = control.get("greedy", True)
        k_destroy = control.get('k_destroy', 0.5)
        k_create = control.get('k_destroy', 0.5)
        k_percolate = control.get('k_destroy', 0.5)
        k_list = [k_destroy, k_create, k_percolate]

        self.count_iter = 0
        if output and output.ready(self.count_iter):
            output.report(self.count_iter, self)
            output.dump(self.count_iter, self)

        while True:

            # print('pqueue', thermostat.pqueue)
            # print('sidx', self.partition.model.sidx)
            # print()

            tup = thermostat.pop_entry()
            if tup is None:
                break # empty queue
            priority, task, index = tup
            T = 1 - priority
            self.T = T

            mv_temp, mvid = thermostat.get_next_move_at(index)
            local_print('iter', self.count_iter, 'T', T, 'sidx', task, 'index', index, 'move', Solver.mv_name[mvid], 'M', self.partition.model.M)
            

            accepted = False
            if mvid == Solver.DESTROY:
                accepted = self.destroy_at(index)
                if accepted:
                    thermostat.destroy_at(index)
                    
            elif mvid == Solver.CREATE:
                accepted = self.create_at(index)
                if accepted:
                    is_better = self.get_score()[0] < self.get_last_score()[0]  
                    thermostat.create_at(index, excite=is_better)
                    
            elif mvid == Solver.PERCOLATE:
                self.percolate_at(index, partition=self.partition)

            if not accepted:
                thermostat.update(index, mvid, k_list[mvid])

            self.count_iter += 1
            self.update_history()

            if output and output.ready(self.count_iter):
                output.report(self.count_iter, self)
                output.dump(self.count_iter, self)

            pair = (T, self.count_iter)
            if T < t_end:
                local_print("exit at critical temperature", pair)
                break

            if tolerance != None and self.check_converged(window=window, tol=tolerance):
                local_print("converged at tolerance after {} iterations".format(self.count_iter))
                break

            if self.count_iter >= maxiter:
                local_print("exit at max iterations ", pair)
                break
        
        return (T, self.count_iter)

    def get_last_score(self):
        return self._score_history[-1]

    def random_accept_score(self, pre_score, aft_score, T):
        # if description length, the first item is DL
        # if chi2, the items are (k, M, P(M))

        # self.accept_worse = dl2 < dl1
        if self.greedy:
            return pre_score < aft_score
        else:
            dl1 = pre_score[0]
            dl2 = aft_score[0]

            # !shift by +1 after accept
            s = np.exp(-(dl2 - dl1 + 1)/T) 

            u = self.rng.uniform()
            return s >= u

    def check_converged(self, window=100, tol=1e-6):
        if self.count_iter < window:
            return False
        loss = list(map(sum, self._loss_history[-window:]))
        dl = [t[0] for t in self._score_history[-window:]]

        return (np.max(np.abs(np.gradient(dl))) < tol
            and np.max(np.abs(np.gradient(loss))) < tol)

    # ----------------------------------------------------------------
    # operations which run the length of the model

    def linear_solve(self):
        # shuffle the nodes
        index_list = list(range(0, self.partition.model.M))
        self.rng.shuffle(index_list)
        for index in index_list:
            si, sf = self.partition.model.expand_index(index, index+1, local=3)
            self.partition.local_optimise_at(si, sf, optconf=self.optconf)

    def create(self):
        # attempt destroy on each segment in turn
        index = 0
        while index < self.partition.model.M-1:
            accepted = self.create_at(index, percolate=True)
            if not accepted:
                index += 1

    def cleanup(self):
        # attempt destroy on each segment in turn
        index = 1
        while index < self.partition.model.M:
            accepted = self.destroy_at(index, percolate=False)
            if not accepted:
                index += 1

    def percolate(self):
        # pick a random node, try to shift the mapping in both directions
        index_list = list(range(1, self.partition.model.M-1))

        self.rng.shuffle(index_list)
        for index in index_list:
            self.percolate_at(index)
            # self.partition = self.percolate_loop_at(index, FORWARD, self.partition)
            # self.partition = self.percolate_loop_at(index, REVERSE, self.partition)

    def percolate_at(self, index, partition=None, local=2):
        partition = self.partition if not partition else partition

        if self.use_binary_percolate:
            partition = self.binary_percolate_at(index, partition=partition, local=local)
            return partition

        else:
            partition = self.percolate_loop_at(index, FORWARD, partition)
            partition = self.percolate_loop_at(index, REVERSE, partition)
            return partition

    # ~
    # def random_accept_lsq(self, pre_lsq, aft_lsq, T):
    #     if self.greedy:
    #         return aft_lsq < pre_lsq
    #     else:
    #         s = np.exp(-(aft_lsq - pre_lsq)/T) 
    #         u = self.rng.uniform()
    #         return s >= u



    def anneal_around(self, index, local=2):
        partition = self.partition
        # use a thermostat to search for the  

        # construct a discrete probability distribution around index
        local_proportion = [0.1, 0.4, 1.0, 0.4, 0.1]
        local_pdf = np.array(local_proportion)/np.sum(local_proportion)
        cumulative = np.cumsum(local_pdf) 
        def sample():
            u = self.rng.uniform()
            local_index = np.searchsorted(cumulative, u)
            return index + (local_index - 2)

        maxiter = 10
        for i in range(maxiter):
            sample_index = sample()
            self.percolate_at(index)
            # accepted = self.binary_percolate_at(sample_index, local=2)
            # instead of this binary search idea, we can stick with percolating except use annealing

    def local_optimise_at(self, index, local=2):
        si, sf = self.partition.model.expand_index(index, index+1, local=local)
        self.partition.local_optimise_at(si, sf)
        
    def binary_percolate_at(self, index, partition=None, local=2):
        partition = self.partition if not partition else partition

        if index < 1 or index >= partition.model.M-1:
            return partition

        si, sf = partition.model.expand_index(index, index+1, local=2)

        # left bound, position, right bound
        px = partition.model.dt[index] 
        L, R = self.min_constraint, px + partition.model.dt[index+1] - self.min_constraint
        # print('bounds', L, px, R)

        pre_lsq = partition.eval_loss_function(si, sf)
        atscore = {px: pre_lsq}

        it = REVERSE
        while (R-L) > 0:
            at = L if it == REVERSE else R
            dn = at - px
            if at not in atscore:
                clone = partition.clone()
                clone.model.move_data_at(index, dn)
                clone.local_optimise_at(si, sf, optconf=self.optconf)
                aft_lsq = clone.eval_loss_function(si, sf) if at not in atscore else atscore.get(at)
                atscore[at] = aft_lsq
            else:
                aft_lsq = atscore[at]

            # print('at ', at, 'scores', pre_lsq, aft_lsq)

            if aft_lsq < pre_lsq:
                # accept
                partition = clone
                pre_lsq = aft_lsq
                px = at
                local_print('accepted move data dn = ', dn, aft_lsq, pre_lsq)
            else:
                # reject
                d = abs(dn)
                d = d//2 if d > 3 else d - 1 
                if it == REVERSE:
                    L = min(L + d, px) if d > 0 else px
                elif it == FORWARD:
                    R = max(R - d, px) if d > 0 else px
                
                it = FORWARD if it == REVERSE else REVERSE
                
            # print('move data at index', index, dn, partition.model.dt)
            # print('R, px, L', L, px, R)
        
        return partition

    # todo: may need a better local solve method than this ...
    def percolate_loop_at(self, index, direction, partition):
        if index < 1 or index >= partition.model.M-1:
            return partition
        loop = True
        while loop:
            si, sf = partition.model.expand_index(index, index+1, local=2)
            pre_lsq = partition.eval_loss_function(si, sf)
            
            # print('loss', partition.get_loss_terms(si, sf))

            clone = self.percolate_one(index, direction, partition)
            aft_lsq = clone.eval_loss_function(si, sf)

            # print('loss', clone.get_loss_terms(si, sf))
            # print('direction', direction, 'loss',pre_lsq, aft_lsq, aft_lsq < pre_lsq) 
            # print()

            if aft_lsq < pre_lsq:
                partition = clone
            else:
                loop = False
        return partition

    def percolate_one(self, index, direction=FORWARD, partition=None, dn=1):
        partition = self.partition if not partition else partition

        si, sf = partition.model.expand_index(index, index+1, local=2)
        clone = partition.clone()

        # todo: we want our initial guess to satisfy this constraint?
        min_constraint = self.min_constraint

        if direction == FORWARD:
            # percolate right 
            if clone.model.dt[index] < dn + min_constraint:
                return clone
            clone.model.dt[index] -= dn
            clone.model.dt[index+1] += dn
        elif direction == REVERSE:
            # percolate left
            if clone.model.dt[index+1] < dn + min_constraint:
                return clone
            clone.model.dt[index+1] -= dn
            clone.model.dt[index] += dn

        # 
        loss = clone.local_optimise_at(si, sf, optconf=self.optconf)
        return clone

    def destroy_at(self, index, local=2, percolate=True):
        # !why not allow deleting the first node?
        if index < 1: 
            return False

        pre_score = self.get_score(self.partition)
        clone = self.partition.clone()
        clone.delete_at(index)

        si, sf = max(0, index-local), min(clone.model.M, index+local)
        clone.local_optimise_at(si, sf, optconf=self.optconf)
        immediate_aft_score = self.get_score(clone)
        if immediate_aft_score < pre_score:
            # already good enough
            local_print('destroy without percolating', pre_score, immediate_aft_score)
            self.partition = clone
            return True

        if percolate:
            clone = self.percolate_at(index, partition=clone)
            # M =  clone.model.M
            # adjacent = [index, max(index-1,0), index, min(index+1,M-1)]
            # for idx in adjacent:
            #     print('percolate after destroy at idx = ', idx)
            #     clone = self.percolate_at(idx, partition=clone)
            # clone = self.binary_percolate_at(index, partition=clone)
            
        aft_score = self.get_score(clone)
        if aft_score < pre_score:
            accepted = True
            self.partition = clone
        else:
            accepted = False
        local_print('destroy', accepted, pre_score, aft_score)
        return accepted

    # def accept_condition(self, pre_score, aft_score):
    #     pre_score = partition.get_score()
    #     aft_score = clone.partition.get_score()
    #     return aft_score < pre_score

    def create_at(self, index, local=2, percolate=True):
        pre_score = self.get_score(self.partition)
        clone = self.partition.clone()
        created = clone.create_at(index)
        if not created: 
            local_print("did not create at", index)
            return False

        si, sf = max(0, index-local), min(clone.model.M, index+local)
        clone.local_optimise_at(si, sf, optconf=self.optconf)
        if percolate:
            # clone = self.percolate_at(index, partition=clone)
            clone = self.binary_percolate_at(index, partition=clone)
            

        aft_score = self.get_score(clone)
        if self.random_accept_score(aft_score, pre_score, self.T):
            if aft_score > pre_score:
                local_print("random accept create with T = ", self.T)

            accepted = True
            self.partition = clone
        else:
            accepted = False
        local_print('create', accepted, pre_score, aft_score)
        return accepted


    def dump_state(self, path):
        if not path.endswith('.pkl'):
            path += '.pkl'
        local_print('writing to ', path)
        with open(path, 'wb') as f:
            pickle.dump(self, f)

    def clone(self):
        return deepcopy(self)

    @classmethod
    def load_state(cls, path):
        if not path.endswith('.pkl'):
            path += '.pkl'
        local_print("load state at ", path)
        with open(path, 'rb') as f:
            new = pickle.load(f)
        return new

    def get_iter_count(self):
        return getattr(self, 'count_iter', None)

    def get_model(self):
        return self.partition.model

    def get_data(self):
        return self.partition.data



# ------------------------------------------------------------------------------
# solver scheme
# an oversimplified percolate-cleanup scheme 


import scipy.stats
import scipy.interpolate

def estimate_r(sigma, p=0.98):
    _p = 1 - (1-p)/2
    normal = scipy.stats.norm(scale=sigma)
    return normal.ppf(_p)

# numerical
# def estimate_r(sigma, p=0.98):
#     normal = scipy.stats.norm(scale=sigma)
#     basis = np.linspace(-3*sigma, 3*sigma, 1000)

#     def p_contained(r):
#         # probability of single normal distributed sample being contained in [-r, r]
#         return 1 - 2 * normal.cdf(-r)

#     # invert numerically
#     rbase = np.linspace(0, 5*sigma, 1000)
#     evf = np.array([p_contained(r) for r in rbase])

#     f = scipy.interpolate.interp1d(evf, rbase)

#     r = float(f(p))



def sim_wavemodel(data, sigma, angle_threshold=0):
    denoised = vary_sigma(data.x, data.y, sigma)
    wavemodel = model_from_denoised(denoised, sigma=sigma, angle_threshold=angle_threshold)
    return wavemodel


def initial_guess(x, y, params={"coarsen_sigma" : None}):
    sig = estimate_error(x ,y)

    coarsen_sigma = params.get("coarsen_sigma", None)
    if coarsen_sigma is None:
        coarsen_sigma = sig

    angle_threshold = params.get("angle_threshold", 5)

    print(f'using sigma = {sig}')
    p_threshold = params.get("p_threshold", 0.98)
    threshold_r = estimate_r(sig, p=p_threshold)
    denoised = vary_sigma(x, y, sig)

    print(f'coarsening with d = {coarsen_sigma}')
    wavemodel = model_from_denoised(denoised, sigma=coarsen_sigma, angle_threshold=angle_threshold)
    meta = {'sigma': sig, 'r': threshold_r}
    lptr = mdl.LPtrack(None, x, y)
    return wavemodel, lptr, meta

def wavelet_guess(track, N=None, head=True):
    xy = track.get_head2d()[:N] if head else track.get_trail2d()[:N]
    x, y = xy.T
    return initial_guess(x, y)

default_loss_conf = {'contour_term': 0.01}
def percolate_cleanup(wavemodel, lptr, loss_conf=default_loss_conf, dump=False, outputdir='', n_iterate=5):

    partition = PartAnneal(wavemodel, lptr, loss_conf=loss_conf)
    rng = np.random.RandomState(0)
    solver = Solver(partition, rng=rng)

    def pathgen():
        form = 'iteration_{:04d}'
        i = 0
        while True:
            yield form.format(i)
            i += 1
    path = pathgen()

    with support.Timer() as t:
        solver.linear_solve()
        solver.cleanup()
        def iterate():
            solver.create()
            solver.percolate()
            solver.cleanup()
            solver.percolate()
            if dump:
                solver.dump_state(join(outputdir, next(path)))

        for i in range(n_iterate):
            iterate()

    
    return solver


# -------------------------------------------------------------------------
# plotting

def plot_wavemodel_on_data(ax, wavemodel, data, config={}):
    match_data = config.get("match_data", False)
    color = itertools.cycle(['#FEC216', '#F85D66', '#75E76E'])
    ptlkw = {"linestyle":'--', "marker":"o", "alpha":0.4}
    

    if match_data:
        # use the color cycle and the time data to plot chunks of points and lines  with matching colors
        time = wavemodel.get_time()
        # x, y = data['x'], data['y']
        x, y = data.x, data.y
        xwave, ywave = wavemodel.x, wavemodel.y
        for i in range(len(wavemodel)-1):
            _col = next(color)
            ti, tf = time[i], time[i+1]
            ax.plot(x[ti:tf], y[ti:tf], color=_col, **ptlkw) 
            ax.plot([xwave[i], xwave[i+1]], [ywave[i], ywave[i+1]], color=_col, lw=2, alpha=0.6)
        ax.plot(xwave, ywave, linestyle='none', marker='D', alpha=0.6)

    else:
        ax.plot(data['x'],  data['y'], label="data", **ptlkw)
        xwave, ywave = wavemodel.x, wavemodel.y
        _label = config.get("denoised_label", "my wavelet")
        ax.plot(xwave, ywave, lw=4, marker='D', label=_label)
        ax.plot(xwave, ywave, linestyle='none', marker='D', markersize=8, c='#FC7B33')
    ax.set_aspect("equal")


def model_plot(solver, data, fs=20, config={}):
    fig, ax = plt.subplots(figsize=(fs,fs))
    is_outlier = solver.get_outliers()
    _config = {'h_outlier': True, 'r' : solver.r}
    _config.update(config)
    mdl.plot_model_on_data(ax, solver.partition.model, data, intermediate={'is_outlier':is_outlier}, config=_config)
    return fig, ax

def simple_model_plot(ax, model, data=None, pwl=None):
    defcolors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    c1, c2, c3 = defcolors[:3] 
    
    truth_style = {"linestyle": '--', 'lw':2, 'alpha':0.5, 'label':'truth', 'color':c1}
    model_style = {"linestyle": '-', 'lw':4, 'alpha':0.5, 'label':'model', 'color':c2}
    ptlkw = {"linestyle":'--', "marker":"o", "alpha":0.3, 'color':c3}

    if pwl:
        ax.plot(pwl.x, pwl.y, **truth_style)

    ax.plot(model.x, model.y, **model_style)

    if data:
        ax.plot(data.x, data.y, label='data', **ptlkw)

    ax.set_aspect('equal')
    ax.legend(fontsize=20, loc=(1.04, 0))

def plot_solver_convergence(solver, labels=None):
    score, loss = solver.get_history()
    fig, axes = plt.subplots(1, 2, figsize=(12,4))
    score_label = "DL" if solver.use_description_length else "score"
    ax = axes[0]
    ax.plot(score)
    if labels:
        ax.legend(labels=labels)
    ax.set_ylabel(score_label)
    ax = axes[1]
    ax.plot(loss.sum(axis=1))
    ax.set_ylabel('loss')
