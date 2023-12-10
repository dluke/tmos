
import numpy as np
import scipy.stats
import itertools
pi = np.pi
norm = np.linalg.norm

import mdl

class Gen(object):
    
    def sample(self, N=1):
        pass

class Constant(Gen):

    def __init__(self, v):
        self.value = v

    def sample(self, N=1):
        if N == 1:
            return self.value
        else:
            return np.full(N, self.value)

class Cycler(Gen):

    def __init__(self, _seq):
        self._seq = _seq
        self.seq = itertools.cycle(_seq)

    def sample(self, N=1):
        if N == 1:
            return next(self.seq)
        else:
            a = np.array([next(self.seq) for _ in range(N)])
            return a



class Uniform(Gen):

    def __init__(self, a, b):
        self.a = a
        self.b = b

    def sample(self, N):
        a, b = self.a, self.b
        return np.random.uniform(a,b,size=N)

# ...
class Exponential(Gen):

    def __init__(self):
        pass

class Normal(Gen):

    def __init__(self, loc=0, scale=1):
        self.loc = loc
        self.scale = scale
        self.normal = scipy.stats.norm(loc, scale)

    def pdf(self, x):
        return self.normal.pdf(x)

    def sample(self,N=1):
        return np.random.normal(self.loc, self.scale, size=N)

# create an offset and mirrored normal distribution to artificially generate trajectories
# where adjacent segments always have a significant angle

class MirrorNormal(Gen):

    def __init__(self, loc=0, scale=1):
        self.loc = loc
        self.scale = scale

        self.right_normal = scipy.stats.norm(loc, scale)
        self.left_normal = scipy.stats.norm(-loc, scale)

    def pdf(self, x):
        return (self.left_normal.pdf(x) + self.right_normal.pdf(x))/2

    def sample(self,N=1):
        x = np.random.normal(self.loc, self.scale, size=N)
        flipmask = np.random.random(N) < 0.5
        x[flipmask] = -x[flipmask]
        return x

class AlternatingMirrorNormal(MirrorNormal):

    # def __init__(self, loc=0, scale=1):
    #     self.loc = loc
    #     self.scale = scale

    def sample(self,N=1):
        x = np.random.normal(self.loc, self.scale, size=N)
        flipmask = np.ones(N, dtype=bool)
        flipmask[::2] = False
        x[flipmask] = -x[flipmask]
        return x


class SegmentDx(Gen):

    def __init__(self, dxseq, pwl):
        self.dxseq = dxseq
        self.cumulative_length = pwl.get_cumulative_length()
        # track the cumulative distance within this object
        self.v = 0

    def next_dx(self):
        # todo efficiency
        # which segment are we on
        index = np.searchsorted(self.cumulative_length, self.v)
        dx = self.dxseq[index]
        self.v += dx
        return dx
        
    def sample(self, N=1):
        if N == 1:
            return self.next_dx()
        else:
            return np.array([self.next_dx() for _ in range(N)])


def new_pwl_process(length, angle, N):
    # piecewise linear process generates a new segment using the direction of the previous segment
    l = length.sample(N)
    a = [0.]
    a.extend(angle.sample(N-1))
    x = [0., l[0]]
    y = [0., 0.]
    for i in range(1, N):
        p = np.array([x[-1]-x[-2], y[-1]-y[-2]])
        u = p/norm(p)
        # rotate
        cosa = np.cos(a[i])
        sina = np.sin(a[i])
        _ux = u[0] * cosa - u[1] * sina
        _uy = u[0] * sina + u[1] * cosa 
        _x = x[-1] + l[i] * _ux
        _y = y[-1] + l[i] * _uy
        x.append(_x)
        y.append(_y)
    return mdl.LPshape(np.array(x), np.array(y))

def new_static_process(length, angle, N):
    # static process generates new segments against e_x
    l  = length.sample(N)
    a = angle.sample(N)

    x = [0.]
    y = [0.]
    for i in range(1, N):
        _x = x[-1] + l[i] * np.cos(a[i])
        _y = y[-1] + l[i] * np.sin(a[i])
        x.append(_x)
        y.append(_y)
    return mdl.LPshape(np.array(x), np.array(y))


def sample_pwl(pwl, dx, error):
    v = 0
    max_v = pwl.get_contour_length()
    measurements = []
    while v < max_v:
        errxy = np.array([error.sample(), error.sample()]).ravel()
        pt = pwl(v) + errxy
        measurements.append(pt)
        v += dx.sample()
    DT = 0.1
    dt = np.full(len(measurements), DT)
    dt[0] = 0
    x, y = np.stack(measurements).T
    return mdl.LPtrack(dt, x, y)

def sample_pwl(pwl, dx, error):
    v = 0
    max_v = pwl.get_contour_length()
    measurements = []
    while v < max_v:
        errxy = np.array([error.sample(), error.sample()]).ravel()
        pt = pwl(v) + errxy
        measurements.append(pt)
        v += dx.sample()
    DT = 0.1
    dt = np.full(len(measurements), DT)
    dt[0] = 0
    x, y = np.stack(measurements).T
    return mdl.LPtrack(dt, x, y)

def new_ideal_synthetic(N, params):
    _l = params.get('l', 1.0)
    sigma = params.get('sigma', 0.10)
    length = Uniform(_l, _l)

    dx = Constant(_l/10)
    error = Normal(scale=sigma)

    # test mirrored normal 
    mnormal = AlternatingMirrorNormal(loc=pi/4, scale=pi/16)
    
    pwl = new_static_process(length, mnormal, N)
    synthdata = sample_pwl(pwl, dx, error)
    return pwl, synthdata

# ----------------------------------------------------------------
# plotting

def exactplot(ax, gen, N=1000):
    basis = np.linspace(-pi, pi, N)
    p = list(map(gen.pdf, basis))
    ax.plot(basis, p)
    def_blue = '#1f77b4'
    ax.fill_between(basis, 0, p, alpha=0.2, hatch="/", edgecolor=def_blue)
    ax.set_xticks([-pi,0,pi])
    ax.set_xticklabels([r'$-\pi$','0',r'$\pi$'])
    ax.set_xlabel(r'$\theta$')
    ax.set_ylabel(r'$P$')
