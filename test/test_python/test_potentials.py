

import tmos
import tmos.surface as surface
# test surface.potentials

import numpy as np
import matplotlib.pyplot as plt

# for finding the height of the force barrier
import scipy
import scipy.optimize
import scipy.interpolate

#print surface.va_force(1., 0.5)

delta = 0.01
eps = 10.
re = 0.49
R = 0.5
inside = 0.04

def test_LJ():
    rspace = np.linspace(R-0.05, R, 1000)
    def lj(r): 
        return surface.wca_energy_abs(eps, R, r)
    def lj_f(r):
        return surface.wca_force_abs(eps, R, r)

    print('energy at 0.04 units inside the particle')
    #r = 0.497
    #r = 0.4999
    #print lj(r)
    #print 'force'
    #print lj_f(r)
    print('Force at the anchor point >> F_stall', lj_f(R - inside))

    #plt.clf()
    fig, (ax1, ax2) = plt.subplots(1,2)
    ax1.plot(rspace, list(map(lj, rspace)), label='lj energy')
    ax1.legend()
    ax2.plot(rspace, list(map(lj_f, rspace)), label='lj force')
    ax2.legend()
    plt.show()
    

def make_get_max_force(a, b):
    def get_max_force(eps):
        sp = surface.Spotential(eps, re, R)
        faeps = lambda r: sp.sforce(r) # gradient
        result = scipy.optimize.fminbound(faeps, re, R)
        maxf = abs(faeps(result))
        return maxf
    return get_max_force

def construct_maxf_functions(get_max_force, basislims=(10, 510), step=10):
    emin, emax = basislims
    eps_space = np.arange(emin, emax, step)
    xmaxf = list(map(get_max_force, eps_space))
    # forward function
    forward = scipy.interpolate.interp1d(eps_space, xmaxf)
    reverse = scipy.interpolate.interp1d(xmaxf, eps_space)
    return forward, reverse

# wrapper which takes ParameterList directly
def maxf_functions(args):
    b = args.cell.R
    a = b - args.Cell3d.attract_width
    return construct_maxf_functions(make_get_max_force(a, b))

def check_force_at_anchor():
    eps_list = np.array([0.01, 0.1, 1, 10, 100])
    re = 0.49
    R = 0.5
    get_max_force = make_get_max_force(re, R)
    for eps in eps_list:
        sp = surface.Spotential(eps, re, R)
        print('for eps = {} force at anchor = {}'.format(eps, sp.sforce(0.46)))


def test_poly():

    eps = 10.
    re = 0.49
    R = 0.5

    eps = 10.
    eps_list = np.array([10., 20., 50., 80., 100., 200., 500.])
    maxf = np.zeros(eps_list.size)

    fig, (ax1, ax2) = plt.subplots(1,2)
    get_max_force = make_get_max_force(re, R)

    def plot_fa(eps):
        sp = surface.Spotential(eps, re, R)
        faeps = lambda r: -sp.sforce(r) # gradient
        eneps = lambda r: sp.senergy(r)
        rspace = np.linspace(re-0.01, R+0.01, 1000)
        #plt.plot(rspace, map(eneps, rspace), label='energy')
        ax1.plot(rspace, list(map(faeps, rspace)), label='gradient')
        ax2.plot(rspace, list(map(eneps, rspace)), label='energy')
        return get_max_force(eps)

    for i, eps in enumerate(eps_list):
        thismaxf = plot_fa(eps)
        maxf[i] = thismaxf

    ax1.legend()
    ax2.legend()
    plt.savefig("polynomial_attractive_potential.png")
    #plt.show()
    plt.clf()

    forward, reverse = construct_maxf_functions(get_max_force)
    eps_space = np.arange(10, 510, 10)

    fig, (ax1, ax2) = plt.subplots(1,2)
    maxf_space = list(map(forward, eps_space))
    ax1.plot(eps_space, maxf_space)
    print(eps_space)
    print(maxf_space)
    ax1.set_xlabel("\epsilon")
    ax1.set_ylabel("max attractive force")
    ax2.plot(maxf_space, eps_space)
    ax2.set_xlabel("max attractive force")
    ax2.set_ylabel("\epsilon")
    plt.tight_layout()
    plt.show()

test_poly()

def test_smooting_function():
    def make_lj(eps, rmin):
        def lj(r):
            sr = (rmin/r)**6
            return eps*(sr*sr - 2*sr)
        return lj

    def make_lj_deriv(eps, rmin):
        def lj_deriv(r):
            sr = (rmin/r)**6
            return -12 * eps * (1/r) * ( sr*sr - sr )
        return lj_deriv


    def make_smoothing_function(a, b):
        def smoothing_function(r):
            if r < a:
                return 1.
            elif r > b:
                return 0.
            else:
                rsas = r*r - a*a
                bsas = b*b - a*a
                return 1 + (rsas/bsas)**2 *( 2 * (rsas/bsas) -3 )
        return smoothing_function

    def make_smoothing_derivative(a, b):
        def smoothing_deriv(r):
            if r < a:
                return 0.
            if r > b:
                return 0.
            rsas = r*r - a*a
            bsas = b*b - a*a
            return 6*(rsas*2*r/bsas)*( rsas/bsas - 1 )
        return smoothing_deriv

    def make_smoothing_deriv2(a, b):
        def sderiv2(r):
            if r < a:
                return 0.
            if r > b:
                return 0.
            bsas = b*b - a*a
            rsas = r*r - a*a
            R2 = rsas/bsas
            R1 = 2*r/bsas
            return 6*( R1*R1*(R2-1) + 2*R2*(R2-1) + R2*R1*R1 )
        return sderiv2

    a = 0.490
    b = 0.5
    smoothing_function = make_smoothing_function(a, b)
    smoothing_deriv = make_smoothing_derivative(a, b)
    smoothing_deriv2 = make_smoothing_deriv2(a, b)
    
    eps = 10.
    rmin = 0.49
    lj = make_lj(eps, rmin)
    lj_deriv = make_lj_deriv(eps, rmin)

    def squeeze_lj(r):
        return lj(r) * smoothing_function(r)
    def squeeze_lj_deriv(r):
        if r < a:
            return lj_deriv(r)
        elif r > b:
            return 0.
        else:
            return lj_deriv(r) * smoothing_function(r) + lj(r) * smoothing_deriv(r)

    print('at b', smoothing_function(b))

    # global settings
    params = {'legend.fontsize': 20, 'axes.titlesize':30, 'font.size':20}
    plt.rcParams.update(params)

    rspace = np.linspace(0.48, 0.51, 1000)
    big_rspace = np.linspace(0.45, 0.55, 1000)
    lj_rspace = np.linspace(0.40, 0.75, 1000)
    fig, axes = plt.subplots(2,3)
    fig.set_size_inches((15, 10))
    ax1, ax2, ax3 = axes[0,:]
    ax1.plot(lj_rspace, list(map(lj, lj_rspace)), label='lj')
    ax1.set_title('LJ')
    ax2.plot(rspace, list(map(smoothing_function, rspace)), label='smoothing')
    ax2.set_title('S(r)')
    ax3.plot(lj_rspace, list(map(squeeze_lj, lj_rspace)), label='Energy')
    ax3.set_title('LJ(r)S(r)')

    # derivatives
    bx1, bx2, bx3 = axes[1,:]
    bx1.plot(lj_rspace, list(map(lj_deriv, lj_rspace)), label='lj deriv')
    bx1.set_title('LJ Deriv')
    bx2.plot(rspace, list(map(smoothing_deriv, rspace)), label='smoothing deriv')
    bx2.set_title('S(r) Deriv')
    bx3.plot(rspace, list(map(squeeze_lj_deriv, rspace)), label='Gradient')
    bx3.set_title('Gradient')
    #fig.suptitle("Constructing Surface Potential")
    
    for ax in axes.flatten():
        ax.axvline(a, c='grey', alpha=0.3)
        ax.axvline(b, c='grey', alpha=0.3)

    plt.tight_layout()
    plt.savefig("smoothing_function_attractive_potential.png")
    #plt.show()

    plt.clf()
    plt.plot(rspace, list(map(smoothing_deriv2, rspace)), label='smoothing second derivative')
    #plt.show()

    def get_second_deriv_roots(a,b):
        ca = 5
        cb = -3*(b**2+a**2)
        cc = a*a*b*b
        sqdet = np.sqrt(cb**2 - 4*ca*cc)
        r1, r2 = (-cb - sqdet)/(2*ca), (-cb + sqdet)/(2*ca)
        return np.sqrt(r1), np.sqrt(r2) # only interested in the +ve roots
    r1, r2 = get_second_deriv_roots(a,b)
    #print squeeze_lj(r2), squeeze_lj_deriv(r2)
    # not roots ? must be mistake in paper calculation of derivative
    # come back to this
    function = lambda r: -squeeze_lj_deriv(r)

    #result = scipy.optimize.fminbound(function, a, b)
    #print function(result)
    #print result
    
    def make_test_poly(a,b):
        def test_poly(r):
            R = r*r
            return 5*R*R - 3*(b**2+a**2)*R + a*a*b*b
        return test_poly
    test_poly = make_test_poly(a,b)
    plt.clf()
    plt.plot(rspace, list(map(test_poly, rspace)), label='smoothing second derivative')
    # plt.show() 




#test_LJ()
#test_poly()
#check_force_at_anchor()

# test_smooting_function()

def plot_smoothing_function():
    a = 1
    b = 3
    show_deriv = True

    def make_smoothing_function(a, b):
        bsas = b*b - a*a
        def smoothing_function(r):
            if r < a:
                return 1.
            elif r > b:
                return 0.
            else:
                rsas = r*r - a*a
                rb = rsas/bsas
                return 1 + rb*rb * ( 2 * (rb) - 3)
        return smoothing_function
    smooth = make_smoothing_function(a, b)

    def smooth_(x):
        return 1 + x**2 * (2*x - 3)

    def new_make_smoothing(a, b):
        def smoothing_function(r):
            t = (r - a)/(b - a)
            return 1 + t**2 * (2*t - 3)
        return smoothing_function
    newsmooth = new_make_smoothing(a, b)

    # edgecase = new_make_smoothing(2,2)
    # assert(edgecase(2) == 0), edgecase(2)

    def new_make_deriv(a, b):
        def deriv(r):
            t = (r - a)/(b - a)
            return 6*t*(t-1)
        return deriv
    newderiv = new_make_deriv(a, b)

    xspace = np.linspace(a,b,1000,True)
    y = [smooth(x) for x in xspace]
    y_ = [smooth_(x) for x in xspace]
    newy = [newsmooth(x) for x in xspace]
    newdy = [newderiv(x) for x in xspace]
    # plt.plot(xspace, y, alpha=0.4)
    # plt.plot(xspace, y_, alpha=0.4, label='')
    plt.plot(xspace, newy, alpha=0.4)
    if show_deriv:
        plt.plot(xspace, newdy, alpha=0.4)
    plt.show()

# plot_smoothing_function()