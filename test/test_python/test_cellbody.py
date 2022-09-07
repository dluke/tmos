
import os
import numpy as np
import scipy.special
from scipy.stats import vonmises

import tmos 
import ctfp3d
import parameters


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


#################################################################################
# helper method

# Why is equal aspect ratio not implemented in Matplotlib for 3D axis!?
#https://stackoverflow.com/questions/8130823/set-matplotlib-3d-plot-aspect-ratio
def axisEqual3D(ax):
    extents = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
    sz = extents[:,1] - extents[:,0]
    centers = np.mean(extents, axis=1) 
    maxsize = max(abs(sz))
    r = maxsize/2
    for ctr, dim in zip(centers, 'xyz'):
        getattr(ax, 'set_{}lim'.format(dim))(ctr - r, ctr + r)
#

#################################################################################
# test the c++ code itself and draw the distribution

def tonp(v):
    return np.array([v.x,v.y,v.z])

def get_one_pt():
    cell = setup_cell(2.0, 1.0)
    X = tonp(cell.pili_body_position())
    print(X)

R = 0.5
def setup_cell(pilivar=R * np.pi/2., polarisation=0.):
    args = parameters.read().apply()
    args.ACell.pilivar = pilivar
    args.ACell.polarisation = polarisation
    cell = ctfp3d.setup_cell(args)
    return cell

def gen_distrib(cell, N):
    def gen_pt(cell):
        return tonp(cell.pili_body_position())
    return np.stack([gen_pt(cell) for _ in range(N)])
    
def scat_distrib(pgen, distrib):
    p = next(pgen)

    fig = plt.gcf() 
    ax = fig.add_subplot(3,3,p, projection='3d')

    xs, ys, zs = np.split(distrib, 3, axis=1)
    
    ax.axis('off')
    ax.scatter(xs, ys, zs)
    setaspect = (1,1,1.5)
    ax.set_box_aspect(setaspect)
    aspect = ax.get_box_aspect()
    
    xlim = 0.6

    ax.set_xlim(-xlim, xlim)
    ax.set_ylim(-xlim, xlim)
    ax.set_zlim(0., 2*xlim*setaspect[2])

    # zdlim = (zs.min(), zs.max())
    # zdz = zdlim[1]-zdlim[0]
    # zmid = zdlim[0] + zdz/2.
    # zlims = zmid - aspect_ratio * zdz/2., zmid + aspect_ratio * zdz/2.
    # ax.set_zlim(*zlims)


def test_pili_body_distrib():
    N = int(1e3)

    fig = plt.figure(figsize=(10,15))
    bigax = fig.add_subplot(111)
    bigax.tick_params(bottom='off', left='off')
    bigax.set_xticklabels([])
    bigax.set_yticklabels([])
    bigax.set_xlabel("polarisation")
    bigax.set_ylabel("Distribution width")
    bigax.set_title("TEST: Generate Pili Anchor Points")

    sq = R * np.pi/2
    #pilivarls = [sq/4., sq, sq + 1.]  # uniform
    pilivarls = [0.5, 2, 7.0] # gaussian
    polls = [1., 0.8, 0.]
    pgen = iter(list(range(1, len(pilivarls)*len(polls)+1, 1)))
    for i, pilivar in enumerate(pilivarls):
        for j, pol in enumerate(polls):
            cell = setup_cell(pilivar, pol)
            distrib = gen_distrib(cell, N)
            scat_distrib(pgen, distrib)

    print("writing to pili_anchors.png")
    plt.savefig("pili_anchors.png")

def view3d():
    from mpl_toolkits.mplot3d import Axes3D
    test_pili_body_distrib()

#################################################################################
# Von mises distribution (2d)

def _probability_in_range(von, xlims):
    xn, xm = xlims
    return (von.cdf(xm) - von.cdf(xn))

def von_deriv(kappa, mu=0.):
    def vd(x):
        return -kappa*np.sin(x-mu)*np.exp( np.cos(x-mu) ) / 2*np.pi*scipy.special.j0(kappa)
    return vd

def _make_stretch_function(L):
    def g(x):
        return np.pi/L  * x + np.pi/2.*(1 - np.pi/L)
    return g

def make_composite_function(kappa, L):
    von = vonmises(kappa)
    g = _make_stretch_function(L)
    def modified_vh(x):
        if x <= np.pi/2.:
            return von.pdf(x)
        elif x > np.pi/2.:
            return von.pdf(g(x))
    def modified_von(x):
        return modified_vh(abs(x))
    return modified_von


def main():
    L = 2.
    fig, axes = plt.subplots(1,2, figsize=(14,5))
    ax1, ax2 = axes
    for ax in axes:
        ax.set_ylim(bottom=0)
        ax.set_xlim(-np.pi, np.pi)
        ax.tick_params(labelsize=15.)
        ax.axvline(-np.pi/2., c='r', linestyle='--')
        ax.axvline(np.pi/2., c='r', linestyle='--')

    lspace = np.logspace(-2, 2, 8)
    ilspace = np.linspace(-2, 2, 8)
    for kappa in lspace:
        draw_modified_von(ax1, kappa, L)
    for kappa in np.logspace(ilspace[3], ilspace[-2], 8):
        draw_modified_von(ax2, kappa, L)
    #draw_von(kappa)
    ax1.legend(fontsize=14.)
    ax2.legend(fontsize=14.)

def main2():
    kappa = 2.0 # shape parameter
    L = 4.0
    fig, axes = plt.subplots(1,2)
    fig.set_size_inches(11.5, 5)
    ax, clax = axes
    ax.tick_params(labelsize=15)
    clax.tick_params(labelsize=15)
    #
    ax.set_ylim(bottom=0)
    lim = np.pi/2. + L/2.
    ax.set_xlim(-lim, lim)
    ax.axvline(-np.pi/2., c='r', linestyle='--')
    ax.axvline(np.pi/2., c='r', linestyle='--')
    #
    clax.set_ylim(0, 0.2)
    clax.set_xlim(np.pi/2.-np.pi/10., lim)
    clax.axvline(np.pi/2., c='r', linestyle='--')

    draw_von(ax, kappa)
    line = draw_modified_von(ax, kappa, L)
    line.set_label('non-normalised modified von-mises')
    #line.set_linewidth(0.5)
    line.set_alpha(0.6)
    ax.legend(fontsize=15)

    draw_von(clax, kappa)
    line = draw_modified_von(clax, kappa, L)
    line.set_alpha(0.6)

    #draw_von(kappa)



def draw_von(ax, kappa): 
    von = vonmises(kappa)
    vd = von_deriv(kappa)
    #xlims = (-np.pi/2., np.pi/2.)
    #x = np.linspace(-np.pi/2., np.pi/2., 100) # angle
    largex = np.linspace(-np.pi, np.pi, 100) # angle
    #ax.plot(x, von.pdf(x), 'r-', lw=5, alpha=0.6, label='vonmises pdf')
    line, = ax.plot(largex, von.pdf(largex), 'r-', lw=1, alpha=0.6, label='vonmises pdf')
    return line
    

#################################################################################
# Von mises-fisher distribution (3d)

# 
def fw(w, kappa):
    return kappa/(2*np.sinh(kappa)) * np.exp(kappa * w) 

def uniform_v():
    theta = np.random.uniform(-np.pi, np.pi)
    return np.array([np.cos(theta), np.sin(theta)])

def rW(kappa):
    # inverse cumulative distrubution function method
    y = np.random.uniform()
    c = 2*np.sinh(kappa)/kappa
    return 1/kappa * np.log(np.exp(-kappa) + kappa * c * y)

def sample_vmf(kappa):
    vx, vy = uniform_v()
    W = rW(kappa)
    wf = np.sqrt(1-W**2)
    return np.array([wf*vx, wf*vy, W])

def _axis_line(ax):
    # create a line

    v = np.array([0,0,1])
    ext = 2
    m= np.linspace(-ext, ext, 200)
    line = np.vstack([mv*v for mv in m])
    ax.plot(line[:,0], line[:,1], line[:,2], c='k', alpha=0.8, linewidth=0.4)

def _sphere_surface(ax):
    # Create a sphere
    r = 1
    pi = np.pi
    cos = np.cos
    sin = np.sin
    density = 40j
    phi, theta = np.mgrid[0.0:pi:density, 0.0:2.0*pi:density]
    x = r*sin(phi)*cos(theta)
    y = r*sin(phi)*sin(theta)
    z = r*cos(phi)

    ax.plot_wireframe(
        x, y, z, rstride=1, cstride=1, color='c', alpha=0.6, linewidth=0.1)

def scat_vmf_distrib(sample):
    N = 1000

    fig = plt.gcf() 
    ax = fig.add_subplot(111, projection='3d')

    xs, ys, zs = np.split( np.vstack([sample() for _ in range(N)]), 3, axis=1 ) 

    r = 1.+0.1
    ax.set_xlim(-r, r)
    ax.set_ylim(-r, r)
    ax.set_zlim(-r, r)
    ax.axis('off')
    ax.scatter(xs, ys, zs, s=1.)
    _sphere_surface(ax)
    _axis_line(ax)


def scat_mod_distrib(ax, sample):
    N = 1000


    xs, ys, zs = np.split( np.vstack([sample() for _ in range(N)]), 3, axis=1 ) 

    r = 1.+0.1
    ax.set_xlim(-r, r)
    ax.set_ylim(-r, r)
    ax.set_zlim(-r, r)
    ax.axis('off')
    ax.scatter(xs, ys, zs, s=1.)
    _sphere_surface(ax)
    #_axis_line(ax)

    axisEqual3D(ax)


def _make_vmf_stretch_function(L):
    def g(x):
        return L/np.pi  * x + (np.pi-L)/2.
    return g


ez = np.array([0,0,1])
def modified_distrib(sample, L):
    #... if theta > pi/2. need to map onto the body length [R*pi/2., R*pi/2. + L/2.] ...
    # L is the body length, we would use L/2. for the half body execpt that R = 0.5 so use L
    X = sample()
    g = _make_vmf_stretch_function(2*L)

    theta = np.arccos(np.dot(X, ez))
    if theta > np.pi/2.:
        # get the distance around the unit sphere
        md = g(theta) - np.pi/2.
        # project X onto unit circle in x-y plane
        X[2] = 0.
        X = X/np.linalg.norm(X) 
        # extend X by md in -z direction
        X -= md * ez
    return X

def draw_this_distribution():
    kappa = 2.0
    sample = lambda: sample_vmf(kappa)

    scat_vmf_distrib(sample)

def draw_modified_distribution():
    kappa = 0.5
    L = 2.
    sample = lambda: sample_vmf(kappa)
    mod_sample = lambda: modified_distrib(sample, L)

    fig = plt.gcf() 
    fig.set_size_inches(11.5, 5)
    ax = fig.add_subplot(121, projection='3d')
    scat_mod_distrib(ax, mod_sample)
    ax2 = fig.add_subplot(122, projection='3d')
    ax2.view_init(90.,0)
    scat_mod_distrib(ax2, mod_sample)

    modfile = 'mod_distribution_kappa_{:08.4f}.png'.format(kappa)
    modfile = os.path.join('plots/', modfile)
    plt.savefig(modfile)


CB91_Blue = '#2CBDFE'
CB91_Green = '#47DBCD'
CB91_Pink = '#F3A0F2'
CB91_Purple = '#9D2EC5'
CB91_Violet = '#661D98'
CB91_Amber = '#F5B14C'

color_list = [CB91_Blue, CB91_Pink, CB91_Green, CB91_Amber,
              CB91_Purple, CB91_Violet]
r_color_list  = list(reversed(color_list))
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=r_color_list)

import pili.publication as pub
import matplotlib as mpl

def public():
    # plt.style.use('ggplot')
    mpl.rcParams['text.usetex'] = True
    L = 2.
    kspace = np.array([0.5,1.0,2.5,7.0,15.0])
    fig, ax  = plt.subplots(figsize=(7,5))
    def _plot(ax, kspace):
        ax.set_ylim(bottom=0)
        ax.set_xlim(-np.pi/2-1, np.pi/2+1)
        ax.tick_params(labelsize=24.)
        glstyle = {'color':'lightgrey', 'linestyle':'--' , 'linewidth':4, 'alpha':1.0}
        lstyle = {'linewidth':3}
        ax.axvline(-np.pi/2., **glstyle)
        ax.axvline(np.pi/2., **glstyle)

        for kappa in kspace:
            draw_modified_von(ax, kappa, L, lstyle)
        pi = np.pi
        # ax.set_xticks([-pi/2-1, -pi/2, 0, pi/2, pi/2 +1])
        ax.set_xticks([-pi/2, 0, pi/2])

        # ax.set_xticklabels([r'$-\pi$', r'$-\pi/2$', r'0', r'$\pi/2$', r'$\pi$'])
        ax.set_xticklabels([r'$-\pi/4$', r'0', r'$\pi/4$'])
        ax.grid(False)
    # plot first figure
    _plot(ax, kspace)
    ax.legend(fontsize=20.)
    notename = "pili/src/test_python/test_cellbody.py"
    pub.save_figure("mod_vonmises", notename)

    # these values are from original crawling ABC (mc4d)
    # triplet = [1.97, 2.99, 4.56]
    # outfile = "predicted_kappa"

    # these values are from crawling ABC with some improvements (N=50) (mc4d)
    # triplet = []
    # outfile = "crawling_predicted_kappa"
    
    # these values are from walking ABC with some improvements (N=50) (mc4d_walking)
    triplet = [1.10, 1.38, 2.22]
    outfile = "walking_predicted_kappa"

    # plot second figure
    dashstyle = {'color':color_list[-2]}
    kspace = np.array([]) # 1.97, 4.56
    fig, ax  = plt.subplots(figsize=(5,5))
    lim = np.pi/2. + L/2.
    _plot(ax, kspace)
    x = np.linspace(-lim, lim, 100)
    y = make_composite_function(triplet[1], L)
    _y = np.array([y(_x) for _x in x])
    handles = []
    l, = ax.plot(x, _y, color=color_list[-1], linewidth=3)
    handles.append(l)
    y = make_composite_function(triplet[0], L)
    y1 = np.array([y(_x) for _x in x])
    y = make_composite_function(triplet[2], L)
    y2 = np.array([y(_x) for _x in x])
    ax.fill_between(x, y1, y2, alpha=0.2)
    l, = ax.plot(x, y1, linestyle='--', **dashstyle)
    handles.append(l)
    l, = ax.plot(x, y2, linestyle='-.', **dashstyle)
    handles.append(l)
    labels = ['$\kappa = {:6.2f}$'.format(_v) for _v in [triplet[1], triplet[0], triplet[2]]]
    ax.legend(handles, labels, fontsize=16.)
    ax.set_xlabel("cylindrical distance $x$ ($\mu m$)", fontsize=24)
    ax.set_ylabel(r"$P(x | \kappa)$", fontsize=24)
    plt.tight_layout()

    pub.save_figure(outfile, notename)

    # ...
    # draw the shapes for 2.99 [1.97, 4.56]
    
def draw_modified_von(ax, kappa, L, lstyle={}): 
    modified_von = make_composite_function(kappa, L)
    lim = np.pi/2.+L/2.
    x = np.linspace(-lim, lim, 100) 
    line, = ax.plot(x, list(map(modified_von, x)), label='$\kappa = {:6.2f}$'.format(kappa), **lstyle)
    return line



if __name__=='__main__':

    #0.
    # view3d()
    #get_one_pt()

    #1.
    #main()
    #plt.savefig('plots/von_mises_distribution.png')
    #plt.show()
    
    # 2.
    #N = 1000
    #for _ in range(N):
        #X = sample()
        #print np.linalg.norm(X)
 
    # 3.
    #draw_this_distribution()
    # main2()
    # plt.savefig('plots/von_mises_modified.png')

    # 4.
    #draw_modified_distribution()

    # 5.
    # publication: plot distribution
    public()
    plt.savefig('publish_distribution.png')


