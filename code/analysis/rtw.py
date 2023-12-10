#!/usr/bin/env python3
import os, sys
import numpy as np
import scipy.stats
import copy
from collections import OrderedDict

import command
from command import defaultsave

import filesystem
import readtrack
import readmat
import readyr
import parameters 
import txtdata

import stats
import plotutils
import twutils
import shapeplot
import twanalyse
import eventanalyse
import astrack
import analyse2d
import _fj
import fjstats
import yrstats
import simstats 

"""matplotlib 3d setup"""
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cmx
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import matplotlib.transforms as mplt


#################################################################################



def tw_walk(ddirs, tw, *args):
    n = len(ddirs)
    cwd = os.getcwd()
    for i, dd in enumerate(ddirs):
        plt.subplot(1,n,i+1)
        print(("Moving to directory", dd))
        os.chdir(dd)
        tw(*args)
        os.chdir(cwd)


import itertools
"""
DataCube finds directories loads and plots data
Break this down into 
->parsing directories
->providing an interface to run analysis on any directory
->storing and plotting the results
"""
#cubedir = os.path.join(command.ddir, 'cube/')
class DataCube():


    def __str__(self):
        # quick overview of the parameters 
        sl = [
            "DataCube inspected directories at --> {}".format(self.targetdir), 
            "parameters: {}".format(self.pnames),
            "with shape: {}".format(self.basisdims),
            "with basis: ",
        ]
        for i, p in enumerate(self.pnames):
            sl.append("{} = {}".format(self.pnames[i], self.basis[i]))
        sl.append('<--')
        return '\n'.join(sl)

    def split_parameter_dir(self, dname):
        # using parameters.meta data
        nvls = filesystem.par_split(dname, self.pnames)
        names = nvls[::2]
        vals = [eval(s.lstrip(' 0')) for s in nvls[1::2]]
        return names, vals

    def __init__(self, target='.', cubedir='plots/data/cube/'):
        """
        On the current directory full of generated simulation directories
        with grid like parameter sets in n (3) dimensions.
        Construct an object which lets us run some analytical routine
        like those defined in twanalyse.

        Those routines should take no arguments and be run from the directory
        with data/

        Members:
            basis - list of lists of parameter values
            ibasis - list of lists of indices which span the parameter space
            tmpbasis - subset basis
            dircube - Numpy array containing the directories of the data
            pbounds - An ordered dictionary of bounds {parameter_name -> (from, to)}
            basisdims - whole shape
            subdims - shape of current subset

        Future:

        generalise to aribtrary dimensions
        sensible output for any dimensions

        """
        self.targetdir = target
        ddirs = filesystem.getdirs(target)
        self.pnames = list(filesystem.read_param_meta(target).keys())
        self.modpnames = list(map(plotutils.fix_tex_label, self.pnames))
        self.dimension = len(self.pnames)

        pvals = np.array([self.split_parameter_dir(dd)[1] for dd in ddirs])
        # [parameter1 basis, parameter2 basis...]
        self.basis = [sorted(list(set(pvalset))) for pvalset in pvals.T]
        self.basisdims = list(map(len, self.basis))

        directory_lookup = OrderedDict(list(zip( list(map(tuple, pvals)), ddirs)))

        #for key, dirname in directory_lookup.items():
        self.ibasis = list(itertools.product( *[list(range(len(x))) for x in self.basis] ))
        self.slice_ibasis = copy.deepcopy(self.ibasis)
        self.slice_ioffset = np.zeros(self.dimension, dtype=int)
        self.slice_basis = copy.deepcopy(self.basis)
        self.pbounds = {name: (0, dim) for name, dim in zip(self.pnames, self.basisdims)}
        self.subdims = copy.copy(self.basisdims)

        dircube = self.new_cube(dtype=object)
        pcube = self.new_cube(dtype=object)
        for ipt, pt in zip(self.ibasis, itertools.product(*self.basis)):
            dd = directory_lookup[tuple(pt)]
            dircube[ipt] = os.path.join(target, dd)
            pcube[ipt] = self.split_parameter_dir(dd)[1]
        self.dircube = dircube
        self.pcube = pcube
        """dircube contains simulation directories"""
        """pcube contains set parameter list"""

        """Setup a directory to save and load numpy array data cubes """
        cubedir = os.path.join(target, cubedir)
        if not os.path.exists(cubedir):
            print(('making directory {}'.format(cubedir)))
            os.makedirs(cubedir)
        self.cubedir = cubedir

        self.trcube = None
        self._lldata = None


    def prettynames(self):
        return [txtdata.prettynames.get(pname, 'no pretty name') for pname in self.pnames] 

    def find_index(self, vset):
        # given a point determined by parameter values find the index
        # the inverse operation is a lookup on self.pcube[vset]
        return [np.searchsorted(self.basis[i], v) for i, v in enumerate(vset)]

    """set bounds on the data we plan to analyse"""
    def set_idx_bounds(self, pbls):
        """Take a list of ranges [(from, to)] so that range(from, to) gives 
        the desired subset
        """
        if pbls == None:
            return # accept None as argument
        self.pbounds = OrderedDict(list(zip(self.pnames, pbls)))
        self.subdims = [t[1]-t[0] for t in list(self.pbounds.values())]
        # sequences of indexes in the full data scope
        self.slice_ibasis = list(itertools.product( *[list(range(*t)) for t in pbls] ))
        self.slice_ioffset = np.array([pb[0] for pb in pbls])
        self.slice_basis = [basis[t[0]:t[1]] for basis, t in zip(self.basis, pbls)]

        # TODO: surely we need to slice more member variables ...
        # TODO: the sliced basis should be calle self.basis and the original basis should be saved 

    """convienience methods"""
    # slice a small part of the data for testing 
    def set_mini(self):
        self.set_idx_bounds([(0,2)]*self.dimension)

    def new_cube(self, dtype=float):
        return np.empty(self.basisdims, dtype=dtype)
    def new_subcube(self, dtype):
        return np.empty(self.subdims, dtype=dtype)
 
    #
    def reload_data(self, dataname):
        """
        load specific plotting data from a simulation in the grid
        """
        def loader(x):
            datafile = os.path.join(x, 'plots/data', dataname)
            return np.loadtxt(datafile)
        datacube = self.new_cube(dtype=object)
        for idx in self.slice_ibasis:
            data = loader(self.dircube[idx])
            datacube[idx] = data
        return datacube

    def sortdir(self, cube):
        """Sort the simulation directories by values in cube"""
        sdargs = np.argsort(cube.flatten())
        return self.dircube.flatten()[sdargs]

    def sort_by(self, cube, metric):
        metriccube = self.new_cube(float)
        for idx in self.slice_ibasis:
            value = metric(cube[idx])
            metriccube[idx] = value
        if hasattr(metric, 'name'):
            saveable = os.path.join(self.cubedir, metric.name)
            np.save(saveable, metriccube)

        # sort on flattended array
        flatidx = np.argsort(metriccube, axis=None)
        return metriccube, flatidx

    def closeparams(self, metriccube, flatidx, cutoff):
        """
        flatidx is the sorted indices of metriccube
        """
        closeidx = [self.ibasis[i] for i in flatidx if metriccube.ravel()[i] < cutoff]
        closeparams = np.array(
                [[self.basis[j][k] for j,k in enumerate(idx)] for idx in closeidx])
        return closeidx, closeparams


   
    # concatenate cubes which contain lists
    def sumcubes(self, *args):
        # create new cube of empty lists
        c0 = args[0]
        newcube = np.zeros_like(c0)
        for idx in self.ibasis:
            newcube[idx] = copy.deepcopy(c0[idx])
        #
        for cube in args[1:]:
            for ipt in self.ibasis:
                newcube[ipt].extend(cube[ipt])

        return newcube

    """Run calculations"""

    def apply(self, method, trcube, dtype=object, mvdir=False):
        """
        trcube is a numpy object array of list(Track) objects
        method is a function which operates on lists of Tracks
        if mvdir is True then run method in the simulation directory
        """
        newcube = self.new_subcube(dtype=dtype)
        cwd = os.getcwd()
        for idx in self.slice_ibasis:
            idx = np.array(idx)
            mname = method.__name__ if hasattr(method, '__name__') else type(method)
            print(("computing method {} for {}".format(mname, idx)))
            if mvdir:
                os.chdir(self.dircube[idx])
                try:
                    # use get_subindex here?
                    res = method(trcube[idx])
                    newcube[self.get_subindex(idx)] = res
                finally:
                    os.chdir(cwd)
            else:
                # use get_subindex here?
                res = method(trcube[idx])
                
                newcube[self.get_subindex(idx)] = res
        return newcube

    def get_subindex(self, idx):
        idx = np.array(idx)-self.slice_ioffset
        idx = tuple(idx) if len(idx) > 1 else idx[0]
        return idx

    def calculate(self, tw_method, dtype=object):
        """
        return a numpy array with values determined by running tw_method 
        in each directory
        tw_method must run with no arguments
        """
        cube = self.new_subcube(dtype)
        cwd = os.getcwd()
        for idx in self.slice_ibasis:
            os.chdir( self.dircube[idx] )
            try:
                # fails if tw_method requires arguments
                cube[self.get_subindex(idx)] = tw_method()
            finally:
                os.chdir(cwd)
        return cube

    def load_local(self):
        if self._lldata is None:
            with command.chdir(self.targetdir):
                self._lldata = self.calculate(stats.load)
        return self._lldata

    def get_local_array(self, getter):
        lldata = self.load_local()
        C = np.zeros(self.subdims)
        for idx in self.slice_ibasis:
            s_idx = self.get_subindex(idx)
            ld = lldata[s_idx]
            try: 
                C[s_idx] = getter(ld)
            except (KeyError, AttributeError):
                C[s_idx] = np.nan
        return C

    def cube_on_trcube(self, method, dtype=float, force=False):
        """
        method argument operates on trs = readtrack.trackset()
        """
        trcube = self.get_trcube(method)
        return self.apply(method, trcube, dtype=dtype)

    def autocalculate(self, method, dtype=object, force=False):
        """
        close to autoloading but method as no arguments and gathers 
        track data or pilus data
        """
        saveable = os.path.join(self.cubedir, method.__name__)
        loadable = saveable+'.npy'
        recalculate = (not os.path.exists(loadable)) or force
        if recalculate:
            valuecube = self.calculate(method, dtype=dtype)
            np.save(saveable, valuecube, allow_pickle=True)
        else:
            valuecube = np.load(loadable, allow_pickle=True)
        return valuecube

    def autoloading(self, method, trcube, force=False, dtype=object):
        """
        save data cubes as .npy files in the current directory unless 
        a .npy file with that name exists in which case load it instead

        autoloading(a, b) signature should mirror apply(a, b)
        """
        saveable = os.path.join(self.cubedir, method.__name__)
        loadable = saveable+'.npy'
        recalculate = (not os.path.exists(loadable)) or force
        if recalculate:
            try:
                #valuecube = self.cube_on_trcube(method, dtype=dtype, force=force)
                valuecube = self.apply(method, trcube, dtype=dtype)
            except TypeError:
                print(('Bad call signature for method', method.__name__))
                raise

            np.save(saveable, valuecube, allow_pickle=True)
        else:
            valuecube = np.load(loadable, allow_pickle=True)
        return valuecube


    """
    Plotting routines for 3d data --- separate to another folder.
    Can subclass DataCube and add them there.
    """
    def scatter_cube(self, cube, colorlims=None):
        if self.dimension == 3:
            self.scatter_3d(cube, colorlims)
        else:
            print('scatter cube method is for 3d scatter plot')

    # one dimensiona
    def linsubplots(self, cube, cube_name=''):
        base = self.basis[0]
        name = self.modpnames[0]
        ax = plt.gca()
        ax.plot(base, cube.flatten(), marker='o', linestyle='None')
        ax.set_xlabel(name)
        ax.set_ylabel(cube_name)
        plt.tight_layout()

    # three dimensions
    def subplots(self, cube, subxy=[0,1,2], cube_name=''):
        axmap = dict(list(zip(list(range(3)),subxy)))
        # swap axes ...
        n1, n2, n3 = self.basisdims
        base1, base2, base3 = self.basis
        name1, name2, name3 = self.modpnames

        fig, axes = plt.subplots(1, n1, figsize=(5*n1,10), sharey=True)
        for ix, ax in enumerate(axes):
            for linen in range(n2):
                line = cube[ix,linen,:]
                ax.plot(base3, line, marker='o', linestyle='None',
                        label='{} = {}'.format(name2,base2[linen]))

            ax.set_xlabel(name3)
            ax.set_ylabel(cube_name)
            box = ax.get_position()
            scale = 0.5
            ax.set_position([box.x0, box.y0 + box.height * scale,
                             box.width, box.height * (1-scale)])
            # Put a legend below current axis
            ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15)) 
            ax.set_title('{} = {}'.format(name1, base1[ix]))


# https://matplotlib.org/2.2.5/tutorials/intermediate/gridspec.html?highlight=customizing%20figure%20layouts%20using%20gridspec%20other%20functions
    def lingrid(self, cube, pfunction, big_xlabel='', pcube=None, 
        sharey=False, prj='rectilinear', sortby=None):
        # figure and large axis
        width = 5.
        fig = plt.figure(figsize=(cube.size * width, 2.0 * width))
        # bbox = mplt.Bbox.from_bounds(1.,1.,cube.size*width, 1.5*width)
        big_grid = fig.add_gridspec(1,1)
        with plt.style.context([plotutils.get_style('bigax')]):
            bigax = fig.add_subplot(big_grid[0])
        subgrid = big_grid[0].subgridspec(1,cube.size,wspace=0.,hspace=0.)
        # plotutils.invisax(bigax)
        bigax.set_xlabel(big_xlabel)
        bigax.yaxis.set_ticklabels([])
        bigax.spines['left'].set_visible(False)
        bigax.spines['right'].set_visible(False)
        bigax.spines['top'].set_visible(False)
        offset = 1./cube.size/2.
        bigax.set_xticks(np.linspace(offset,1-offset,cube.size,True))

        if pcube is not None:
            bigax.set_xticklabels([p[0] for p in self.pcube])

        subplot_kw={'projection':prj}
        lincube = cube.flatten() if sortby is None else cube.flatten()[sortby]
        with plt.style.context(['fivethirtyeight', plotutils.get_style('ft8')]):
            linax = subgrid.subplots(sharey=sharey, subplot_kw=subplot_kw)
            for i, ax in enumerate(linax):
                ax = pfunction(ax, lincube[i])
        return bigax, linax

    # two dimensions
    def subplots2d(self, cube, cube_name=''):
        pass 

    # three dimension and might be generalised for 2 and 4 dimensions
    def gridspec(self, cube, pfunction, prj='rectilinear', 
            ccube=None, scm=None):
        """
        DataCube cube containing objects which can be used as arguments to a plotting function
        pfuntion a plotting function taking matplotlib axis (ax, cube.element) -> ax
        """
        n1, n2, n3 = cube.shape
        base1, base2, base3 = self.basis
        name1, name2, name3 = self.modpnames

        size = 10
        fig = plt.figure(figsize=(n3*size+12.0, size+5.0))
        innergls = []

        outer_grid = gridspec.GridSpec(1, n3, bottom=0.25, wspace=0.2, hspace=0.5)

        def setup_container(ax):
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.tick_params(axis='both', which='major', labelsize=40)
            ax.patch.set_alpha(0.)

        inner_ax = []
        for i in range(n3):
            inner_grid = gridspec.GridSpecFromSubplotSpec(n1, n2,
                    subplot_spec=outer_grid[i], wspace=0.15, hspace=0.15)
            in_rect = inner_grid[:,:].get_position(fig)
            in_ax = fig.add_axes(in_rect)
            in_ax.patch.set_alpha(0.)
            inner_ax.append(in_rect)
            innergls.append(inner_grid)
            for j, (c, d) in enumerate(itertools.product(list(range(n1)),list(range(n2)))):
                #ax = plt.Subplot(fig, inner_grid[j])
                args = parameters.thisread(directory = self.dircube[c,d,i])
                ax = plt.subplot(inner_grid[c, d], projection=prj)
                ax.patch.set_alpha(0.)
                cubecdi = cube[c, d, i]
                linekw = {} if ccube is None else {'color':ccube[c,d,i]} 
                pfunction(ax, cubecdi, linekw=linekw, args=args)
                fig.add_subplot(ax)
                #grl = plotutils.Grlabels(outer_grid, innergls, cube.shape)


        #inner_ax[0].set_ylabel(name1)d
        quad = np.linspace(1./8., 1.-1./8., n1, True)
        inaxrect = inner_ax[0]
        laxrect = inner_ax[-1]
        for ax_rect in inner_ax:
            in_ax = fig.add_axes(ax_rect)
            setup_container(in_ax)
            in_ax.set_xlabel(name2, fontsize=40)
            in_ax.spines['left'].set_visible(False)
            in_ax.set_xticks(quad)
            in_ax.set_xticklabels(base2)
            in_ax.tick_params(axis='x', which='major', pad=40)
            in_ax.set_yticklabels([])

        rect = outer_grid[:].get_position(fig)
        rect.y0 /=  2.
        outer_ax = fig.add_axes(rect)
        outer_ax.tick_params(axis='x', which='major', pad=20)
        outer_ax.tick_params(axis='y', which='major', pad=120)
        outer_ax.spines['left'].set_visible(False)
        outer_ax.set_ylabel(name1, fontsize=40)
        outer_ax.set_xlabel(name3, fontsize=40)
        outer_ax.set_xticks(np.linspace(0.1, 1.-0.1, n3, True))
        outer_ax.set_yticks(quad + rect.y0/2.)
        outer_ax.set_xticklabels(base3)
        outer_ax.set_yticklabels(base1[::-1])
        setup_container(outer_ax)

        if ccube is not None:
            x1off = 1-rect.x1
            cbar_rect = [rect.x1+x1off/4., rect.y0, x1off/4., rect.height]
            cbar_ax = fig.add_axes(cbar_rect)

            fig.colorbar(scm, cax=cbar_ax)

        return fig

    # used where?
    def get_basispts(self, idxarr):
        basispts = []
        for row in idxarr:
            basispts.append([self.basis[i][idx] for i, idx in enumerate(row)])
        return np.vstack(basispts)

        
    def fmethod(self, method, trcube, ylabel=None, force=False, nameex=''):
        """
        Run the method on the trcube and plot everywhere
        """
        dc = self
        if ylabel == None:
            ylabel = method.__name__.replace('_', '\_') # try to convert to raw-string
        cube = dc.autoloading(method, trcube, force=force)
        # dimension agnostic
        kw = {'cube':cube, 'cube_name':ylabel}
        if cube.ndim == 1:
            dc.linsubplots(**kw)
        elif cube.ndim == 3:
            dc.subplots(**kw)
        else:
            print('no plotting function for {} dimensional data'.format(cube.ndim))
            sys.exit()

        fullname = '_'.join([method.__name__, nameex]) +'.png'
        sf = os.path.join(command.pdir, fullname)
        plt.savefig(sf)

    # yet another plotting method ...
    def plotmethod(self, ax, pmethod, lldata=None, idx_order=[0,1], lstyle=None, pkw={}):
        # lldata = self.load_local() if lldata==None else lldata
        lldata = self.load_local()
        if lstyle==None:
            if self.dimension == 2:
                cm = default_cmap
                ncolors = len(self.slice_basis[idx_order[1]])+1
                colors = [cm(float(i)/ncolors) for i in range(ncolors)]
                lstyle = [{'color':c} for c in colors]
            elif self.dimension == 1:
                lstyle = [{'color':"blue"}]

        primary_idx, secondary_idx = idx_order
        print('Primary axis is {}'.format(self.pnames[primary_idx]))
        ax.set_xlabel(txtdata.longnames[self.pnames[primary_idx]])

        style = itertools.cycle(lstyle)
        if self.dimension == 1:  
            pmethod(ax, self, 0, lldata, linekw=style, **pkw)
        elif self.dimension == 2:
            second_parameter = self.pnames[secondary_idx]
            sec_basis = self.slice_basis[secondary_idx]
            # loop
            for i, basis_val in enumerate(sec_basis):
                name = txtdata.prettynames.get(second_parameter, second_parameter)
                label = '{} = {}'.format(name, sec_basis[i])
                ldata = lldata[i] if primary_idx == 1 else lldata[:,i]
                pmethod(ax, self, primary_idx, ldata, label=label, linekw=style, **pkw)
            ax.legend()


default_cmap = plt.get_cmap('gist_rainbow')

#################################################################################
# utilities


def var_or_nan(ldata, get):
    var = []
    for ld in ldata:
        try:
            var.append(get(ld))
        except KeyError:
            var.append(np.nan)
    return np.array(var)


##########

# plotutils.default_style()

# local definition of when to use log scale on x axis
uselog = {
    # 'pilivar': True
    # 'k_sp': True
}

# TODO: @ method for applying idx_bounds

# TODO: plot velocity decay times
@defaultsave()
def velocity_decay_time(idx_bounds=None):
    # plotutils.default_style()
    dc = DataCube()
    if idx_bounds != None:
        dc.set_idx_bounds(idx_bounds)

    fig = plt.figure()
    ax = plt.gca()
    lstyle = [{'marker':'D'}]
    idx_order = [1,0]
    dc.plotmethod(ax, _velocity_decay_time, idx_order=idx_order, lstyle=lstyle)
    ax.set_ylim(ymin=0.)

def _velocity_decay_time(ax,  dc, param_idx, ldata, label='', linekw=None):
    decay_time = [ld['velocity_decay_time'] for ld in ldata]
    basis = dc.slice_basis[param_idx]
    ax.plot(basis, decay_time, label=label, **next(linekw))
    ax.set_ylabel(r'velocity autocorrelation time $(s)$')
    
# 
@defaultsave()
def kmsd():
    dc = DataCube()
    def _kmsd(ax, dc, param_idx, ldata, label='', linekw=None):
        mean_kmsd = np.array([ld['kmsd']['mean'] for ld in ldata])
        kmsd_std_error = np.array([ld['kmsd']['std_error'] for ld in ldata])
        basis = dc.slice_basis[param_idx]
        this_style = next(linekw)
        ax.plot(basis, mean_kmsd, label=label, marker='D', **this_style)
        ax.errorbar(basis, mean_kmsd, 1.96*kmsd_std_error, ls='none', elinewidth=1, **this_style)
        #
        ax.set_ylim(ymax=2.0)

    
    fig = plt.figure()
    idx_order = [0,1]
    ax = fig.gca()
    if uselog.get(dc.pnames[idx_order[0]], False):
        ax.set_xscale('log')
    dc.plotmethod(ax, _kmsd, idx_order=idx_order)

    idx_order = [1, 0]
    fig = plt.figure()
    ax = fig.gca()
    if uselog.get(dc.pnames[idx_order[0]], False):
        ax.set_xscale('log')

    fig._plot_name = 'swapped'
    dc.plotmethod(ax, _kmsd, idx_order=idx_order)

@defaultsave()
def meanvel():
    # plotutils.default_style()
    dc = DataCube()
    # TODO
    # dc.set_idx_bounds([(0,11),(1,13)])
    def _meanvel(ax, dc, param_idx, ldata, label='', linekw={}):
        basis = dc.slice_basis[param_idx]
        velocity = var_or_nan(ldata, lambda ld: ld['l_mean_velocity'])
        style = next(linekw)
        style.update({'marker':'D'})
        ax.plot(basis, velocity, label=label, **style)
        ax.set_ylabel(r'$\langle \mathrm{velocity} \rangle \mu m/s$')
        # ax.set_xscale('log')
    idx_order = [0, 1]
    idx_order = [1, 0]
    dc.plotmethod(plt.gca(), _meanvel, idx_order=idx_order)
    plt.tight_layout()

# mean/median velocities
@defaultsave()
def localvel(idx_bounds=None):
    # plotutils.default_style()
    dc = DataCube()
    if idx_bounds != None:
        dc.set_idx_bounds(idx_bounds)
    #
    check_requirements = False
    if check_requirements:
        lldata = dc.load_local()
        def has_required(ld):
            return 'vel' in ld
        satisfied = np.vectorize(has_required)(lldata)
        unsatisfied = np.argwhere(~satisfied)
        if np.any(unsatisfied):
            raise RuntimeError('Unsatisfied Requirements:\n{}'.format(str(unsatisfied)))

    idx_order = [1,0]
    # load FJ statistics
    fjld = stats.load(os.path.join(_fj.dataset, stats.localjson))
    def overlay_exp_data(ax, fjld, rescale_max):
        track_ct_stat = twutils.trim(fjld['vel']['each_'+ct_stat], 0.01)
        #
        res = 100
        ybasis, pdf_shape = plotutils.kde_eval(track_ct_stat, res)
        # rescale to x-data
        pdf_shape *= rescale_max/np.max(pdf_shape)
        # plot on yaxis
        print('x', pdf_shape.min(), pdf_shape.max())
        print('y', ybasis.min(), ybasis.max())
        style = {'color':'grey', 'alpha':0.5, 'linewidth':4, 'marker':None}
        ax.plot(pdf_shape, ybasis, **style)

    for ct_stat in ['mean', 'median']:
        #
        idx_order = [0,1]
        fig = plt.figure()
        fig._plot_name = ct_stat
        ctkw = {'field':'lvel.'+ct_stat}
        dc.plotmethod(fig.gca(), _localparam, idx_order=idx_order, lstyle=None, pkw=ctkw)
        
        if dc.dimension > 1:
            #
            idx_order = [1,0]
            fig = plt.figure()
            fig._plot_name = '_'.join([ct_stat, 'swapped'])
            ax = fig.gca()
            dc.plotmethod(ax, _localparam, idx_order=idx_order, lstyle=None, pkw=ctkw)
            #
            rescale_max = np.max(dc.basis[idx_order[0]])/3.
            overlay_exp_data(ax, fjld, rescale_max)


def local1d(field):
    dc = DataCube()
    ctkw = {'field':field}
    ax = plt.gca()
    dc.plotmethod(ax, _localparam, pkw=ctkw)
    command.saveplt(field)


def _localparam(ax, dc, param_idx, ldata, label='', linekw={}, **kw):
    # 1d
    parameter = dc.pnames[param_idx]
    field = kw['field'].split('.')

    if len(field) == 1:
        param, ct_statistic = field[0], None
        linedata = var_or_nan(ldata, lambda ld: ld[param])
    if len(field) == 2:
        param, ct_statistic = field
        linedata = var_or_nan(ldata, lambda ld: ld[param][ct_statistic])
        std_err = var_or_nan(ldata, lambda ld: ld[param]['std_error'])
    basis = dc.slice_basis[param_idx]
    #
    style = {'marker':'D', 'markersize':4}
    this_linestyle = next(linekw)
    style.update(this_linestyle)
    ax.plot(basis, linedata, label=label, **style)
    if ct_statistic == 'mean':
        ax.errorbar(basis, linedata, 1.96*std_err, ls='none', elinewidth=1, **this_linestyle)
    ylabel = kw['field']
    ax.set_ylabel(txtdata.prettynames.get(ylabel, ylabel))
    # TODO: log scale is not parameter agnostic. How do we provide this information?
    if uselog.get(parameter, False):
        ax.set_xscale('log')
        ax.set_xlim(xmin=0.9*10**0)
    xlabel = txtdata.longnames.get(parameter, parameter)
    ax.set_xlabel(xlabel)
    plt.tight_layout()

@defaultsave(autoshow=False)
def rnbt(idx_bounds=None):
    # plotutils.default_style()
    mpl.rcParams["legend.fontsize"] = "medium"
    dc = DataCube()
    dc.set_idx_bounds(idx_bounds)
    fig = plt.figure()
    # order = [1,0] if dc.dimension > 1 else [0,1]
    order = [0,1] 
    dc.plotmethod(fig.gca(), _rnbt, idx_order=order, lstyle=None)
    plt.tight_layout()

def _rnbt(ax, dc, param_idx, ldata, label='', linekw={}, **kw):
    parameter = dc.pnames[param_idx]
    nbound = var_or_nan(ldata, lambda ld: ld['nbound']['mean'])
    nbound_std_err = var_or_nan(ldata, lambda ld: ld['nbound']['std_error'])
    ntaut = var_or_nan(ldata, lambda ld: ld['ntaut']['mean'])
    ntaut_std_err = var_or_nan(ldata, lambda ld: ld['ntaut']['std_error'])
    #
    basis = dc.slice_basis[param_idx]
    style = next(linekw)
    b_style = {'marker': 's', 'linestyle':'-', 'linewidth':1}
    t_style = {'marker': 'x', 'linestyle':'-', 'linewidth':1}
    b_style.update(style)
    t_style.update(style)
    init = 0
    nb_label = r'$\langle N_\mathrm{bound}\rangle$'
    nt_label = r'$\langle N_\mathrm{taut}\rangle$'
    b_label = ' '.join([label, nb_label])
    t_label = ' '.join([label, nt_label])
    ax.errorbar(basis[init:], nbound[init:], 1.96*nbound_std_err, label=b_label , **b_style)
    ax.errorbar(basis[init:], ntaut[init:], 1.96*ntaut_std_err, label=t_label, **t_style)
    ax.set_ylabel(r'$\langle N\rangle$')
    ax.set_xlabel(txtdata.longnames.get(parameter, parameter))
    if uselog.get(parameter, False):
        ax.set_xscale('log')
        ax.set_xlim(xmin=0.9*10**0)


#################################################################################
# inspect timing and data usage
import bookkeeping

# for reference bookkeeping.summary() sets 
# total_clocktime, total_simtime, total_data_size, ptrdatasize, etrdatasize, trdatasize

@defaultsave()
def clocktime():
    # use line graph or image for this?
    dc = DataCube()
    ax = plt.gca()
    def _clocktime(ax, dc, param_idx, ldata, label='', linekw={}, **kw):
        linevals = [ ld['total_clocktime'] for ld in ldata ]
        lstyle = {'marker': 'D'}
        ax.plot(linevals, **lstyle)
        ax.set_ylabel('total clocktime (s)')

    idx_order = [1,0]
    lldata = dc.load_local()
    dc.plotmethod(ax, _clocktime, lldata, idx_order)
    plt.gcf().tight_layout()

# def plotmethod(self, ax, pmethod, lldata=None, idx_order=[0,1], lstyle=None, pkw={}):
#     pmethod(ax, self, 0, lldata, linekw=style, **pkw)


#################################################################################
# attachment distribution on plane


import attachment

@defaultsave()
def att_plane_distrib():
    dc = DataCube()
    def read_att():
        return readtrack.Track('attachment.dat')
    attcube = dc.autocalculate(read_att, force=False)
    pcube = dc.calculate(parameters.thisread)

    pdata = parameters.thisread().params[dc.pnames[0]]
    bigax, linax = dc.lingrid(attcube, attachment._plane,
        big_xlabel=pdata.symbol, pcube=pcube, sharey=True)
    bigax.xaxis.label.set_size(80)
    plt.tight_layout()



#################################################################################
# FJ reorientation distribution on 3d data on a grid

def rmvalues(dc, fastpercent):
    import _analysefj, matdef
    def rmodes():
        Ftracking = readmat.this_fjtrack()
        return _analysefj.frmodes(Ftracking, fastpercent, col=matdef.LEAD)
    valuecube = dc.autocalculate(rmodes, dtype=object, force=False)
    return valuecube

#define a datacube method which lays out n1 * (n2 * n3) plots on a grid
import _analysefj
@defaultsave()
def grid_rmodes(fastpercent=99):
    dc = DataCube()

    valuecube = rmvalues(dc,fastpercent)
    # tmp
    #if valuecube.ndim == 2:
        #valuecube = valuecube[:,0]
    def pfunction(*args):
        return _analysefj.polar(*args, leg_fontsize=10)

    if valuecube.ndim > 1:
        dc.gridspec(valuecube, pfunction, prj='polar')
    else:
        dc.lingrid(valuecube, pfunction, prj='polar')


@defaultsave()
def grid_pili_anchors():
    dc = DataCube()
    mapped_pts = dc.autocalculate(eventanalyse.compute_pili_anchors)
    dc.gridspec(mapped_pts, eventanalyse.pili_anchors)


# helper methods for comparing metrics on a grid
def get_colorcube(dc, frcube, cmap, metric):
    ksdcube, flatidx = dc.sort_by(frcube, metric)
    ksdr = np.percentile(ksdcube.ravel(), 50)
    norm = mpl.colors.Normalize(vmin=ksdcube.min()-(0.05*ksdr), vmax=ksdr, clip=False)
    #cmcube = (ksdcube - ksdcube.min())/(mcut_color-ksdcube.min()) #rescale
    cmcube = norm(ksdcube)
    colorcube = dc.apply(cmap, cmcube, dtype=object)
    # the ScalarMappable object handles normalisation followed by colormapping
    scm = cmx.ScalarMappable(cmap=cmap, norm=norm)
    scm._A = []
    scm.set_clim(ksdcube.min(), ksdcube.max())
    return scm, colorcube

@defaultsave()
def rangle(fastpercent=99):
    dc = DataCube()
    def fr(st):
        return yrstats.frtheta(st, fastpercent)
    # make this autoloading
    frcube = dc.apply(fr, get_exptrack_cube(), dtype=object) # trying to store arrays

    # get metrics
    def get_metrics():
        yt = readyr.inithres()
        yr_theta = yrstats.frtheta(yt, fastpercent)
        # make the basis space independing of the data by setting xlims
        thout = plotutils.kdeplot(yr_theta, xlims=(0,180))
        thspace, yrpde = thout['space'], thout['pde']

        def metric1(sim_theta):
            ksd, pv = scipy.stats.ks_2samp(sim_theta, yr_theta)
            return ksd
        def metric2(sim_theta):
            simout = plotutils.kdeplot(sim_theta, xlims=(0, 180))
            simpde = simout['pde']
            return 180 * simstats.wlsq(yrpde, simpde)
        return {'ksd':metric1, 'wlsq':metric2}
    metrics = get_metrics()

    # for coloring
    cmap = cmx.get_cmap('plasma_r')
    #carr = np.linspace(0.2,1,100)
    #lscmap = mpl.colors.ListedColormap(cmap(carr))
    lscmap = cmap
    metricname = 'ksd'
    scm, colorcube = get_colorcube(dc, frcube, lscmap, metrics[metricname])

    def pfunction(ax, angles, linekw={}, args=None):
        plotutils.ax_kdeplot(ax, angles, linekw)
        ax.set_yticklabels([])
        ax.set_xticklabels([0,90])
        ax.set_xticks([0,90])
        return ax

    #meta = {'fastpercent':fastpercent}
    dc.gridspec(frcube, pfunction, ccube=colorcube, scm=scm)
    plt.suptitle("\"Reorientations\" Fast {}%".format(100-fastpercent), fontsize=40)

    plt.savefig('plots/rangle/rangle_color_by_{}.png'.format(metricname))

@defaultsave()
def best_stau():
    dc = DataCube()
    yt = readyr.inithres()
    yr_step_d = 0.12
    step_d = 0.4

    yrmeanvel = _analysefj.meanvel(yt)
    yrdist = fjstats.calc_timedist(yt, step_d)/step_d * yrmeanvel

    def calc_timedist(st):
        scaling = _analysefj.meanvel(st)
        return fjstats.calc_timedist(st, step_d)/step_d * scaling
    stau_cube = dc.apply(calc_timedist, get_exptrack_cube(), dtype=object)

    # setup metric 
    def metric(simdist):
        ksd, pv = scipy.stats.ks_2samp(yrdist, simdist)
        return ksd
    metric.name = 'stau_ks_metric'

    cmap = cmx.get_cmap('plasma_r')
    scm, colorcube = get_colorcube(dc, stau_cube, cmap, metric)
    #ksdcube, flatidx = dc.sort_by(stau_cube, metric)

    arbcut = 98
    def pfunction(ax, stau, linekw={}, args=None):
        vcut  = np.percentile(stau, arbcut)
        stau = stau[stau < vcut]
        plotutils.ax_kdeplot(ax, stau, linekw)
        ax.set_yticklabels([])
        return ax

    dc.gridspec(stau_cube, pfunction, ccube=colorcube, scm=scm) 

    plt.suptitle("Step Time distribution d = {}".format(step_d), fontsize=40)

# order parameters by ks similarity of R angle distributions and YR data 
def best_rangle(fastpercent):
    dc = DataCube()
    yt = readyr.inithres()
    yr_theta = yrstats.frtheta(yt, fastpercent)

    def fr(st):
        return yrstats.frtheta(st, fastpercent)
    frcube = dc.apply(fr, get_exptrack_cube(), dtype=object) # trying to store arrays

    # ks metric
    def metric1(sim_theta):
        ksd, pv = scipy.stats.ks_2samp(sim_theta, yr_theta)
        return ksd
    metric1.name = 'rangle_ks_metric'

    # make the basis space independing of the data by setting xlims
    thout = plotutils.kdeplot(yr_theta, xlims=(0,180))
    thspace, yrpde = thout['space'], thout['pde']
    # least squares metric
    def metric2(sim_theta):
        simout = plotutils.kdeplot(sim_theta, xlims=(0, 180))
        simpde = simout['pde']
        return 180 * simstats.wlsq(yrpde, simpde)
    metric2.name = 'rangle_wlsq_metric'

    metric = metric1

    ksdcube, flatidx = dc.sort_by(frcube, metric)
    
    cmap = cmx.get_cmap('plasma')
    cmcube = (ksdcube - ksdcube.min())/(ksdcube.max()-ksdcube.min()) #rescale
    colorcube = dc.apply(cmap, cmcube, dtype=object)

    #cutoff = np.percentile(ksdcube, 30)
    #closeidx, closeparams = dc.closeparams(ksdcube, flatidx, cutoff)
    #for i, idx in enumerate(closeidx):
        #print closeparams[i], ksdcube[idx]





# sort by avg number of taut pili
@defaultsave()
def sortrmodes(fastpercent=99):
    ntcube = dc.autoloading(eventanalyse.avgnt, force=False)
    #sdirs = dc.sortdir(ntcube)
    valuecube = rmvalues(fastpercent)
    def pfunction(*args, **kw):
        return _analysefj.polar(*args, leg_fontsize=10, **kw)
    sortby = np.argsort(ntcube.flatten())
    dc.lingrid(valuecube, pfunction, prj='polar', 
            sortby=sortby, svalues=ntcube.flatten())



#################################################################################
# basic numeric analysis
def allft():
    """
    Just print out which parameter sets didn't complete all their simulation runs.
    """
    dc = DataCube()
    #dc.set_mini()
    method = astrack.list_ftime

    lftc = dc.autocalculate(method, force=True)
    i_incomplete = lftc < 1.
    inc_v = lftc[i_incomplete]
    inc_par = dc.pcube[i_incomplete]
    par_v = list(zip(inc_par, inc_v))
    def pretty_pairs():
        for par, v in par_v:
            print(('{: <30s} {}'.format(str(par), v)))
    print((list(dc.pbounds.keys())))
    pretty_pairs()

#################################################################################

# get numpy array where each element is a list of Track objects
def get_Track_cube(force=False):
    dc = DataCube()
    trcube = dc.autocalculate(readtrack.trackset, force=force)
    return trcube

def get_pTrack_cube(force=False):
    dc = DataCube()
    trcube = dc.autocalculate(readtrack.eventset, force=force)
    return trcube

def get_exptrack_cube(force=False):
    dc = DataCube()
    trcube = dc.autocalculate(readmat.this_fjtrack, force=force)
    return trcube

@defaultsave()
def lt_cube():
    dc = DataCube()
    trcube = get_Track_cube()
    dc.gridspec(trcube, shapeplot.apply_surface_annotation()(shapeplot.longtracks))

@defaultsave()
def a_shapelt():
    dc = DataCube()
    trcube = get_Track_cube(force=False)
    dc.gridspec(trcube, shapeplot.ltdraw)

#  primary observable measures

# velocity distributions

# linearised velocity distribution
# @defaultsave()
def rvdistrib(ymax=None):
    dc = DataCube()
    assert(dc.dimension == 1)
    trcube = get_Track_cube(force=True)
    cm = plt.get_cmap('gist_ncar')
    N = trcube.size
    auxfig, auxaxes = plt.subplots(N,sharex=True,figsize=(8,5*N))

    def _ltrvel(ax, linvel):
        h = None
        if (linvel.size > 200):
            try:
                outd = plotutils.ax_kdeplot(ax, linvel, res=50)
                h = outd['handle']
            except:
                print(linvel)
                raise 
        else:
            pass # warn
        return h

    def _indiv(ax, linvel):
        outd = plotutils.ax_kdeplot(ax, linvel, res=50, hist=True)
        return outd['handle']

    fig, ax = plt.subplots()
    if ymax:
        ax.set_ylim((0,ymax))
    ax.set_prop_cycle(color=[cm(1.*i/N) for i in range(N)])
    ax.set_ylabel('P')
    ax.set_xlabel(r'velocity $(\mu m s^{-1})$')
    handles = []
    labels = []
    base = dc.basis[0]
    for i, trs in enumerate(trcube):
        ltrs = [_fj.linearize(tr) for tr in trs]
        linvel = np.concatenate([np.linalg.norm(ltr.get_step_velocity(),axis=1) for ltr in ltrs])
        h = _ltrvel(ax, linvel)
        aux = auxaxes[i]
        aux.xaxis.set_tick_params(which='both', labelbottom=True)
        ih = _indiv(aux, linvel)
        aux.legend([ih],[base[i]])
        if h:
            handles.append(h)
            labels.append(base[i])
    ax.legend(handles, labels)
    fig.tight_layout()
    command.savefig(auxfig, 'rvdistrib_indiv.png')
    command.savefig(fig, 'rvdistrib.png')


##################################
# periodgram

# def rspectrum():
#     pass


def rradgy():
    dc = DataCube()
    dc.fmethod(astrack.radgy,
            get_Track_cube(),
            ylabel='$R_g$',
            force=False)

def ravgvel(step=1):
    dc = DataCube()
    def a_avgvel(trs):
        return twanalyse._avgvel(trs, step)
    args = parameters.thisread()
    dt = step * args.system.deltat
    dc.fmethod(a_avgvel,
            get_Track_cube(force=False),
            ylabel=r'$\delta t = {}$'.format(dt),
            nameex='dt={}'.format(dt),
            force=True)

def rtilt():
    dc = DataCube()
    dc.fmethod(astrack.avgtilt,
            get_Track_cube(),
            ylabel='tilt angle',
            force=True
            )

def rncontact():
    dc = DataCube()
    dc.fmethod(astrack.avgncontact,
            get_Track_cube(),
            ylabel='No. Contacts',
            force=True
            )

# divider dimension 
def div_dimension():
    dc = DataCube()
    dc.fmethod(fjstats._r_div_dimension,
            get_exptrack_cube(),
            force=True)


##

def rksdisp():
    """ 
    1. compute the K-S measure as compared to an experimental distribution 
    2. find the best fit parameters
    """
    print('loading experimental displacements')
    #expdisps = np.loadtxt(os.path.join(twanalyse.analysis_dir, 'fjdata/disps.npy'))
    expvels = np.loadtxt(yrstats.load_vel())
    def ks_measure(trs):
        measure = twanalyse.ks_measure(expvels, trs)
        return measure.statistic
    dc = DataCube()
    #kscube = dc.autoloading(ks_measure,
            #get_Track_cube(),
            #force=True)
    
    kscube = dc.fmethod(ks_measure,
            get_Track_cube(),
            force=True)
    

# compare using mean kmsd
# use a simple comparison for simulations based on a cutoff kmsd value

def bestkmsd():
    dc = DataCube()
    kmsdcube = dc.autocalculate(astrack._kmsd)
    def metric(kmsd):
        return abs(kmsd-readyr.meankmsd)
    metric.name = 'kmsd_metric'
    metriccube,  flatidx = dc.sort_by(kmsdcube, metric)
    print('parameters')
    cutoff = 0.1 # for studying close values of the metric
    closeidx = [dc.ibasis[i] for i in flatidx if metriccube.ravel()[i] < cutoff]
    closeparams = np.array([[dc.basis[j][k] for j,k in enumerate(idx)] for idx in closeidx])
    for params in closeparams:
        print(params)
    print(('means', np.mean(closeparams, axis=0)))
    print(('medians', np.median(closeparams, axis=0)))

    # Avg No. bound
    print()
    print("Mean no. bound pili for good kmsd")
    nboundcube = dc.autocalculate(eventanalyse.avgnb)
    closenbound = [nboundcube[idx] for idx in closeidx]
    print(closenbound)
    print(('mean', np.mean(closenbound)))

    # Avg. No. taut
    print()
    print("Mean no. taut pili")
    nboundcube = dc.autocalculate(eventanalyse.avgnt)
    closenbound = [nboundcube[idx] for idx in closeidx]
    print(closenbound)
    print(('mean', np.mean(closenbound)))



# compare the velocity distributions at the minimum of ks_measure
@command.defaultsave()
def best_compare():
    dc = DataCube()
    expdispline = np.loadtxt(os.path.join(twanalyse.analysis_dir, 
        'yrplots/data/view_vel.dat'))
    exp_space, exp_line = expdispline.T

    fig, axes = plt.subplots(1,2, figsize=(15, 6))
    ax1, ax2 = axes

    ax1.plot(exp_space, exp_line, label='experiment')

    #tau_00.500_pilivar_02.00000_free_tau_eovers_0000.500/plots/data/
    data = dc.reload_data('vels.dat')
    space, line = data[2,0,1].T # reading the best distribution from ks statistic
    ax1.plot(space, line, label='simulated')
    ax1.legend()

    plot_cdf(exp_space, exp_line, ax2, label='experiment')
    plot_cdf(space, line, ax2, label='simulated')

    ax1.set_ylabel('probability density')
    ax2.set_ylabel('cumulative probability density')
    ax1.set_xlabel('frame velocity at 0.1 seconds')
    plt.tight_layout()

def plot_cdf(space, line, ax, **linekw):
    cdf = np.cumsum(line)
    cdf /= cdf[-1]
    ax.plot(space, cdf, **linekw)

def rstats():
    dc = DataCube()
    statscube = dc.autoloading(twanalyse.getstats, get_Track_cube())

#  pili basic measures

def ravgnb():
    dc = DataCube()
    dc.fmethod(eventanalyse.avgnb, get_pTrack_cube(), force=False)

def ravgnt():
    dc = DataCube()
    dc.fmethod(eventanalyse.avgnt, get_pTrack_cube(), force=False)

def rplengths():
    dc = DataCube()
    dc.fmethod(astrack._plengths, get_Track_cube(), force=False)

# generalise this concept for plotting subplots using and registered function

def lin_ravgnb():
    dc.set_trcube(get_pTrack_cube())
    method = eventanalyse.avgnb
    cube = dc.autoloading(method, force=True)
    dc.linsubplots(cube, cube_name=method.__name__)
    sf = os.path.join(command.pdir, method.__name__+'.png')
    plt.savefig(sf)

def lin_ravgnt():
    dc.set_trcube(get_pTrack_cube())
    method = eventanalyse.avgnt
    cube = dc.autoloading(method, force=True)
    dc.linsubplots(cube, cube_name=method.__name__)
    sf = os.path.join(command.pdir, method.__name__+'.png')
    plt.savefig(sf)


###############################################################################
# surface contacts

@defaultsave()
def lin_avgcontacts():
    dc = DataCube()
    def thismethod():
        mean, _ = twanalyse.anymeanvar('ncontacts')
        return mean
    cube = dc.calculate(thismethod)
    ax = plt.gca()
    ax.plot(dc.basis[0], cube, marker="D", linestyle=None)
    ax.set_xlabel(dc.pnames[0])
    ax.set_ylabel("Average No. contacts")

def final_time(tr):
    return tr.track['time'][-1]
def avg_final_time(trs):
    return np.mean(list(map(final_time, trs)))

#################################################################################
#  
import scipy
import scipy.optimize
import scipy.interpolate
import tmos.surface as surface

@defaultsave()
def to_walking_time():
    def make_get_max_force(a, b):
        def get_max_force(eps):
            sp = surface.Spotential(eps, a, b)
            faeps = lambda r: sp.sforce(r) # gradient
            result = scipy.optimize.fminbound(faeps, a, b)
            maxf = abs(faeps(result))
            return maxf
        return get_max_force

    def construct_maxf_functions(get_max_force, basislims=(0.1, 510)):
        emin, emax = basislims
        eps_space = np.linspace(emin, emax, 1000)
        xmaxf = list(map(get_max_force, eps_space))
        def wrap_to_single_value(f): 
            def wrapf(arg):
                return float(f(arg))
            return wrapf
        # forward function
        forward = wrap_to_single_value(scipy.interpolate.interp1d(eps_space, xmaxf))
        # reverse function
        reverse = wrap_to_single_value(scipy.interpolate.interp1d(xmaxf, eps_space))
        return forward, reverse

    def construct_maxf(args):
        b = args.cell.R
        a = b - args.Cell3d.attract_width
        return construct_maxf_functions(make_get_max_force(a, b))

    args = parameters.thisread()
    R = args.cell.R
    get_maxf = make_get_max_force(R-args.pget('attract_width'), R)
    dc = DataCube()
    
    def thismethod():
        return avg_final_time(readtrack.trackset())

    ## setup to output the data
    outd = OrderedDict()
    outd['epsilon'] = dc.basis[0]

    #cube = dc.calculate(thismethod)
    cube = dc.autoloading(thismethod, force=False)
    # 1st dimension is eps
    # 2nd dimension is npili
    f_lines = np.split(cube, cube.shape[1], axis=1)
    basis = list(map(get_maxf, dc.basis[0]))
    outd['maxf'] = basis
    ax = plt.gca()
    lines =  []
    for i, f_line in list(enumerate(f_lines)):
        npili = dc.basis[1][i]
        line = ax.plot(basis, f_line, label='{:d}'.format(int(npili)), marker='D')
        lines.append(line)
        outd['npili={}'.format(npili)] = f_line

    # vertical guide lines
    conf = {'linestyle':'--', 'alpha':0.6}
    lc = (i[0].get_color() for i in lines)
    ax.axvline(25, c=next(lc), **conf)
    ax.axvline(50, c=next(lc), **conf)
    ax.axvline(75, c=next(lc), **conf)
    ax.axvline(100, c=next(lc), **conf)

    ax.set_ylabel("Switching time (s)")
    #ax.set_xlabel("Depth of Potential $\epsilon$ ($pN \mu m$)")
    ax.set_xlabel("Force barrier ($pN$)")

    ax.legend(title='No. Pili', frameon=False)
    plt.tight_layout()

    return outd

@defaultsave()
def to_crawling_time():
    args = parameters.thisread()
    dc = DataCube()

    def thismethod():
        return avg_final_time(readtrack.trackset())
    cube = dc.autoloading(thismethod, force=False)

    ax = plt.gca()
    conf = {'linestyle':'--', 'alpha':0.6}
    ax.axvline(0.5*np.pi/2., **conf)
    ax.axvline(0.5*np.pi/2. + 1., c='k', **conf)
    ax.plot(dc.basis[0], cube, marker='D')

    ax.set_ylim(0, 1.1*args.system.simtime)
    ax.set_ylabel('Switching Time (s)')
    ax.set_xlabel('Pili Distribution Width ($\mu m$)')
    plt.tight_layout()


def stime():
    dc = DataCube()
    dc.apply(simstats._steptime, get_exptrack_cube(), mvdir=True)

#### test spherical geometry

# compute displacements

import twsphere
def main_sphere(modk=1):
    dc = DataCube()
    cube = dc.apply(twsphere.msdlike, get_Track_cube(force=False))
    print((cube.shape))

    #2d array, sph
    radii = [5., 10., 20., 50., 100., 1000.]
    pvar = [0.1, 0.5, 2.0, 7.0]
    curvatures = modk * np.array([1./r for r in radii])
    style = {'marker':'o', 'linestyle':'None'}
    #label='pvar = {}'.format(pv),
    plt.plot(curvatures[:3], cube[:3,:].mean(axis=1),  **style)
    print(('linedata', curvatures[:3], cube[:3,:].mean(axis=1)))
    plt.legend()

    plt.xlabel('curvature (1/radius) (\mu m^-1)')
    plt.ylabel('mean msd after 1000 seconds')

    plt.show()


#### Hexagonal grid geometry

import twgrid
def visitmap(force=True):
    dc = DataCube()
    print('compute/loading visitmap')
    dstcube = dc.autoloading(twgrid._visitmap, 
            get_Track_cube(force=False),
            force=force)
    # saves plot locally
    print('plotting')
    dc.apply(twgrid.hexvisitmap, dstcube, mvdir=True)


#### Using Imshow 

# plot any parameter

## MAIN FUNCTION ##
def summary():
    image_summary()
    jobs = ['lvel.mean', 'qhat.estimate', 'ahat.estimate']
    for par in jobs:
        print('attempt to plot local variable', par)
        param_image(par)

# pass a string like "lvel.mean" and return a getter function
make_get = twutils.make_get

# tell param_image how to format a data field
anform = { 
    'lvel.mean' : '{:6.4f}',
    'lvel.std' : '{:6.4f}',
    'lvel.median' : '{:7.5f}',
    'qhat.estimate' : '{:5.3f}',
    'ahat.estimate' : '{:5.3f}'
}

def param_image(localvar):
    dc = DataCube()
    plt.style.use(plotutils.get_style('image'))
    jd = stats.load()
    getter = make_get(localvar)
    ax = plt.gca()
    an_form = anform.get(localvar, '{:5.2f}')
    _param_image(ax, dc, getter, c_label=localvar, annotate=True, annotate_form=an_form)

    pdir = os.path.join('plots/image', localvar +  ".png")
    print('saving to ', pdir)
    plt.savefig(pdir)
    plt.clf()

def image_summary():
    dc = DataCube()
    # plot a few useful simulation parameters in image format
    plt.style.use(plotutils.get_style('image'))

    def _wrap_angle(s):
        return ' '.join([r'$\langle$',s,r'$\rangle$'])

    extra = ['nbound', 'ntaut', 'bound_pili_participation']
    c_labels = [_wrap_angle('nbound'), _wrap_angle('ntaut'), 'bound_pili_participation'.replace('_', ' ')]
    getter = [
        lambda ld: ld['nbound']['mean'],
        lambda ld: ld['ntaut']['mean'],
        lambda ld: ld['bound_pili_participation']
    ]
    space = 'index'
    for i, c_label in enumerate(c_labels):
        fig = plt.figure()
        ax = fig.gca()
        _param_image(ax, dc, getter[i], c_label=c_label, annotate=True, space=space)
        pdir = {'index':'plots/image/', 'data':'plots/d_image'}.get(space)
        filesystem.safemkdir(pdir)
        name = os.path.join(pdir, extra[i])
        out = '{}.png'.format(name)
        print('saving to ', out)
        plt.savefig(out)
        plt.clf()

def pivot_image_summary():
    dc = DataCube()
    plt.style.use(plotutils.get_style('image'))
    getter = [
        lambda ld: ld['pivot_freq'],
        lambda ld: ld['median_action_trail_dx'],
        lambda ld: ld['median_pivot_trail_dx']
    ]
    extra = ['pivot_freq', 
        'action_trail_dx', 'pivot_trail_dx']
    c_labels = [r'pivot frequency', 
        r'action trail dx', r'pivot trail dx']
    for i, c_label in enumerate(c_labels):
        fig = plt.figure()
        ax = fig.gca()
        _param_image(ax, dc, getter[i], c_label=c_label, annotate=True, 
            annotate_form='{:5.4f}', use_lognorm=False)
        pdir = 'plots/image/'
        filesystem.safemkdir(pdir)
        name = os.path.join(pdir, extra[i])
        out = '{}.png'.format(name)
        print('saving to ', out)
        plt.savefig(out)


# imshow parameter space
# -> ticks (yes)
# -> annotations (yes)
# -> spacing (tricky to make this look good)
# -> colormap (yes)
# -> TODO: annotate confidence interval 

def _param_image(ax, dc, _getter, c_label=None, 
        annotate = False, annotate_form='{:5.2f}', space='index',
        use_lognorm=True, _getter_an=None
        ):
    # plt.style.use(plotutils.get_style('image'))
    if _getter_an is None:
        _getter_an = _getter
    lldata = dc.load_local()
    # iterate over the cube
    C = dc.get_local_array(_getter)
    C_an = dc.get_local_array(_getter_an)
    _data_image(ax, dc, C, c_label=c_label, annotate=annotate,
        annotate_form=annotate_form, space=space, use_lognorm=use_lognorm, C_an=C_an)

viridis = mpl.cm.get_cmap('viridis', 256)
colors = [viridis(x) for x in np.linspace(0, 0.9, 256)]
vir_no_yellow = mpl.colors.LinearSegmentedColormap.from_list('vir_no_yellow', colors)

def _data_image(ax, dc, C, c_label=None, 
        annotate = False, annotate_form='{:5.2f}', space='index',
        use_lognorm=True, C_an=None
        ):
    cmap = vir_no_yellow
    assert(dc.dimension == 2) # only makes sense for two dimensional data
    if C_an is None:
        C_an = C
    X, Y = np.meshgrid(*dc.slice_basis, indexing='ij')
    #
    if space == 'index':
        ax.set_xticks(np.arange(len(dc.slice_basis[1])))
        ax.set_yticks(np.arange(len(dc.slice_basis[0])))
    elif space == 'data':
        ax.set_xticks(dc.slice_basis[1])
        ax.set_yticks(dc.slice_basis[0])
    ax.set_xticklabels(dc.slice_basis[1])
    ax.set_yticklabels(dc.slice_basis[0])
    ax.set_xlabel(txtdata.prettynames.get(dc.pnames[1]))
    ax.set_ylabel(txtdata.prettynames.get(dc.pnames[0]))
    #
    
    # Loop over data dimensions and create text annotations.
    if annotate:
        for idx in dc.slice_ibasis:
            i, j = idx
            val = C_an[i,j]
            text = ax.text(j, i, annotate_form.format(val),
                    ha="center", va="center", color="w")
    vmin = np.nanmin(C)
    if use_lognorm:
        if vmin < 1e-10:
            print("WARNING: min value = 0, use linear color ")
            norm = mpl.colors.Normalize(vmin=np.nanmin(C), vmax=np.nanmax(C))
        norm = mpl.colors.LogNorm(vmin=np.nanmin(C), vmax=np.nanmax(C))
    else:
        norm = mpl.colors.Normalize(vmin=np.nanmin(C), vmax=np.nanmax(C))
    if space == 'index':
        Image = ax.imshow(C, norm=norm, cmap=cmap, origin='lower')
    elif space == 'data':
        Image = ax.pcolormesh(Y, X, C, shading='nearest', norm=norm)
    if not (vmin < 1e-10):
        bar = plt.colorbar(Image)
        if c_label:
            bar.set_label(c_label)
    plt.tight_layout()


if __name__=='__main__':
    try:
        ff, thisargs = command.process(locals())
        ff(*thisargs)
    except RuntimeError as e:
        dc = DataCube()
        print()
        print(str(dc))
        raise e
