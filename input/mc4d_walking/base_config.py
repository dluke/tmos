#!/usr/bin/env python
# this is the configuration file.

# Just to reiterate that. This IS the configuration file, not part of the code.

# python implementation
import sys,os
import itertools
import numpy as np

# C++ implementation
import range_tfp
import ctfp3d as tfp

from math import pi, exp

from collections import OrderedDict


#################################################################################
# Construct ParameterList object

# setup defaults
import runtfp_global
args = runtfp_global.setup()


#################################################################################
# Configuration
#################################################################################

# write all frames to trackxy.dat
args.system.track = False

args.system.deltat = 0.1 # frequency of collecting data
args.system.vtkdeltat = 1.0
args.system.simtime = 2000
args.system.repeats = 1
args.system.vtkwrite = False
args.system.debug = False 

#################################################################################

# Energy minimisation type

# styles are fgrad, NM, LBFGS
args.system.minstyle = 'myfire'
#args.system.minstyle = 'pelefire'

ks = 10**4.
args.fire.maxsteps = int(1e4)
args.fire.dt = 1e-6  # reasonable for fixed timestep
args.fire.dt_max = 1e-4
args.fire.Nmin = 5
#
# pele
args.fire.force_target = 1.*np.sqrt(6.)
args.fire.f_dec = 0.5
args.fire.f_inc = 1.1
args.fire.maxstep = 0.001

#args.Pili.force_threshold = args.fire.target_target
args.Pili.force_threshold = 1.

#################################################################################

# investigate this parameter for higher velocity peaks
args.Pili.ks = ks

# tug of war hypothesis
args.Pili.f_stall = 100.
args.Pili.enforce_stalling = True

# elementary step for pili
args.Pili.d_free = 0.004
args.Pili.d_bound = 0.004


# pili model
args.cell.pili_model = 'wlc'
args.cell.pili_generator = 'kpgen'
args.pili.a = 0.2
args.Pili.force_max_length = True
args.Pili.max_length = 19.99

# 
free_extension = 0.28 # \mu s^-1
free_retraction = 0.17
bound_retraction = 0.09
args.Pili.kb_sh = bound_retraction/args.Pili.d_bound
args.Pili.kb_ex = free_extension/args.Pili.d_bound
args.Pili.kf_sh = free_retraction/args.Pili.d_free
args.Pili.kf_ex = free_extension/args.Pili.d_free


# 
args.Pili.k_ext_off = 0.625


# limit extension of pili (no bending)
args.Pili.allow_bound_extension = False
# Tala-like model (i)
# assume that surface detection forces immediate unbinding of ext_motor 
# and prevents rebinding while bound
args.Pili.allow_bound_ext_motor = False
args.Pili.force_bound_retraction = False

# increase dwell time approx 1s (Tala)
# based on force sensing hypothesis 
# (this is the next parameter we should consider investigating)
args.Pili.dwell_time = 1.0
# demand that pili detach on time (default)
args.Pili.force_detachment = True

# 
args.Pili.detach_grace_length = 0.
args.Pili.detach_grace_time = 0.

#################################################################################
# number of pili used to construct cells
args.ACell.k_spawn = 2.0
args.cell.npili = 0

R = 0.5;
length = 2.
args.cell.R = R
args.cell.length = length

# use guassian distribution
args.Cell3d.distrib_type = 'gaussian'
args.ACell.pilivar = 2.0

# cell surface interaction
# fix in crawling mode
args.Cell3d.eps = 200.
args.Cell3d.repulsive_only = True

# Surface configuration
args.surface.shape = 'plane'

# crawling mode
args.surface.initial_contact = False

################################################################################
# End Configuration
#################################################################################

# SET RANDOM SEED. ONLY FOR DEBUGGING
import tmos.base as base
base.init_generator(0)
# Using the production run loop overrides

#################################################################################
