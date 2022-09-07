#!/usr/bin/env python

# this is the configuration file.

# Just to reiterate that. This is the configuration file, not part of the code.

# python implementation
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
args.system.simtime = 200
args.system.repeats = 1
args.system.vtkwrite = True

#################################################################################

ks = 10**4.
args.fire.maxsteps = 1000
args.fire.dt = 0.000005
args.fire.dt_max = 10 * args.fire.dt
args.fire.Nmin = 5
#
args.fire.target = 5./ks

#################################################################################

#
args.Pili.ks = ks

#
args.Pili.f_stall = 100.
# or a little less than f_stall, i.e. 0.9 * f_stall
args.Pili.f_release =  args.Pili.f_stall

# elementary step for pili
args.Pili.d_free = 0.004
args.Pili.d_bound = 0.004

# pili force parameter
args.Pili.min_length = 0.0
# pili geometry parameters
args.Pili.inside_length= 10 * args.Pili.d_bound

# pili model
args.cell.pili_model = 'wlc'
args.pili.a = 1.0
#args.pili.ka = -1 # only used in pili configuration simulation
args.Pili.force_max_length = True
args.Pili.max_length = 8.9

pv = 0.5
args.Pili.kb_sh = pv/args.Pili.d_bound
args.Pili.kb_ex = pv/args.Pili.d_bound
args.Pili.kf_sh = pv/args.Pili.d_free;
args.Pili.kf_ex = pv/args.Pili.d_free;

# this arbitrary timescale sets the steepness of the release rate as a function of force
args.Pili.tau = 0.1
# Major timescale controlling step size (switching frequency) for retract/extend dynamics
args.Pili.rtau = 2.

# check attachment 10 times per second
args.Pili.k_adh = 10.

# rtau sets the retraction timescale, but we also need a parameter to set the 
# rate of switching from extension back to retraction 

#// the r_over constant is a stand-in used to describe free pili elongation 
# My hypothesis is that free pili must predominantly extend so free elong_to_sh_rate < 1/rtau
# setting to 1 means that free pili e -> r rate = r -> e rate
# setting simplify_free_pili_etor = False, turns us back to the original system 
# where free pili e->r rate is much less than r -> e rate
r_over = 1.
args.Pili.simplify_free_pili_etor = True
args.Pili.free_elong_to_sh_rate = r_over * 1./args.Pili.rtau


#################################################################################
# number of pili used to construct cells
args.cell.npili= 4

R = 0.5;
length = 2.
args.cell.R = R
args.cell.length = length
#args.ACell.pilivar = np.pi * R/2.
args.ACell.pilivar = np.pi * R/10.

# cell surface interaction
args.Cell3d.eps = 100.
args.Cell3d.repulsive_only = False
args.Cell3d.eps_attract = 1.

#################################################################################
# End Configuration
#################################################################################

# SET RANDOM SEED. ONLY FOR DEBUGGING
import tmos.base as base
base.init_generator(0)

#################################################################################


def single_run():
    args.system.repeats = 1
    args.system.simtime = 500
    singlepool = range_tfp.buildsinglepool(args)

    cores = 10
    range_tfp.runpool(singlepool, args, cores)

# shortened names
snames = [('pilivar', 'pilivar'), ('k_adh', 'kadh'), ('free_elong_to_sh_rate', 'fetosh')]

def pool():
     import range_tfp

    # construct parameter dictionary

    pilivar_setup = ('pilivar', [pi/2., pi/4., pi/6., pi/8., pi/40])
    k_adh_setup = ('kadh', [1/16., 1/8., 1/4., 1/2., 1., 2., 5.])
    r_over = np.array([0.02, 0.1, 0.3, 1.0])
    r_over_setup = ('fetosh', r_over * Pili.sh_to_elong_rate)

    allparams = OrderedDict([pilivar_setup, k_adh_setup, r_over_setup])

    maxcores = 16
    pool = range_tfp.buildpool(args, allparams)
    range_tfp.runpool(pool, args, maxcores)

def main():
    #cell = tfp.initialise_cell_for_step()
    #cell = tfp.initialise_wlc_cell()
    #tfp.minrun(args, cell)

    single_run()
    #pool()


if __name__=='__main__':
    main()

