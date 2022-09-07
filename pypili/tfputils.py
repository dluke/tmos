
"""
Define the standard simulate method for a cell which is not slowed down by 
collecting debugging information.
"""

import sys, os, time
import glob
import logging

import tmos
import tmos.base 
import tmos.pili as pili
import tmos.mdynamics as md
import tmos.surface as surface
if 'vtkwriter' in dir(tmos):
    import tmos.vtkwriter as vtkwriter

# using
from tmos.base import Vector3d

import numpy as np
import scipy
import wrinterval as wr

# needs to be an object not a module
import datapaths


# exceptions
class MDnoconverge(Exception):
    pass

class MDexplode(Exception):
    pass

# give this null method to exit_condition as default 
def no_exit_condition(cell):
    return False

""" for sequencing output """
class Sampler(object):
    def __init__(self, simtime, deltat):
        outt = np.arange(0., simtime+(1.1*deltat), deltat)
        self.outi = 0
        self.outt = iter(outt)
        next(self.outt) # we write out the time 0 configuration elsewhere?
        self.outtarget = next(self.outt)
        self.passed_targets = []

    # it is possible to step forward by a time greater than self.deltat
    # In this case we need to step the iterator multiple times to catch up
    # WORKING : force an output on every step of self.deltat even if nothing changes
    def check(self, tme):
        self.passed_targets = []
        if tme > self.outtarget:
            self.outi += 1
            while True:
                self.passed_targets.append(self.outtarget)
                try:
                    self.outtarget = next(self.outt)
                except StopIteration:
                    break # kmc step jumped over the end of the simulation time which can happen for low pili production
                if self.outtarget > tme:
                    break
            return True
        else:
            return False

    @classmethod
    def from_args(cls, args):
        return cls(args.system.simtime, args.system.deltat)

# read
WLCconfigurations_path = os.path.join(os.path.dirname(__file__), '../single_pilus/CPDfunctions/Cdata_*.npy')
def read_configuration_npy(path=WLCconfigurations_path):
    """
    >>> len(read_configuration_npy())
    20
    """
    Cfiles = sorted( glob.glob(WLCconfigurations_path) )
    return [np.load(cfile).astype('float64') for cfile in Cfiles]

################################################################################
# need to replace parameter handling with a serios objected oriented parameter handling
# system in python

### write parameters for the run to keep track 
from collections import OrderedDict

# utility

# function constructs a simple counting generator
def make_counter(step=1):
    x = 0
    while True:
        yield x
        x += step

# function construction a generator which returns True on x % step == 0 else Fals
def make_interval(step=1):
    x = 0
    while True:
        yield x % step == 0
        x += 1

################################################################################
# debugging output utilities
def get_bound_pili(cell):
    bound = []
    for pilus in cell.pili:
        if pilus.isbound == 1:
            bound.append(pilus.idx)
    return bound
def get_pili_avg_len(cell):
    return sum([pilus.pl for pilus in cell.pili])/12.
def get_pili_lens(cell):
    return [pilus.pl for pilus in cell.pili]

def get_pili_dtheta(cell):
    for pilus in cell.pili:
        angle = pilus.axis.angle(pilus.axisEq)
        print(angle)


####################################################################################
# debugging only

# utility for break_fgrad
from matplotlib import pyplot as plt
def dumpplot(arr, name):
    plt.clf()
    plt.plot(arr, label=name)
    plt.legend()
    plt.savefig(name+'.png')


"""part of debugging loop for my minimisation implementation """
def break_fgrad(args, cell, mdmin):
    mdmin.reset(cell)
    ft = open('touching.txt', "w")  

    arrrms = np.zeros(mdmin.maxsteps)
    arrdt = np.zeros(mdmin.maxsteps)
    arren = np.zeros(mdmin.maxsteps)
    arrcontacts = np.zeros(mdmin.maxsteps)
    arrdx = np.zeros((mdmin.maxsteps, 6))
    arrgrad = np.zeros((mdmin.maxsteps, 6))
    arrstate =  np.zeros((mdmin.maxsteps, 6))
    arrxax =  np.zeros((mdmin.maxsteps, 6))
    rms = 0
    failed = False
    for istep in range(mdmin.maxsteps):
        try: 
            rms = mdmin.one_step(cell)
        except Exception as e:
            failed = True
            break
        arrrms[istep] = rms
        arrdt[istep] = mdmin.dt
        arren[istep] = cell.energy()
        arrcontacts[istep] = cell.get_num_contacts()
        arrdx[istep] = mdmin.get_dx()
        state = cell.get_state()
        grad =  cell.grad()
        arrstate[istep] = state # state vector is [body_center, p_rotation]
        cpt = cell.real_centre()
        axis = cell.body.get_axis()
        arrxax[istep] = np.array([cpt.x, cpt.y, cpt.z, axis.x, axis.y, axis.z])
        arrgrad[istep] = grad
        ft.write(cell.report_touching() + "\n")
        if np.isnan(state).any() or np.isnan(grad).any():
            print(cell)
            sys.exit("found nan")
        if rms < mdmin.target:
            break
    ft.close()

    rms_condition = rms > args.fire.force_target
    step_condition = istep == args.fire.maxsteps -1
    if rms_condition or step_condition or failed:
        print('FAILED CONDITION rms = {}'.format(rms))
        print('steps = {}'.format(istep))
        dumpplot(arrrms, 'rms')
        dumpplot(arrdt, 'dt')
        dumpplot(arren, 'energy')
        dumpplot(arrcontacts, 'ncontacts')
        np.savetxt('arrdx.npy', arrdx)
        np.savetxt('arrstate.npy', arrstate)
        np.savetxt('arrgrad.npy', arrgrad)
        np.savetxt('arrxax.npy', arrxax)
        raise MDnoconverge("\nrms = {}\ntarget = {}\nsteps = {}\ntime = {}".format(
            rms, 'unknown', istep, 'unknown'))
    
    return istep, rms

#################################################################################
# describe the main simulation loop


# string to print on simulation unexpected stop
def write_exit(exit_string):
    with open('final.txt', 'w') as f:
        f.write(exit_string)

dx_threshold = 2 * 0.004
length = 2.0
dot_theshold = np.cos(2 * 0.004/((length+1)/2))
# ...
def md_events(mdev, data):
    #
    mdev.track_mdevent(data)


def _handle_md(args, cell, mdmin, tme):
    savestate = cell.get_state()
    try:
        istep, rms = mdmin.step(cell)
    except Exception as e:
        cell.set_state(savestate)
        istep, rms = break_fgrad(args, cell, mdmin)
        raise e

    # For debugging
    force_condition = rms > args.fire.force_target
    step_condition = istep == args.fire.maxsteps
    nan_condition = np.any(cell.get_state() == np.nan)
    failed_condition = force_condition or step_condition or nan_condition
    if args.system.debug and failed_condition:
        ## can try again with debugging
        print(args)
        print("debugging is ON, rerunning minimisation")
        cell.set_state(savestate)
        istep, rms = break_fgrad(args, cell, mdmin)

    if force_condition or step_condition:

        excp = MDnoconverge(
                "\nforce = {}\ntarget {}\nnsteps {}\ntime {}".format(
            rms, args.fire.force_target, istep, tme))
        write_exit(str(excp))
        raise excp
    if nan_condition:
        excp = MDexplode('Found nan in cell state. Exit\n{}'.format(
            cell.get_state))
        write_exit(str(excp))
        raise excp

    # return the output of mdmin.step()
    return istep, rms


# cell object
# args
# integrator
# fpili is the name form of the pili binding and unbinding file.
# name_out name form of track output file to store in /data 
def simulate_cell(cell, args, mdmin, datfile, exit_condition): 

    simtime = args.system.simtime
    deltat = args.system.deltat
    # setup output tracker
    outputwhen = Sampler.from_args(args)
    vtkwhen = Sampler(args.system.simtime, args.system.vtkdeltat)


    fout = datfile['track'].format(cell.idx)
    bout = datfile['event'].format(cell.idx)
    pout = datfile['pili'].format(cell.idx)
    mdout = datfile['md'].format(cell.idx)

    txy = wr.Track3d.setup_track_3d(fout)
    bind = wr.Track3d.setup_event_3d(bout)
    if args.system.output_pilidata:
        tpili = wr.PiliTrack3d(pout)
        tpili.include_unbound = True
    mdev = wr.Track3d.setup_mdevent(mdout)

    # 
    kmc = pili.Kmc(1./args.Pili.k_adh)

    # timing
    start_t = time.time()

    # setup vtk output
    if args.system.vtkwrite:
        vtkbodypath = datapaths.vtkit_body()
        vtkpilipath = datapaths.vtkit_pili()
        vtkwriter.write_cell3d(cell, next(vtkbodypath))
        vtkwriter.write_pili3d(cell, next(vtkpilipath))

    exit_early = False
    try:
        loop = 0
        tme = 0.
        after_tme = 0.
        while True:

            if args.system.debug:
                #print "time is {}".format(tme)
                pass

            kmcstep =  kmc.kmcstep(cell, tme)
            after_tme, process, need_md, pidx = kmcstep

            # only perform the equilibrate step if process is detach/extend/retract
            istep = 0
            rms = 0.
            if need_md:
                istep, rms = _handle_md(args, cell, mdmin, tme)
                md_events(mdev, (tme, pidx, process, cell, mdmin.get_summary()) )
                if not exit_early and exit_condition(cell):
                    print("Early exit condition satisfied after time", tme)
                    exit_early = True
                kmc.postmd(cell, tme)
            tme = after_tme

            for ev in kmc.events:
                bind.track_event(ev)

            # primarily debugging tool
            if args.system.track: 
                state = wr.Simstate(kmcstep, cell, istep, rms)
                line = txy.track(state)
                #
                if args.system.output_pilidata:
                    tpili.track_pili(tme, cell)
            else:
                # regular text output
                if outputwhen.check(tme):
                    for output_target in outputwhen.passed_targets:
                        # a loop is used to hangle the edge case where kmc update time > args.system.deltat (output time)
                        _kmcstep = (output_target, process, need_md, pidx)
                        state = wr.Simstate(_kmcstep, cell, istep, rms)
                        # body tracking
                        line = txy.track(state)
                        # pili tracking
                        if args.system.output_pilidata:
                            tpili.track_pili(output_target, cell)

            # vtk output for debugging and visualisation
            if args.system.vtkwrite and vtkwhen.check(tme):
                for output_target in vtkwhen.passed_targets:
                    vtkwriter.write_cell3d(cell, next(vtkbodypath))
                    vtkwriter.write_pili3d(cell, next(vtkpilipath))

            # exit condition
            if tme > args.system.simtime or exit_early:
                break

            loop += 1

    except:
        # write error which ended the program to a file 
        with open('final.txt', 'a') as ff:
            ff.write("cell {} raised exception at time = {}\n".format(cell.idx, tme))
            ff.write(repr(sys.exc_info()) + '\n')
        raise
    finally:
        # if not closed we may reach open file limit on cluster
        txy.close()
        bind.close()
        mdev.close()
        if args.system.output_pilidata:
            tpili.close()

        end_t = time.time()
        exec_t = end_t - start_t
        with open("timing.txt", 'a') as ft:
            ft.write("cell {:04d} simtime {:11.4f} clock_time {}\n".format(cell.idx, tme, exec_t))


