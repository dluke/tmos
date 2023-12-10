
"""
1. Standard run method for constructing cell object from configuration arguments
+ utilities for constructing cell object

2. Minrun simulation loop which collects heavy debugging formation

"""

import os
import math
import numpy as np
import matplotlib.pyplot as plt

# TMOS setup

import tmos.base
import tmos.surface as surface 
import tmos.pili as pili
import tmos.mdynamics as md
if 'vtkwriter' in dir(tmos):
    import tmos.vtkwriter as vtkwriter

# using
from tmos.base import Vector3d

import txtdata
import wrinterval as wr
import datapaths
import tfputils

from tfputils import MDnoconverge, MDexplode

import time

class Timer(object):
    def __enter__(self):
        self.start_t = time.clock()

    def __exit__(self, type, value, tb):
        self.end_t = time.clock()
        dt = self.end_t - self.start_t
        with open('timing.txt', 'a') as ft:
            ft.write('execution time: {}s\n'.format(dt))


pi = np.pi
# hide/active print statements
verbose = False


def seqtov3d(seq):
    assert(len(seq) == 3), 'variable seq {} should be length 3'.format(seq)
    x, y, z = seq
    return Vector3d(x,y,z)

def v3dtoseq(v):
    #return np.array([v.x,v.y,v.z])
    return [v.x,v.y,v.z]


def setup_kpgen(args):
    import scipy.constants
    kb = scipy.constants.Boltzmann
    T = 273.0 + 30.0 # absolute temperature
    kbT = (kb * T) * 10**18. 
    # Lp = 5. # persistence length
    Lp = args.pili.Lp
    a = args.pili.a
    ka = Lp * kbT
    nmax = int(args.Pili.max_length/a)
    kpgen = pili.KPgenerator(ka, a, nmax, T)
    return kpgen

def setup_wlc_generator(args):
    if args.cell.pili_generator == 'kpgen':
        return setup_kpgen(args)
    elif args.pili.generator == 'sample':
        cls = tfputils.read_configuration_npy()
        # doesn't exist yet
        #args.pili.ka 
        wlcg = pili.WLCgenerator(-1., args.pili.a, len(cls))
        assert(args.Pili.force_max_length)
        assert(args.Pili.max_length < len(cls) * args.pili.a), 'values: {} {} {}'.format( args.Pili.max_length, len(cls), args.pili.a)
        for i, cdata in list(enumerate(cls)):
            wlcg.append(cdata)
        return wlcg

# just return the Equilibrate object initialised with defaults 
def setup_fire(args):
    # add params
    alpha_start = 0.1
    f_alpha = 0.99
    # setup
    mdmin = md.Equilibrate3d(
            args.fire.dt,
            args.fire.maxsteps,
            args.fire.force_target,
            args.fire.Nmin,
            args.fire.f_inc,
            args.fire.f_dec,
            alpha_start,
            f_alpha,
            args.fire.dt_max,
            args.fire.maxstep)
    return mdmin

def setup_md3d(args):
    # default
    return setup_fire(args)

def setup_cell_surface(args):
    # add surface types ...
    if args.surface.shape == 'plane':
        sf = surface.Plane()
        # flipped
        #sf = surface.Plane(Vector3d(), Vector3d(0,0,-1))
    elif args.surface.shape == 'hexspheregrid':
        Normal = Vector3d(0,0,1)
        sf = surface.HexSphereGrid(args.surface.sphere_radius, Vector3d(), Normal)
        # write the surface to the current working directory
        if args.surface.write_surface and not os.path.exists("HexGrid.vtk"):
            print("Writing the surface to file...")
            vtkwriter.write_hgrid(surface.HexGrid(args.surface.sphere_radius)) # file name hardcoded
    elif args.surface.shape == 'none':
        sf = surface.NullPlane()
        print(str(sf))
    elif args.surface.shape == 'infsteps':
        smallr = args.cell.R/2.
        sep = args.surface.sep
        height = args.surface.height
        sf = surface.InfSteps(height, sep, smallr)
        if args.surface.write_surface and not os.path.exists("InfSteps.vtk"):
            print("Writing the surface to file...")
            vtkwriter.write_infsteps(sf) 
    elif args.surface.shape == 'segplane':
        sf = surface.SegPlane(args.surface.planeangle, args.surface.panewidth)
        if args.surface.write_surface and not os.path.exists("SegPlane.vtk"):
            print("Writing the surface to file...")
            vtkwriter.write_segplane(sf)
    elif args.surface.shape == 'sphere':
        radius = args.surface.sphere_radius
        print("Creating spherical surface at origin with radius {}".format(radius))
        sf = surface.Sphere(radius)
    elif args.surface.shape == 'shell':
        radius = args.surface.sphere_radius
        print("Creating spherical shell at origin with radius {}".format(radius))
        sf = surface.Shell(radius)
    elif args.surface.shape == 'sineplane':
        frame = tmos.base.Frame()
        sf = surface.SinePlane(frame, args.surface.A, args.surface.B, args.surface.nsign)
        print("Creating SinePlane at {} with paramters {} {} {}".format(
                frame.get_origin().__str__(), sf.A, sf.B, args.surface.nsign))
        if args.surface.write_surface and not os.path.exists("SinePlane.vtk"):
            print("Writing the surface to file...")
            vtkwriter.write_sineplane(sf)

    return sf

def setup_cell_body(args, sf):
    # cell body connected to surface

    if args.surface.shape == 'hexspheregrid':
        cellaxis, cellcenter = setup_on_hexspheregrid(args, sf, args.surface.setup_style)

    elif args.surface.shape == 'plane': 
        # logic for plane
        delta = 0.01
        if args.surface.initial_contact:
            # for horizontal axis
            cellaxis = Vector3d(1,0,0)
            height = args.cell.R-delta
            cellcenter = Vector3d(0.,0.,height)
        else:
            cellaxis = Vector3d(0,0,-1)
            height = args.cell.length/2. + args.cell.R - delta
            cellcenter = Vector3d(0.,0.,height)
    elif args.surface.shape == 'infsteps':
        # use position set in input file
        cellaxis = seqtov3d(args.cell.initaxis)
        cellcenter = seqtov3d(args.cell.initcenter)
    elif args.surface.shape == 'segplane':
        # use position set in input file
        cellaxis = seqtov3d(args.cell.initaxis)
        cellcenter = seqtov3d(args.cell.initcenter)
    elif args.surface.shape == 'sineplane':
        cellaxis, cellcenter = setup_on_sineplane(args, sf, args.surface.setup_style)
    else:
        cellcenter = seqtov3d(args.cell.initcenter)
        cellaxis = seqtov3d(args.cell.initaxis)
    
    if verbose:
        print("Setting up cell body with \ncenter {}\norientation{}".format(cellcenter, cellaxis))
    # update args so that config.txt is correct
    args.cell.initaxis = v3dtoseq(cellaxis)
    args.cell.initcenter = v3dtoseq(cellcenter)

    body = surface.Capsule(
            cellcenter,
            cellaxis,
            args.cell.R,
            args.cell.length)
    return body

def setup_on_hexspheregrid(args, sf, word):
    # place the cell vertically in contact with the surface
    delta = 0.001
    normal = 1.
    if word == 'vertical':
        height = normal*(args.cell.length/2.+args.cell.R+args.surface.sphere_radius-delta)
        cellcenter = Vector3d(0.,0.,height)
        cellaxis = Vector3d(0,0,-normal)
        return cellaxis, cellcenter
    else:
        print("Warning: no cell body setup style defined on hexspheregrid")
        return args.cell.initaxis, args.cell.initcenter


def setup_on_sineplane(args, sf, word): 
    delta = 0.001
    R = args.cell.R
    length = args.cell.length
    if word == 'contact':
        n = sf.normal_form(0)
        initaxis = Vector3d(0,0,-1)
        initcenter = sf.get_origin() + (R-delta)*n - (length/2)*initaxis
        return initaxis, initcenter
    if word == 'minimum':
        troughat = -pi*sf.B/2.
        center = sf.sp_form(troughat) + Vector3d(0,0,R-delta)
        axis = Vector3d(0,1,0)
        return axis, center
    if word == 'perpcrest':
        crestat = pi*sf.B/2.
        center = sf.sp_form(crestat) + Vector3d(0,0,R-delta)
        axis = Vector3d(1,0,0)
        return axis, center
    else:
        print("Warning: no cell body setup style defined")
        return args.cell.initaxis, args.cell.initcenter


def setup_cell(args, idx=0):
    sf = setup_cell_surface(args)
    body = setup_cell_body(args, sf)
    if args.cell.pili_model == "wlc":
        wlcg = setup_wlc_generator(args)
        cell = pili.CellWLC3d(idx, args.cell.npili, body, sf, wlcg)
    elif args.cell.pili_model == "rod":
        cell = pili.Cell3d(idx, args.cell.npili, body, sf)
    # call once
    cell.common_init()
    return cell

#################################################################################
# mimic ctfp.run 
# and use generic tfputils.simulate_cell mainloop

# Setup and run simulation defined by ParameterList
import time, sys, os
from multiprocessing import Process
def run(args, idx, exit_condition):
    # initialise the random number generator
    # print('using random seed ', args.system.seed)
    args.init_generator(args.system.seed)

    cell = setup_cell(args, idx)

    class MD: pass
    mdmin = MD() # dummy class

    # initial MD step
    if args.system.minstyle == 'myfire':
        # print('Using C++ implementation of Fire')
        mdmin = setup_md3d(args)
    elif args.system.minstyle == 'pelefire':
        # print('Using Pele implementation of Fire')
        mdmin.step = setup_pelefire(args)

    assert(hasattr(mdmin, 'step')) # assert minimiser exists

    # output the parameters before the run starts
    fconf = 'config.txt'
    with open(fconf, 'w') as fc:
        args.write_config(fc)
    fseeds = 'random_seeds.txt'
    with open(fseeds, 'a') as fs:
        fs.write("{:02d}  {:d}\n".format(idx, args.system.seed))

    # setup data directory and output file
    datdir = wr.datdir
    if not os.path.exists(datdir):
        os.mkdir(datdir)
    datfile = {}
    # track output
    datfile['track'] = os.path.join(datdir, 'bacterium_{:05d}.dat')
    # pili binding/unbinding and other events
    datfile['event']= os.path.join(datdir, 'event_{:05d}.dat')
    # pili tracking
    datfile['pili'] = os.path.join(datdir, 'pili_{:05d}.dat')
    # mdstep 
    datfile['md'] = os.path.join(datdir, 'mdevent_{:05d}.dat')


    # vtk output directory
    if args.system.vtkwrite:
        import datapaths
        datapaths.setup_vtkdir()

    from tfputils import simulate_cell
    onecellargs = (cell, args, mdmin, datfile, exit_condition)
    p = Process(target=simulate_cell, args=onecellargs)
    return p

#################################################################################

def run_attachment(args, cell=None, n_attach_stop=None):
    if n_attach_stop < 0: # according to parameters.xml default value is -1
        n_attach_stop = None
    if args.system.debug:
        print('Using <attachment> trial')

    # init_generator for reruns
    tmos.base.init_generator(args.system.seed)

    cell = setup_cell(args)

    fconf = 'config.txt'
    with open(fconf, 'w') as fc:
        args.write_config(fc)

    pilitrackingwhen = tfputils.Sampler.from_args(args)
    vtkwhen = tfputils.Sampler(args.system.simtime, args.system.vtkdeltat)

    # setup vtk output
    datapaths.setup_vtkdir()
    vtkbodypath = datapaths.vtkit_body()
    vtkpilipath = datapaths.vtkit_pili()
    vtkwriter.write_cell3d(cell, next(vtkbodypath))
    vtkwriter.write_pili3d(cell, next(vtkpilipath))

    # 
    kmc = pili.Kmc(1./args.Pili.k_adh)

    # setup tracking
    pout = 'pilitracking.dat'
    att_out = 'attachment.dat'
    tpili = wr.PiliTrack3d(pout)
    tevent = wr.Track3d('event.dat', header=wr.event_3d)
    tattach = wr.Track3d(att_out, header=wr.attachment_event)

    # spawn the first pilus
    pilus = cell.spawn_pilus()
    cell.add_pilus(pilus)
    # create a CellEvent just holding cell data and populate it by hand 
    ev = cell.create_event()
    ev.pidx = pilus.idx
    ev.process = 'None'
    ev.trigger = 'spawn'
    tevent.track_event(ev)

    # only output the inital config
    fout = 'trackxy.dat'
    txy = wr.Track3d.setup_track_3d(fout)
    state = wr.Simstate((0,'None',0,0), cell, -1, -1)
    
    txy.track(state)

    # use one of these numbers to decide when to end the simulation
    n_attach = 0
    tme = 0.
    #
    while True:
        kmcstep =  kmc.kmcstep(cell, tme)
        tme, process, need_md, pidx = kmcstep

        # vtkoutput
        vtkwrite_condition = args.system.vtkwrite and vtkwhen.check(tme)
        if vtkwrite_condition:
            vtkwriter.write_cell3d(cell, next(vtkbodypath))
            vtkwriter.write_pili3d(cell, next(vtkpilipath))

        if args.system.track or pilitrackingwhen.check(tme):
            pilus = cell.pili[0]
            if False:
                print('MC process {} at time {}, leq={}'.format(
                    process, tme, pilus.leq))
            tpili.track_pili(tme, cell)

        for ev in kmc.events:
            if ev and ev.trigger == 'spawn':
                # track spawn events which happen in response to dissolving pili in kmc step
                tevent.track_event(ev)

            if ev and ev.process == 'attach':
                if args.system.debug:
                    print('time {:>12.3f}  n_attach {:>8d}'.format(tme, n_attach))
                tattach.track_attachment(ev)
                n_attach += 1
                # the last thing we do is destroy this pilus and create a fresh one
                cell.dissolve_pilus(pilus)
                pilus = cell.spawn_pilus()
                cell.add_pilus(pilus)
                # if we already tracked an event on this step, track another one
                # have to get an updated event here because we modified the cell
                ev = cell.create_event()
                ev.pidx = pilus.idx
                ev.trigger = 'spawn'
                ev.process = 'attach'
                tevent.track_event(ev)
        
        # kmc.events goes out of scope

        # exit condition
        if tme > args.system.simtime:
            break
        if n_attach_stop and n_attach >= n_attach_stop:
            break

#################################################################################
# debugging main loop


printform = 'performed MC process {} at time {}'.format(
        txtdata.hform['process'], txtdata.hform['time'])

def pyrun(args, cell=None):
    print('Using <pyrun> main loop')

    # init_generator for reruns
    tmos.base.init_generator(args.system.seed)

    cell = setup_cell(args)

    fconf = 'config.txt'
    with open(fconf, 'w') as fc:
        args.write_config(fc)

    outputwhen = tfputils.Sampler.from_args(args)
    pilitrackingwhen = tfputils.Sampler.from_args(args)
    vtkwhen = tfputils.Sampler(args.system.simtime, args.system.vtkdeltat)

    # setup vtk output
    datapaths.setup_vtkdir()
    vtkbodypath = datapaths.vtkit_body()
    vtkpilipath = datapaths.vtkit_pili()
    vtkwriter.write_cell3d(cell, next(vtkbodypath))
    vtkwriter.write_pili3d(cell, next(vtkpilipath))

    #
    md3d = setup_md3d(args)
    args.Cell3d.cost_anchor_intersection = False

    # 
    kmc = pili.Kmc(1./args.Pili.k_adh)

    # setup tracking
    fout = 'trackxy.dat'
    bout = 'event.dat'
    pout = 'pilitracking.dat'
    mdout = 'mdev.dat'
    txy = wr.Track3d.setup_track_3d(fout)
    bind = wr.Track3d.setup_event_3d(bout)
    tpili = wr.PiliTrack3d(pout)
    mdev = wr.Track3d.setup_mdevent(mdout)

    tme = 0.
    after_tme = 0.
    while True:
        kmcstep =  kmc.kmcstep(cell, tme)
        after_tme, process, need_md, pidx = kmcstep

        istep = 0
        force = -1
        if need_md:
            istep, force = md3d.step(cell)
            summary = md3d.get_summary()
            tfputils.md_events(mdev, (tme, pidx, process, cell, summary))
            kmc.postmd(cell, tme)
        tme = after_tme
           
        for ev in kmc.events:
            bind.track_event(ev)

        # vtkoutput
        vtkwrite_condition = args.system.vtkwrite and vtkwhen.check(tme)
        if vtkwrite_condition:
            for target in vtkwhen.passed_targets:
                vtkwriter.write_cell3d(cell, next(vtkbodypath))
                vtkwriter.write_pili3d(cell, next(vtkpilipath))

        if args.system.track:
            state = wr.Simstate(kmcstep, cell, istep, force)
            txy.track(state)
            tpili.track_pili(tme, cell)
            if outputwhen.check(tme):
                    print(printform.format(process, tme))
        else:
            if outputwhen.check(tme):
                for output_target in outputwhen.passed_targets:
                    print(printform.format(process, output_target))
                    kmcstep = (output_target, process, need_md, pidx)
                    state = wr.Simstate(kmcstep, cell, istep, force)
                    txy.track(state)
                    tpili.track_pili(output_target, cell)

        # exit condition
        if tme > args.system.simtime:
            break


# construct the c++ Equilibrate.step which iteratively applies an integrator step 
# Use the given Integrate object (Bacterium/mdynamics.hpp)
# record force and energies
def stepper(md3d):
    def step(cell):
        print('Start Minimisation')
        md3d.reset(cell)
        ens_p = np.zeros(md3d.maxsteps)
        ens_s = np.zeros(md3d.maxsteps)
        ens = np.zeros(md3d.maxsteps)
        xyz = np.zeros((md3d.maxsteps, 3))
        axx = np.zeros((md3d.maxsteps, 3))
        fmag = np.zeros(md3d.maxsteps)
        dt = np.zeros(md3d.maxsteps)
        for istep in range(md3d.maxsteps):
            force = md3d.one_step(cell)
            ens_p[istep] = cell.energy_pili()
            ens_s[istep] = cell.energy_surface()
            print('step {} has energy'.format(istep))
            print('pili energy ', ens_p[istep])
            print('surface energy ', ens_s[istep])
            print('step {} has force'.format(istep))
            print('force', force)
            ens[istep] = ens_p[istep] + ens_s[istep]
            fmag[istep] = force
            dt[istep] = md3d.dt
            cc = cell.real_centre()
            xyz[istep,:] = np.array([cc.x,cc.y,cc.z])
            ax = cell.get_axis()
            axx[istep,:] = np.array([ax.x,ax.y,ax.z])

            if force < md3d.target:
                break
        return istep, force, (ens_p, ens_s, ens, fmag, dt, xyz, axx)
    return step


# define a minimal simulation loop that prints on each step. For debugging.
# ... "minrun" turned into a full debug run here. With vtk output and MD loop output
def minrun(args, cell):
    print('Using <minrun> main loop')
    args.system.track = True
    print(args)

    print('Initial Cell')
    print(cell)

    md3d = setup_md3d(args)
    #istep, force = md3d.step(cell)

    # debugging
    mdstep = stepper(md3d)

    fconf = 'config.txt'
    with open(fconf, 'w') as fc:
        args.write_config(fc)

    # track the numerical error in the rotation matrix
    otherror = []
    # track the overlap with the surface
    overlap = []

    # setup vtk output
    datapaths.setup_vtkdir()
    vtkbodypath = datapaths.vtkit_body()
    vtkpilipath = datapaths.vtkit_pili()
    vtkwriter.write_cell3d(cell, next(vtkbodypath))
    vtkwriter.write_pili3d(cell, next(vtkpilipath))

    # 
    kmc = pili.Kmc(1./args.Pili.k_adh)

    # setup tracking
    fout = 'trackxy.dat'
    txy = wr.Track3d.setup_track_3d(fout)
    tme = 0.
    while True:
        print('performing monte carlo step at time', tme)
        kmcstep =  kmc.kmcstep(cell, tme)
        tme, process, need_md, pidx = kmcstep
        
        ###################### DEBUGGING #################################
        print('this step pilus {} performs process {}'.format(pidx, process))
        print('rates ------------')
        for pilus in cell.pili:
            #print pilus.get_rates(), "pl, leq, diff", pilus.pl, pilus.leq, abs(pilus.pl - pilus.leq)
            print(pilus.get_rates())
        print('---------- end rates')
        ##################################################################

        istep = 0; force = -1
        if need_md:
            print('performing molecular dynamics step at time', tme)
            #istep, force = md3d.step(cell)
            istep, force, extra = mdstep(cell)
            kmc.postmd(cell, tme)
            print('number of steps {}'.format(istep))
            print('final force {}'.format(force))
            # track the numerical error in the rotation matrix
            oth = cell.body.frame.orthogonal_error()
            eun = cell.body.frame.unit_error()
            #otherror.append(oth)
            print('error in frame othogonality', oth)
            print('error in frame normality', eun)

            force_condition = force > args.fire.force_target
            step_condition = istep == args.fire.maxsteps
            if force_condition or step_condition:
                # check is called
                #cell.get_body().frame.normalise();
                #
                ens_p, ens_s, ens, fmag, dt, xyz, axx = extra
                if force_condition:
                    print("FAILED FORCE CONDITION force = {}".format(force))
                if step_condition:
                    print("FAILED TO FIND EQUILIBRIUM AFTER {} STEPS".format(istep))
                plt.plot(fmag, label='Total Force')
                plt.ylabel("force")
                print('saving to {}'.format("last_md_force.svg"))
                plt.savefig("last_md_force.svg")
                plt.show()
                plt.clf()

                plt.plot(dt, label='Timestep')
                plt.ylabel("dt")
                print('saving to {}'.format("last_md_dt.svg"))
                plt.savefig("last_md_dt.svg")
                #plt.show()
                #plt.clf()

                plot_xyz(xyz, ['x','y','z'], 'xyz')
                plot_xyz(axx, ['ax','ay','az'], 'axis')

                plt.plot(ens_p, label='Pili Energy')
                plt.plot(ens_s, label='Surface Energy')
                plt.plot(ens, label='Total Energy')
                plt.legend()
                plt.ylabel("energy")
                plt.savefig("last_md_energy.svg")
                raise MDnoconverge(
                        "MD step failed to converge. Force is {} target is {}"
                        .format(force, args.fire.force_target)
                        )

            # check for nan
            if math.isnan(cell.body.get_origin().x):
                print("MD STEP EXPLODED")
                print("CELL STATE")
                print(cell)
                print("FAILED")
                raise MDexplode(
                        "Found nan during MD step. This usually means the integration exploded"
                        )

            # vtkoutput
            if args.system.vtkwrite:
                vtkwriter.write_cell3d(cell, next(vtkbodypath))
                vtkwriter.write_pili3d(cell, next(vtkpilipath))


        print('------------------------------------------------------------')
        state = wr.Simstate(kmcstep, cell, istep, force)
        line = txy.track(state)

        # exit condition
        if tme > args.system.simtime:
            break

def plot_xyz(xyz, labels, name):
    plt.clf()
    x, y, z = xyz.T
    xl, yl, zl = labels
    plt.plot(x, label=xl)
    plt.savefig(xl+'.svg')
    plt.clf()
    plt.plot(y, label=yl)
    plt.savefig(yl+'.svg')
    plt.clf()
    plt.plot(z, label=zl)
    plt.savefig(zl+'.svg')
    plt.clf()


