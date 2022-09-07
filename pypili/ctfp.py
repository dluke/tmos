
# dep
# This was for 2d simulation, throw it out


# The tfp.py program defines the main loop of the simulation 
# The runtfp program just defines the parameters so it is an input file
# The runtfp program can also choose the python or C++ backend

# debugging
#from test_utils import do_cprofile
import cProfile

# 
import sys, os, time
import logging


# used for initalising Cells
import tmos.base 
#import tmos.surface as surface
import tmos.pili as pili
import tmos.mdynamics as md
import tmos.vtkwriter as vtkwriter

# using
from tmos.base import Vector3d

# used for data I/O and analysis
import numpy as np
import wrinterval as wr
from math import *
from matplotlib import pyplot as plt

import datapaths as dpath

# duplicated from utils
plist = ['elongation','retraction','switching','rotation','attach','release']
proc = dict(list(zip(plist, list(range(len(plist))))))

from collections import OrderedDict

tfputils import *

# setup a basic simulation but with the state of some pili pre-set
def detailed(args):

    initpt = Vector3d(0.,0.,0.)
    initaxis = Vector3d(1.,0.,0.)
    cell = pili.Cell(0, initpt, initaxis)

    tme = 0.

    p0 = cell.pili[0]
    p1 = cell.pili[1]
    #p2 = cell.pili[2]
    unit = 0.05 # \micro m
    p0.pl = 5.; p0.leq = 5. - 2 * unit
    p0.axis = Vector3d(0.,1.,0.)
    p0.axisEq = Vector3d(0.,1.,0.)
    p1.axis = Vector3d(0.,-1.,0.)
    p1.axisEq = Vector3d(0.,-1.,0.)
    p1.pl = 3.; p1.leq = 3. + 5 * unit

    cell.attach(p0, tme)
    cell.attach(p1, tme)
    #cell.attach(p2, tme)

    mdmin = setup_MD(args)

    # vtk data output setup
    dpath.setup_vtkdir()

             
    vtkwriter.write_cell(cell, dpath.next_vtk(), 0, 1, 0.)


    # setup tracking
    fout = 'trackxy.dat'
    txy = wr.Trackxy.setup_track(fout)
    while True:
        tup =  pili.kmcstep(cell, tme)
        tme, process, need_md, pidx = tup
        print('process ', process, ' pidx ', pidx)

        istep = 0
        force = -1 
        if need_md:
            istep, force = mdmin.step(cell)
        path = dpath.next_vtk()
        procid = proc[process]
        vtkwriter.write_cell(cell, path, procid, pidx, tme)

        # args.track
        line = txy.track(tme, process, cell, istep, force)

        # exit condition
        if tme > args.simtime:
            break


def md_senario(args):

    # setup a directory name of storing the vtk output for a suspicious event
    mdvtkdirform = 'suspicious_event_t_{:0.4f}'


    # standard intialisation
    initpt = Vector3d(0.,0.,0.)
    initaxis = Vector3d(1.,0.,0.)
    cell = pili.Cell(0, initpt, initaxis)

    mdmin = setup_MD(args)
    dpath.setup_vtkdir()

    fout = 'trackxy.dat'
    txy = wr.Trackxy.setup_track(fout)

    # standard vtk outputs 
    path = dpath.next_vtk()
    maxsteps = 1000

    # initialise arrays for md step data
    tetas = np.zeros(maxsteps)
    pls = np.zeros(maxsteps)
    ens = np.zeros(maxsteps)

    # output a vtk fileset representing a single MD equilibration
    fullmdvtk = False
    # for this condition we would like to keep track of the cell axis
    lastaxis = cell.axis
    thisaxis = cell.axis

    tme = 0.
    while True:
        tme, process, need_md, pidx = pili.kmcstep(cell, tme)

        istep = 0
        force = -1 
        if need_md:

            # copy the cell so re can try forking the simulation state
            holdcell = cell.copy()

            ###################################################
            mdmin.reset()  
            ens = np.zeros(maxsteps)
            fns = np.zeros(maxsteps)
            x = np.zeros(maxsteps)
            y = np.zeros(maxsteps)
            for istep in range(args.maxsteps):
                force = mdmin.one_step(cell)
                fns[istep] = force
                ens[istep] = cell.energy()
                head = cell.hexpt
                x[istep] = head.x
                y[istep] = head.y
                # output vtk file

                if force < mdmin.target:
                    break

            ###################################################
            print(tme, process, istep)


            # Find out why the MD step didn't converge
            if (istep > 998):
                print('tme ', tme)
                path = dpath.next_vtk()
                procid = proc[process]
                vtkwriter.write_cell(cell, path, procid, pidx, tme)

                #plt.plot(ens[:i])
                plt.plot(fns)
                plt.show()
                plt.scatter(x, y, c=list(range(maxsteps)))
                plt.gca().set_aspect('equal')
                plt.locator_params(nticks=3)
                plt.show()
                
        # args.track
        line = txy.track(tme, process, cell, istep, force)

        # exit condition
        if tme > args.simtime:
            break

        istep += 1

#dep
def ka_senario(args):

    initpt = Vector3d(0.,0.,0.)
    initaxis = Vector3d(1.,0.,0.)
    cell = pili.Cell(0, initpt, initaxis)

    mdmin = setup_MD(args)
    dpath.setup_vtkdir()

    fout = 'trackxy.dat'
    txy = wr.Trackxy.setup_track(fout)

    pilus = cell.pili[0]
    #pilus.axis = pilus.axis.xyrotate(pi/8.)
    pilus.pl += -1.00
    pilus.isbound = 1

    path = dpath.next_vtk()
    procid = 0; pidx = 0
    tme= 0.
    vtkwriter.write_cell(cell, path, procid, pidx, tme)
    print('before')
    print('isbound', cell.pili[0].isbound)
    print('axisEq ', cell.pili[0].axisEq)
    print('axis ',cell.pili[0].axis)
    #istep, force = mdmin.step(cell)
    maxsteps = 1000

    tetas = np.zeros(maxsteps)
    pls = np.zeros(maxsteps)
    #target = 1 * 10**-5
    target = pili.Pili.ks / 10**3.
    try: 
        for i in range(maxsteps):
            force = mdmin.one_step(cell)
            print(cell.hexpt)
            #tetas[i] = cell.pili[0].teta() - cell.pili[0].tetaEq()
            pls[i] =  cell.pili[0].pl - cell.pili[0].leq
            if force < target:
                break
    finally:
        #plt.plot(tetas[:i])
        plt.plot(pls[:i])
        plt.show()

    print('after')
    print('isbound', cell.pili[0].isbound)
    print('axisEq ', cell.pili[0].axisEq)
    print('axis ',cell.pili[0].axis)
    path = dpath.next_vtk()
    vtkwriter.write_cell(cell, path, procid, pidx, tme)

def testvtk(args, vtkall=False):
    print('Using <testvtk> mainloop')


    initpt = Vector3d(0.,0.,0.)
    initaxis = Vector3d(1.,0.,0.)
    cell = pili.Cell(0, initpt, initaxis)


    mdmin = setup_MD(args)

    # vtk data output setup
    dpath.setup_vtkdir()


    simtime = args.simtime
    deltat = 0.1
    # setup when to ouput
    outt = np.arange(0., simtime+(1.1*deltat), deltat)

    outi = 0
    outtarget = outt[outi]
    
    fout = 'trackxy.dat'
    txy = wr.Trackxy.setup_track(fout)
    tme = 0.
    while True:
        tme, process, need_md, pidx = pili.kmcstep(cell, tme)

        istep = 0
        force = -1 
        if need_md:
            istep, force = mdmin.step(cell)
            #get_pili_dtheta(cell)

            procid = proc[process]
            if vtkall:
                #print 'step {}, write cell'.format(istep)
                print('tme ', tme)
                path = dpath.next_vtk()
                vtkwriter.write_cell(cell, path, procid, pidx, tme)
            else:
                if tme >= outtarget:
                    print('tme ', tme)
                    path = dpath.next_vtk()
                    vtkwriter.write_cell(cell, path, procid, pidx, tme)

                    outi += 1
                    outtarget = outt[outi]


        # args.track
        line = txy.track(tme, process, cell, istep, force)

        # exit condition
        if tme > args.simtime:
            break

        istep += 1



# define a minimal simulation loop that prints on each step 
# for debugging
allsteps = False
def minrun(args):
    print('Using <minrun> main loop')

    initpt = Vector3d(0.,0.,0.)
    initaxis = Vector3d(1.,0.,0.)
    cell = pili.Cell(0, initpt, initaxis)

    mdmin = setup_MD(args)

    # setup tracking
    fout = 'trackxy.dat'
    txy = wr.Trackxy.setup_track(fout)
    tme = 0.
    while True:
        tup =  pili.kmcstep(cell, tme)
        tme, process, need_md, pidx = tup

        pilus = cell.pili[pidx]
        if allsteps:
            print('time, tme')
            print('rates')
            print(pilus.get_rates())
            print('pidx', pilus.idx)
            print() 

        istep = 0
        force = -1 
        #print 'step'
        if need_md:
            istep, force = mdmin.step(cell)

            print('tme, process, need_md, pidx')
            print(tme, process, need_md, pidx)
            perf = 'performing equilibration, istep = {}, force  = {}'
            print(perf.format(istep, force))

            print('all rates')
            for pilus in cell.pili:
                print(pilus.get_rates())
            print('pili average length is {}'.format(get_pili_avg_len(cell)))
            print('pili all lengths')
            print(get_pili_lens(cell))
            print() 


        # args.track
        line = txy.track(tme, process, cell, istep, force)

        # exit condition
        if tme > args.simtime:
            break


# main simulation loop 
def run(args, idx, cores=1):
    # time the program execution
    start_t = time.time()

    # initialise the random number generator
    import random
    seed = random.randrange(sys.maxsize)
    #print 'using random seed ', seed
    pili.init_generator(seed)

    # setup logging before we do anything else
    logging.basicConfig(filename='out.log', level=logging.DEBUG)
    # disable logging
    #logging.disable(logging.CRITICAL)

    # output
    fconf = 'config.txt'

    # initialize
    # create some cells, with pili and then start the update cycle
    initpt = Vector3d(0.,0.,0.)
    initaxis = Vector3d(1.,0.,0.)
    cell = pili.Cell(idx, initpt, initaxis)

    # output the parameters before the run starts
    params = compile_params(args)
    write_params(params, fconf)

    # start the update cycle
    # MD variables
    mdmin = setup_MD(args)

    # setup data directory and output file
    datdir = wr.datdir
    if not os.path.exists(datdir):
        os.mkdir(datdir)
    # track output
    name_out = os.path.join(datdir, 'bacterium_{:05d}.dat')
    # pili binding/unbinding
    fpili = os.path.join(datdir, 'pili_{:05d}.dat')

    # vtk output directory
    dpath.setup_vtkdir()

    if args.profile:
        profile = cProfile.Profile()
        profile.enable()

    onecellargs = (cell, args, mdmin, fpili, name_out)
    p = Process(target=simulate_cell, args=onecellargs)
    return p
    #p.start()


    if args.profile:
        profile.disable()
        profile.print_stats()

    end_t = time.time()
    exec_t = end_t - start_t
    print('Exited normally after {} seconds'.format(exec_t))



