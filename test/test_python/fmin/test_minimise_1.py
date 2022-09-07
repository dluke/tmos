
import sys

import tmos.base
import tmos.surface as surface 
import tmos.pili as pili
import tmos.mdynamics as md
import tmos.vtkwriter as vtkwriter

from tmos.base import Vector3d

import math
import numpy as np
import matplotlib.pyplot as plt

import parameters


import tfputils
import ctfp3d

# situation 1
# no pili attached, surfaces overlapping


def senario_1():

    cfile = '/data/dan/run/sept_18/inp/testsinglepool/config.txt'
    args = parameters.thisread(cfile)
    R = args.cell.R # 0.5
    length = args.cell.length # 2

    # pele minimisation 
    trypele = ctfp3d.setup_pelefire(args)
    sf = ctfp3d.setup_cell_surface(args)

    wlcg = ctfp3d.setup_wlc_generator(args)

    # cell stuff
    cellaxis = Vector3d(1.,0.,-1.).unit()
    r2 = np.sqrt(2.5**2/2)

    delta = 0.004
    #cellmid = Vector3d(-r2,0,r2)
    cellmid = Vector3d(-r2+delta,0,r2-delta)

    body = surface.Capsule(
            cellmid,
            cellaxis,
            args.cell.R,
            args.cell.length)

    cell = pili.CellWLC3d(0, args.cell.npili, body, sf, wlcg)
    # call once
    cell.common_init()

    print 'head at ', cell.get_headpt()

    nsteps, ngrad = trypele(cell)
    print 'nsteps ', nsteps
    print 'ngrad ', ngrad

# situation 2
# vertical with pili attached, surfaces overlapping

# cell stuff
def senario_2():

    cfile = '/data/dan/run/sept_18/inp/testsinglepool/config.txt'
    args = parameters.thisread(cfile)
    R = args.cell.R # 0.5
    length = args.cell.length # 2

    # pele minimisation 
    trypele = ctfp3d.setup_pelefire(args)
    sf = ctfp3d.setup_cell_surface(args)

    wlcg = ctfp3d.setup_wlc_generator(args)

    delta = 0.004
    full_length = length + 1.
    cellmid = Vector3d(0,0,1.+full_length/2.-delta)
    cellaxis = Vector3d(0,0,-1)

    body = surface.Capsule(
            cellmid,
            cellaxis,
            args.cell.R,
            args.cell.length)

    cell = pili.CellWLC3d(0, args.cell.npili, body, sf, wlcg)
    # call once
    cell.common_init()
    print cell

    nsteps, ngrad = trypele(cell)
    print 'nsteps ', nsteps
    print 'ngrad ', ngrad


def senario_3():
    # step surface 
    cfile = '/data/dan/run/dec_18/stepped_surface/config.txt'
    args = parameters.thisread(cfile)

    hl = args.cell.length/2. + args.cell.R
    delta = 0.001
    args.cell.initcenter = (args.surface.sep/2. - hl + delta, 0, args.cell.R-delta)
    cell = ctfp3d.setup_cell(args)
    md3d = ctfp3d.setup_md3d(args)
    #print cell.get_num_contacts()
    cell.body.init_Rk()
    print cell.grad()

    md3d.step(cell)

def corner_senario():
    # step surface 
    cfile = '/data/dan/run/dec_18/stepped_surface/step_test_bottom/config.txt'
    args = parameters.thisread(cfile)
    surface = args.surface

    delta = 0.001
    # corner center
    corner = np.array([surface.sep/2. + surface.smallr,0,surface.height-surface.smallr])
    e_z = np.array([0,0,1])
    e_x = np.array([1,0,0])
    e_xz = (e_z - e_x)/np.sqrt(2.)

    args.cell.initcenter = corner + (surface.smallr + args.cell.R - delta) * e_xz
    #print 'center',args.cell.initcenter
    args.cell.initaxis = (e_x + e_z)/np.sqrt(2.)
    #print 'axis', args.cell.initaxis

    cell = ctfp3d.setup_cell(args)
    md3d = ctfp3d.setup_md3d(args)
    print cell.get_num_contacts()
    cell.body.init_Rk()
    print cell.grad()

    #md3d.step(cell)

def corner_senario_2():
    cfile = '/data/dan/run/dec_18/stepped_surface/step_test_bottom/config.txt'
    args = parameters.thisread(cfile)
    head = np.array([14.755906,    -4.682189,     5.344189])
    tail = np.array([13.310338,    -3.307736,     5.489843])
    args.cell.initcenter = (head+tail)/2.
    args.cell.initaxis = [0.722784,    -0.687226,    -0.072827]

    cell = ctfp3d.setup_cell(args)
    md3d = ctfp3d.setup_md3d(args)
    #print cell.get_num_contacts()
    #cell.body.init_Rk()
    print cell.grad()

def corner_senario_3():
    # check 4 contacts is possible
    cfile = '/data/dan/run/dec_18/stepped_surface/step_test_bottom/config.txt'
    args = parameters.thisread(cfile)
    delta = 0.001
    args.cell.initcenter = [-4.5-delta, 0, 0.5-delta]
    args.cell.initaxis = np.array([0., 1., 0.])

    cell = ctfp3d.setup_cell(args)
    
    print cell.get_num_contacts()
    sys.exit()


def main():
    #senario_3()
    #corner_senario()
    #corner_senario_2()
    corner_senario_3()

if __name__=='__main__':
    main()
