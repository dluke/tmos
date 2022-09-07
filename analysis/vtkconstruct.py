#!/usr/bin/env python3

import os, collections, math
import numpy as np

import tmos
import parameters
import datapaths

import tmos.base as base
import tmos.surface as surface
import tmos.vtkwriter as vtkwriter

import readtrack

"""
1. Reconstruct a Capsule object as well as the attached pili positions from text output.

2. Use tmos to write this as a vtk file.
"""



args = parameters.thisread()

bodytrackform ='data/bacterium_{:05d}.dat'
pilitrackform ='data/bacterium_{:05d}.dat'

basedir = 'revtk'
odirform = 'revtk/track_{:05d}/'

bodyform = os.path.join('revtk', 'track_{:05d}/revtkbody_{:06d}.vtp')
piliform = os.path.join('revtk', 'track_{:05d}/revtkpili_{:06d}.vtp')


def slice_to_range(item, stop):
    ifnone = lambda a, b: b if a is None else a
    return list(range(ifnone(item.start, 0), stop, ifnone(item.step, 1)))

def tov3d(arr):
    x, y, z = arr
    return base.Vector3d(x,y,z)

def get_track_id(trackfile):
    return int(os.path.splitext(os.path.basename(trackfile))[0].split('_')[-1])


def reconstruct_cell3d(tr):
    #tr.step_to_resolution(stepsize)
    xyz = tr.get_head()
    txyz = tr.get_trail()
    tridx = get_track_id(tr.trackfile)
    for idx, hpt, tpt in zip(slice_to_range(tr.slice, len(tr)), xyz, txyz):
        cpt = 0.5 * (hpt + tpt)
        ax =  (1./args.cell.length) * (hpt - tpt)
        body = surface.Capsule(tov3d(cpt), tov3d(ax), args.cell.R, args.cell.length)
        bfile = bodyform.format(tridx, idx)
        print('writing to ', bfile)
        vtkwriter.write_cell3d_from_body(body, bfile)

def get_pili_history(tr, pr):
    # Construct track-like object for each pilus
    pilitracks = [collections.OrderedDict() for _ in range(args.cell.npili)]
    cols = ['time', 'isbound', 'anchor', 'target'] 
    time = tr['time']
    tstep = (tr.slice.step or 1) * args.system.deltat
    def get_next_idx(offtime):
        return (tstep/args.system.deltat * np.ceil(offtime/tstep)).astype(int)
    _nbound = np.insert(pr['nbound'], -1, pr['nbound'][-1]) # extend by 1 buffer element
    binding = np.argwhere(_nbound[:-1] < _nbound[1:]).ravel()
    unbinding = np.argwhere(_nbound[1:] < _nbound[:-1]).ravel()
    bindidx = get_next_idx(pr['time'][binding])
    unbindidx = get_next_idx(pr['time'][unbinding])
    assert(binding.size - unbindidx.size + pr['nbound'][-1])

    ...


def reconstruct_pili(tr, pr):
    phistory = get_pili_history(tr,pr)
    #


##################################################################################
# TESTING

exdir = "/data/dan/run/dec_18/stepped_surface/low_high_pair/high_step_3x4x5/eps_00020.0000_pilivar_02.00000_height_001.000/data/"
mainf = "bacterium_00000.dat"
pilif = "pili_00000.dat"
def test_reconstruct_cell3d():
    tr = readtrack.Track(os.path.join(exdir, mainf))
    reconstruct_cell3d(tr)

def test_reconstruct_pili():
    tr = readtrack.Track(os.path.join(exdir, mainf))
    tr.slice = slice(0, None, 10)
    pr = readtrack.Track(os.path.join(exdir, pilif))
    reconstruct_pili(tr, pr)


if __name__=='__main__':


    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('-n', default=[0], action='store', nargs='+', type=int,
                        help='Use which tracks to build vtk')
    parser.add_argument('-s', '--step', default=1, action='store', type=int,
                        help='slice.step')
    parser.add_argument('--body-only', action='store_true',
                        help='Only compute the body vtk file')

    aparse = parser.parse_args()

    for n in aparse.n:
        odir = odirform.format(n)
        print("making directory ", odir)
        datapaths.force_mkdir(odir)
        bform = bodytrackform.format(n)
        print("reading data from {}".format(bform))
        tr = readtrack.Track(bform)
        tr.slice = slice(0, None, aparse.step)
        reconstruct_cell3d(tr)
        if not aparse.body_only:
            pass
            #reconstruct_pili3d ...  

