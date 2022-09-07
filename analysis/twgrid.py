#!/usr/bin/env python3


"""
Need some methods to operate on tracks on the hexgrid surface.
"""

from tmos.base import Vector3d
import tmos.surface as surface
from tmos.surface import Hexc, HexSphereGrid 


import sys, os, math
import numpy as np
import matplotlib.pyplot as plt

import datapaths
import command
import parameters
from command import defaultsave
import readtrack

# for kdeplot
import fjstats


### utilities

def tov3d(arr):
    x, y, z = arr
    return Vector3d(x,y,z)

def tonpy(v):
    return np.array([v.x, v.y, v.z])

# information I Need. spherer_radius, bacteria Radius

# prints warning if defaults are returned!
def local_args():
    args = parameters.thisread()
    # 
    if args.is_default:
        print('WARNING: did not find local config.txt, default params loaded')
    R = args.cell.R
    length = args.cell.length
    sphere_radius = args.surface.sphere_radius
    return args

def local_sphere_radius():
    return local_args().surface.sphere_radius

args = local_args()
#print 'sphere_radius = {}'.format(args.surface.sphere_radius)

# reconstruct the surface (c++ object) from sphere_radius
#hsphgrid = HexSphereGrid(sphere_radius)

def visitmap():
    trs = readtrack.trackset()
    trsdst = _visitmap(trs)
    hexvisitmap(trsdst)

def _visitmap(trs):
    """
    From this set of tracks. compute the visit map. 
    Visitmap is is in relative hex coordinates
    where distance is relative to the centre of the closest hex.
    """
    hsphgrid = HexSphereGrid(local_sphere_radius())

    # trajectory
    def pertrack(tr):
        # check distance to closest circle
        head = tr.get_head()
        # throw away the first 5% because of starting config
        head = head[len(head)//20:]
        # if we use C++ data structure, need Vector3d object, # cannot vectorise
        disp = np.zeros((tr.size, 3))
        for i, pt in enumerate(head):
            x,y,z = pt
            vpt = Vector3d(x,y,0)
            vsphere = hsphgrid.get_xyz( hsphgrid.get_rq( vpt ) )
            disp[i] = tonpy(vpt - vsphere)
            #print vsphere, vpt, np.linalg.norm(disp[i])
        # compute distances
        dst = np.linalg.norm(disp, axis=1)
        return dst # 1d array of distances
    trsdst = list(map(pertrack, trs))
    return trsdst

def _hexvisitmap(ax, trsdst):
    # append the distributions for each track
    alldst = np.concatenate(trsdst)
    assert(alldst.ndim == 1)
    # plot smooth histogram
    plt.clf()
    plt.hist(alldst)
    plt.axvline(local_sphere_radius(), c='k')
    plt.xlabel('distance (\mu m)')
    plt.ylabel('weight')
    plt.tight_layout()
    plt.plot() 

@command.defaultsave(autoshow=False)
def hexvisitmap(trsdst):
    _hexvisitmap(plt.gca(), trsdst)



if __name__=='__main__':

    ff, thisargs = command.process(locals())
    ff(*thisargs)

