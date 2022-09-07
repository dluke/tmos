
import sys

import numpy as np


import tmos
from tmos.base import Vector3d
import tmos.surface as surface 
import tmos.vtkwriter as vtkwriter


import readtrack 


"""
# directory structure
simdir/data -- data directory 
simdir/vtk  
simdir/vtk/id_0001/vtkbody_{:06d}.vtp
simdir/vtk/id_0001/vtkpili_{:06d}.vtp

simdir/vtk/clips/id_000x_frames_{:06d}_{:06d}/

"""
# define forms 
clipdir = "vtk/clips/"
idform = "id_{:04d}"
clipform = "id_{:04d}_frames_{:06d}_{:06d}"

exdir = "/data/dan/run/dec_18/plane_set_parameters/crawling_4x4x5/tau_02.000_pilivar_07.00000_free_tau_eovers_0001.000/data"

# Get rotation matrices M_i, apply rotation and translation to attached pili anchors and then use pili length to get the bound pili position

def txt_to_vtk(idx, frslice=slice(0,None), exdir=''):
    # load the track
    tr = readtrack.Track(readtrack.find_track_files(ddir=exdir)[idx])
    pr = readtrack.Track(readtrack.find_track_files(form='pili_*.dat', ddir=exdir)[idx])
    tr.slice = frslice
    R = 0.5
    length = 2.

    # shape = (T, 3)
    cxy, axis, e1 = tr.get_frame()
    translate = np.diff(cxy, axis=0)
    translate = np.insert(translate, 0, np.zeros(3), axis=0)

    # construct capsules
    for i, xy in enumerate(cxy):
        #  construct c++ vectors
        vxy = Vector3d(*xy)
        vaxis = Vector3d(*axis[i])
        ve1 = Vector3d(*e1[i])
        caps = surface.Capsule(vxy, vaxis, ve1, R, length)


txt_to_vtk(0, exdir=exdir)

# ..

NOT FINISHED
