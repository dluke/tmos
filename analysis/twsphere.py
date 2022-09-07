
# spherical distance calculation


import numpy as np
import parameters
import readtrack

import matplotlib.pyplot as plt

args = parameters.thisread()

def get_params():
    R = args.cell.R
    sphere_radius = args.surface.sphere_radius 
    params = {'R':R, 'sphere_radius':sphere_radius}
    params['shape'] = args.sphere.shape

    if args.surface.shape == 'sphere':
        hshell_r = sphere_radius + R 
    elif args.surface.shape == 'shell':
        hshell_r = sphere_radius - R 
    else:
        raise RuntimeError("no spherical shape found.")

    return params

e_z = np.array([0,0,1])

#https://en.wikipedia.org/wiki/Great-circle_distance
def _disp(track):
    time = track['time']
    refn = e_z
    try:
        hptn = track.get_head()
    except:
        print(track['x'].size)
        print(track['y'].size)
        print(track['z'].size)
        raise
    # normalise
    hptn = hptn/np.sqrt(np.sum(hptn**2, axis=1))[:,np.newaxis]
    # distances
    tile_ez = np.tile(e_z, track.size).reshape((-1,3))
    disps  = hshell_r * np.arccos(np.einsum('ij,ij->i', hptn, tile_ez))
    return disps


def msdlike(trs):
    enddisp = []
    for tr in trs:
        enddisp.append( _disp(tr)[-1]**2 )
    x = np.mean(enddisp)
    return x

def main():
    trs = readtrack.trackset()
    msdlike(trs)

if __name__=='__main__':

    main()


