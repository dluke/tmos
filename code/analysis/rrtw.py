#!/usr/bin/env python3

import sys, os
import command

import matplotlib.pyplot as plt

import numpy as np
import rtw
import twstep


"""
Crossing analysis for multiple sets of simulation data
"""


# analyse the low and high step pair of identially shaped data

pairdirs = ['/data/dan/run/dec_18/stepped_surface/low_high_pair/high_step_3x4x5',
            '/data/dan/run/dec_18/stepped_surface/low_high_pair/low_step_3x4x5']


def get_pair():
    thisdir = os.getcwd()
    cubes = []
    for pd in pairdirs:
        os.chdir(pd)
        cube = rtw.get_Track_cube()
        cubes.append(cube)
    os.chdir(thisdir)
    high_cube, low_cube = cubes

    return high_cube, low_cube

def build_crossing_data(dc):
    cubes = get_pair()
    scube = dc.sumcubes(*cubes)
    crosscube = dc.autoloading(twstep._utouchstate, scube, dtype=object)
    return crosscube


@command.defaultsave()
def _project_counts(dc, counts, upcounts, downcounts):
    """for each dimension, reduce crossing counts to that dimension and plot"""
    fig, axes= plt.subplots(1,3, sharey=True)
    ax0 = axes[0]
    ax0.set_ylabel( 'No. Crossings/track' )
    fig.set_size_inches(5*3, 5.5)
    ntracks = 20.
    c, uc, dwc = [dc.apply(np.mean, cu) for cu in [counts, upcounts, downcounts]]

    for i, ax in enumerate(axes):
        ax.set_xlabel( dc.pnames[i] )
        tosum = list(range(3))
        tosum.remove(i)
        cc = np.mean(c, axis=tuple(tosum))
        upcc = np.mean(uc, axis=tuple(tosum))
        downcc = np.mean(dwc, axis=tuple(tosum))
        # standard error
        # ...

        ax.plot(dc.basis[i], cc, marker='D', lw=0, label='total')
        ax.plot(dc.basis[i], upcc, marker='D', lw=0, label='up step')
        ax.plot(dc.basis[i], downcc, marker='D', lw=0, label='down step')

    axl = axes[-1]
    axl.legend()
    plt.suptitle(
            '20 tracks, {:d} seconds per track (~12 cores, 10 hours)'.format(2000))
    plt.tight_layout()

def crossing_counts():
    dc = rtw.DataCube(pairdirs[0])
    crosscube = build_crossing_data(dc)
    #counts = dc.autoloading(twstep._crosscount, crosscube)
    counts = dc.apply(twstep._crosscount, crosscube)
    upcounts = dc.apply(twstep._upcrosscount, crosscube)
    downcounts = dc.apply(twstep._downcrosscount, crosscube)
    _project_counts(dc, counts, upcounts, downcounts)


if __name__=='__main__':
    ff, thisargs = command.process(locals())
    ff(*thisargs)









