#!/usr/bin/env python

## Spitting out the state of the pele minimisation into a temporary file 
## Read this to understand how the minimisation failed

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  

from command import defaultsave

statefile = 'pele_min_state.out'
pfile = 'pgrad.out'
sfile = 'sgrad.out'

ll = ['x', 'y', 'z']

def read_statefile(statefile=statefile):
    data = np.loadtxt(statefile, dtype=float)
    state = data[:,:6]
    energy = data[:,6]
    rms = data[:,7]
    return state, energy, rms

@defaultsave()
def plot_each_xyz(state):
    xyz = state.T[:3]
    fig, axes = plt.subplots(1,3, figsize=(15,10))
    for zi, ax in enumerate(axes):
        ax.plot(xyz[zi])
        ax.set_xlabel(ll[zi])

@defaultsave()
def plot_xyz(state):
    x, y, z = state.T[:3]
    plt.plot(x, label='x')
    plt.plot(y, label='y')
    plt.plot(z, label='z')
    plt.legend()

@defaultsave()
def plot_pxyz(state):
    x, y, z = state.T[3:6]
    plt.plot(x, label='px')
    plt.plot(y, label='pz')
    plt.plot(z, label='py')
    plt.legend()

def pmag(state):
    p = state[:,3:6]
    return np.linalg.norm(p, axis=1)

@defaultsave()
def plot_pmag(state):
    plt.plot(pmag(state))

@defaultsave()
def plot_energy(energy):
    plt.plot(energy)
    ymax = min(1000, energy.max())
    plt.ylim([0,ymax])

@defaultsave()
def plot_rms(rms):
    plt.plot(rms)
    ymax = min(1000, energy.max())
    plt.ylim([0,ymax])

@defaultsave()
def plot_spz():
    pdata = np.loadtxt(pfile, dtype=float)
    sdata = np.loadtxt(sfile, dtype=float)
    fig, axes = plt.subplots(1,3, figsize=(15,10))
    char = ['x','y','z']
    for zi, ax in enumerate(axes):
        ax.plot(pdata[:,zi], label='pili {} gradient'.format(char[zi]))
        ax.plot(sdata[:,zi], label='surface {} gradient'.format(char[zi]))
        ax.legend()

def xyz(state):
    #plot_energy(energy)
    plot_xyz(state)
    plot_pxyz(state)

def xyzpath(state):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    x, y, z = state[::10].T[:3]
    ax.scatter(x,y,z,c=np.arange(x.size))

if __name__=='__main__':

    state, energy, rms = read_statefile() 

    plt.clf()
    plot_xyz(state)

    plt.clf()
    plot_each_xyz(state)

    plt.clf()
    #xyzpath(state)
    plt.clf()
    plot_spz()
    plt.clf()
    plot_energy(energy)
    plt.clf()
    plot_rms(rms)

    plt.clf()
    plot_pmag(state)
