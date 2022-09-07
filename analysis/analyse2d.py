#!/usr/bin/env python3

import sys, os
import numpy as np
import matplotlib.pyplot as plt
import command
from command import defaultsave

import readtrack
from readtrack import Track


"""
dep -- but update these routines and move them to twanalyse.py or astrack.py
"""


@defaultsave(True)
def x(tmeslice=None):
    # get x, y coordinates
    tr= Track()
    track = tr.slice(tmeslice)

    # just plot the track
    x = track['x']
    time = track['time']

    for time, x in zip(time, x):
        if x > 20.:
            print(time)
    plt.plot(time, x)


@defaultsave(True)
def xandxdot(tmeslice=None):
    # get x, y coordinates
    tr= Track()
    track = tr.slice(tmeslice)

    # just plot the track
    x = track['x']
    y = track['y']
    time = track['time']

    xdot = x[1:] - x[:-1]
    tme = time[:-1]
    plt.plot(tme, xdot)



# """velocities in 2d dimensions
# """
# @defaultsave(True)
# def debug_vel():
#     """
#     plot velocity except draw vertical lines on release events 
#     """
#     tr = Track()
#     track = tr.track

#     time = track['time']
#     x = track['x']
#     y = track['y']
#     pl = track['pl_avg']
#     plt.plot(time, pl)
#     plt.savefig('plots/average_pilus_length.png')
#     plt.clf()

#     def plot_release():
#         proc = track['process']
#         tid = np.argwhere(proc == 'release')
#         for tt in time[tid-1]:
#             plt.axvline(x=tt, c='r')
#     plot_release()

#     velocity= np.sqrt( (x[1:]-x[:-1])**2 + (y[1:]-y[:-1])**2 )
#     tme = time[:-1]
#     plt.plot(tme, velocity)


@defaultsave(True)
def vel():
    tr = Track()
    time = tr['time']
    x = tr['x']
    y = tr['y']

    disp = np.sqrt( (x[1:]-x[:-1])**2 + (y[1:]-y[:-1])**2 )
    deltat = time[1:] - time[:-1]
    vel = np.divide(disp, deltat)
    plt.plot(time[:-1], vel)


if __name__=='__main__':
    ff, thisargs = command.process(locals())
    ff(*thisargs)
