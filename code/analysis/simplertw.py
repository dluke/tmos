#!/usr/bin/env python3

"""
move one time plotting routines for rtw here 
"""

import sys, os, math
import numpy as np
import matplotlib.pyplot as plt

import command
import rtw
from rtw import DataCube

def plot_kmsd():
    dc = rtw.DataCube()
    metacube = dc.load_meta()
    base = dc.basis[0]
    kmsd = [meta['kmsd'] for meta in metacube.ravel()]
    kmsd_std = [meta['kmsd_std'] for meta in metacube.ravel()]
    print(kmsd_std)
    style = {'linestyle':'None', 'marker':'o'}
    plt.errorbar(base, kmsd, yerr=1.96*kmsd_std, **style)
    plt.show()

plot_kmsd()

