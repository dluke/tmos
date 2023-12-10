#!/usr/bin/env python3

import os
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
from collections import OrderedDict

import command
import astrack
import _fj
import _analysefj
import fjstats

sample = list(range(1))
def trd():
    #ft = _fj.npyloadall()
    timestep = 0.1
    ft = _fj.npyloadall(sample)
    disps = _analysefj.track_disps(ft, timestep)
    for disp in disps:
        plt.clf()
        #plt.hist(disp, bins=20)
        plotutils.kdeplot(disp)
        plt.show()
        

if __name__=='__main__':
    ff, thisargs = command.process(locals())
    ff(*thisargs)

