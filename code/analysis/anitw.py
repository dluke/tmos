#!/usr/bin/env python3

import command
import twanimation
import readtrack
from readtrack import trackset
import matplotlib.pyplot as plt

#################################################################################
## animations

def shapeidx(idx=0, sample=10, savefile=''):
    # sample = 1
    twanimation.outline(plt.gcf(), [trackset()[idx]], sample, camera='follow', savefile=savefile)

def shape(sample=10, to_save=False):
    #sample = 1 # 10 seconds between frames
    twanimation.outline(plt.gcf(), trackset(), sample, savefile=to_save)


def lt():
    twanimation.ani_longtracks(trackset())



if __name__=='__main__':

    ff, thisargs = command.process(locals())
    ff(*thisargs)
