#!/usr/bin/env python3

import readtrack
import plotutils

import matplotlib.pyplot as plt
import numpy as np

import command
from command import defaultsave



@defaultsave()
def plane():
    att = readtrack.Track('attachment.dat')
    _plane(plt.gca(), att)

def _plane(ax, att, geom=(0.5,2.0)):

    x = att.track['attach_x']
    y = att.track['attach_y']
    z = att.track['attach_z']

    with plt.style.context(['fivethirtyeight', plotutils.get_style('ft8')]):
        kw = {'alpha':0.1}
        ax.scatter(x,y,**kw)

        # annotate the bacteria position
        pill = plotutils.pillgeometry(*geom)
        pill_x, pill_y = pill.T
        ax.plot(pill_x, pill_y, color='black')

        ax.set_aspect('equal')

        # limits ...
        ax.set_xlim(-1,7)
        ax.set_ylim(-5,5)
        ax.set_xticks([0,4])
        ax.set_yticks([-4,0,4])


if __name__=='__main__':

    ff, thisargs = command.process(locals())
    ff(*thisargs)



