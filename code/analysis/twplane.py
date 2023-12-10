


import sys, os, math
import numpy as np
import matplotlib.pyplot as plt

from tmos.surface import SegPlane


# copy of twstep.draw_step_line

def draw_period_line(ax, xlims):
    args = parameters.thisread()
    template_surface = surface.SegPlane(args.surface.planeangle,args.surface.panewidth)
    # sep repeating period (x direction)
    sep = template_surface.xperiod
    xn, xm = xlims
    xdist = xm - xn
    nsteps = math.ceil(xdist/sep) + 1
    # draw this many steps in both directions
    xgrid = np.arange(-nsteps*sep -5., nsteps*sep+5.+sep/2., sep)
    lines = []
    for i, x in enumerate(xgrid):
        line = ax.axvline(x, alpha=0.4, color='k', linestyle='--')
        lines.append(line)
    return lines

