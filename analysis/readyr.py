
"""Read yowren tracks and use ExpTracks object to store since this is the master object

The datatypes controlled by ExpTracks object will then be
1. FanJin data
2. Simulated Data
3. YowRen data
"""
import os
import collections
import numpy as np
import itertools
import command

import matdef
import readmat
import shapeplot
import fjstats

from matplotlib import pyplot as plt

timestep = 30.
yrdir = '/home/dan/twitching/yr/original'
yrnames = ['0um_180812.txt',   '0um_180828.txt',   '0um_180901.txt']
yrfiles = [os.path.join(yrdir, x) for x in yrnames]


# 
partial = '/home/dan/twitching/pypili/src/analysis/yrdata'
plots = '/home/dan/twitching/pypili/src/analysis/yrplots'

#############################################
# TESTING 
@command.defaultsave()
def plot_yrcaps(yr):
    ax = plt.gca()
    ax.set_aspect('equal')
    ax = shapeplot.exp_capsdraw(ax, yr)

def init():
    yrtracks = readmat.ExpTracks.useyr(yrfiles[0], timestep=30.)
    return yrtracks

# load the high resolution YR data
def inithres():
    highresfile = '/home/dan/twitching/yr/190314/tr.txt'
    pix_to_mu = 0.042
    yrtracks = readmat.ExpTracks.useyr(highresfile, timestep=0.1,
            colmap=matdef.yrcolmap2, yrpix_to_mu=pix_to_mu)
    yrtracks.override_orientation_using_lead_trail()
    return yrtracks

def tracklike():
    yt = readmat.ExpTracks.useyr(yrfiles[0], timestep=30.)
    yrtrs = yt.get_tracklike()
    # check the number of frames in tracks, 
    for tr in yrtrs:
        print(len(tr))
    return yrtrs

def check_orient(yr):
    ors = yr.get_columns(['orientation'])
    #for otrack in ors:
        #plt.plot(otrack)
        #plt.show()
    
    def diffor(otrack):
        dor = otrack[1:] - otrack[:-1]
        return dor

    ddor = [diffor(otrack) for otrack in ors]
    for i, dor in enumerate(ddor):
        print(i, np.max(dor))
    wor = yr.get_whole_col('orientation')
    print(np.max(wor))
    print(np.min(wor))

#############################################

def yrgetstats():
    print("Loading YR track data")
    yt = readmat.ExpTracks.useyr(yrfiles[0], timestep=30.)

    # aspect ratio
    print("Aspect Ratio") 
    whl = yt.get_columns(['width', 'length'])
    whs = np.array([np.mean(length/width) for width, length in whl])
    # 

    yrtrs = yt.get_tracklike()

    stats = collections.OrderedDict([('aspect_ratio',whs)])
    stats.update(fjstats.getstats(yrtrs))
    return stats


def yranalyse():
    stats = yrgetstats()
    for statitem in list(stats.items()):
        fjstats.process(statitem, plots)

# pulled from the output of yrgetstats
meankmsd = 1.3815582973230018


if __name__=='__main__':


    #tracklike()
    #yranalyse()

    yt = inithres()

