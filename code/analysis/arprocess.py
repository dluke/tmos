
import os
import _fj
import shapeplot
import filesystem
import numpy as np
import readtrack

import matplotlib.pyplot as plt
"""

Prep for analysing data using autoregression model according to
https://www.nature.com/articles/ncomms8516#Sec19
"""

# setup data directory for this analysis module
datadir = 'velocitydata/'
filesystem.safemkdir(datadir)
def _namify(idx):
    basename  ='head_velocities_{:04d}.txt'.format(idx)
    return os.path.join(datadir, basename)

### load data
# input track to process here
trackidx = 1104
fjdata = _fj.npyload(trackidx)
# print the data columns
print('columns', fjdata.colmap)

# save the velocity data for AR superstatistical analysis
trackidx = [0]
vel = fjdata.get_velocity(trackidx)
nidx = [fjdata.nidx[i] for i in trackidx]
name = _namify(nidx[0])
print('saving to ', name)
np.savetxt(name, vel)


# vis
def quick_show(ft):
    trs = ft.get_tracklike()
    ax = plt.gca()
    shapeplot.longtracks(ax, trs)
    plt.show()
