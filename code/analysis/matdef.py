
import numpy as np

#################################################################################
# FanJin data format information

# trackkeys
TIMEID = 'time'
XLEAD = 'x'
YLEAD = 'y'
XTRAIL = 'trail_x'
YTRAIL = 'trail_y'
CENTROIDX = 'center_x'
CENTROIDY = 'center_y'
ORIENTATION = 'orientation'
LENGTH = 'length'
WIDTH = 'width'
ECCENTRICITY = 'eccentricity'

# compare this with simulated track
#     track_3d = [
#         'time', 'process', 'x', 'y', 'z', 'trail_x', 'trail_y', 'trail_z',
#         'ax_x', 'ax_y', 'ax_z', 'rax_x', 'rax_y', 'rax_z',
#         'npili', 'nbound', 'ntaut', 'ncontacts', 'pl_avg', 
#         'nsteps', '|force|'
#         ]

# common columns are time, x, y, trail_x, trail_y 
# orientation is bit weird in FJ data so we recompute it from the position data
# In addition we can consider length and width of experimental bacteria which are variables

def make_dtype(cols):
    dtype = { 'names': cols, 'formats': ['f8' for _ in cols]}
    return np.dtype(dtype)

#
cn = [TIMEID, XLEAD, YLEAD, XTRAIL, YTRAIL, CENTROIDX, CENTROIDY, ORIENTATION, LENGTH, WIDTH, ECCENTRICITY]
fjcolmap = dict(list(zip(cn, list(range(len(cn))))))

# trackdatakeys
ORIGINAL = 0
DENOISED = 1 
VELOCITY = 2

# property of experimental data which we know
TIMESTEP = 0.1

PLEAD = [TIMEID, XLEAD, YLEAD]
PTRAIL = [TIMEID, XTRAIL, YTRAIL]
PLEADTRAIL = [TIMEID, XLEAD, YLEAD, XTRAIL, YTRAIL]
PCENTROID = [TIMEID, CENTROIDX, CENTROIDY]

CENTROID = [CENTROIDX, CENTROIDY]
# note that leading and trailing poles must be determined by velocity
# but that means if the bacteria switch directions then these would need to be updated
LEAD= [XLEAD, YLEAD]
TRAIL = [XTRAIL, YTRAIL]
AXES = ['ax_x', 'ax_y']


#################################################################################
# Yow Ren data file format information

yrfile = '/home/dan/twitching/gdrive_sync/Yow-Ren experimental data/drive-download-20180411T032311Z-001/tr_0um_180601.txt'
yrtimestep = 1.

yrpix_to_mu = 0.36


yrnames = ['x', 'y', 'length', 'width', 'orientation', 'frame number', 'identification number']
constnames = ['center_x', 'center_y', 'length', 'width', 'orientation', 'time', 'cellid']
yrcolmap = dict(list(zip(constnames, list(range(len(constnames))))))



### For high time resolution data 
yrdistance_names = ['center_x', 'center_y', 'length', 'width', 'trail_x', 'trail_y', 'x', 'y']
#
"""
tr file is 11 column matrix with columns that read:

x center, y center, length, width, orientation, x lead, y lead, x tail, y tail, frame, ID

positions and lengths are in units of pixels. the conversion is 0.042 um per pixel
frame conversion is 1 frame per 0.1 s
orientation is angle in degrees with respect to x axis. 
"""
constnames2 = ['center_x', 'center_y', 'length', 'width', 'orientation', 'x', 'y', 'trail_x', 'trail_y',
        'time', 'identification number']
yrcmolmap2 = dict(list(zip(constnames2, list(range(len(constnames2))))))
