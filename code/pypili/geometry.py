# helpful vector geometry functions

import numpy as np 

pi = np.pi
exp = np.exp
sin, cos = np.sin, np.cos
ex, ey = np.identity(2)


def theta_axis(theta):
    return np.array([cos(theta), sin(theta)])

# the signed angle in [-pi, pi]
def dotangle(u, v):
    dot = np.dot(u,v)
    # floating point errors but this is an ugly fix
    if dot > 1.: dot = 1. 
    if dot < -1.: dot = -1.
    return np.sign(np.cross(u,v)) * np.arccos(dot)

def norm(arr):
    #return np.sqrt(np.sum(np.square(arr)))
    return np.sqrt(np.sum(arr**2.))

# 2dimensions specific
def perp(v):
    x, y = v
    return np.array([-y, x])

# rotation in two dimensions 
def rotate2d(axis, angle):
    cosd = np.cos(angle)
    sind = np.sin(angle)
    x, y = axis
    return np.array([cosd*x - sind*y, sind*x + cosd*y])


