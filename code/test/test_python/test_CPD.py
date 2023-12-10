
import glob

import matplotlib.pyplot as plt
import numpy as np
sin = np.sin
cos = np.cos
norm = np.linalg.norm
def unit(v):
    return v/norm(v)

WLCconfigurations_path = '../single_pilus/CPDfunctions/Cdata_*.npy'
def read_configuration_npy(path=WLCconfigurations_path):
    Cfiles = sorted( glob.glob(WLCconfigurations_path) )
    n = 15
    return [np.load(cfile).astype('float64') for cfile in [Cfiles[n]]]

def skew_matrix(v):
    skv = np.roll(np.roll(np.diag(v.flatten()), 1, 1), -1, 0)
    return skv - skv.T

def Rmatrix(v, u):
    vu = np.cross(v, u)
    skm = skew_matrix(vu)
    cosvu = np.dot(v,u)
    sfactor = 1./(1 + cosvu)
    return np.identity(3) + skm + sfactor * np.dot(skm, skm)


#Matrix3d Rmatrix(Vector3d v, Vector3d u)
#{
  #Vector3d vu = v.cross(u);
  #Matrix3d skm = skew_matrix(vu);
  #double cosvu = v.dot(u);
  #//cout << "Rmatrix cosvu == " << cosvu << endl;
  #if (cosvu == -1) {
    #cout << "Rmatrix with cosvu == -1" << endl;
    #return rx180; // can't do this
  #}
  #double sfactor = 1/(1 + cosvu);
  #return Matrix3d::I() + skm + sfactor * skm*skm;
#}

def spherical(theta, phi):
    return np.array([sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)])
  
import tqdm
e_z = np.array([0,0,1.])
def check_persistance():
    cfdtheta = read_configuration_npy()[-1]
    print(cfdtheta.shape)
    cphi = np.random.rand(*cfdtheta.shape)*(2*np.pi)
    cthetaij = np.zeros_like(cphi)
    for i in tqdm.tqdm(list(range(cfdtheta.shape[0]))[:10000]):
        M = np.identity(3) # keep track of rotation matrix
        s_ax = e_z
        for j in range(cfdtheta.shape[1]):
            theta = cfdtheta[i,j]; phi = cphi[i,j]
            ax = spherical(theta, phi)
            M_dash = Rmatrix(e_z, ax)
            M = np.matmul(M_dash, M)
            tj = np.dot(M, e_z)
            thetaj = np.dot(tj, e_z)
            cthetaij[i,j] = thetaj

    corr = np.mean(cthetaij, axis=0)
    x = np.arange(1, len(corr)+1)
    y = -np.log(corr)
    p = np.polyfit(x, y, 1)
    #print 1./np.mean(np.gradient(y))
    print(1./np.mean((y[-1]-y[0])/(x[-1]-x[0])))
    #print p[0]*1. + p[1]
    print("P", p)
    plt.plot(x, y)
    plt.xlabel('contour length')
    plt.show()

    #cfdtheta[:,s:] - cfdtheta[:-s]

##

check_persistance()

##

def check_rmatrix():
    e_x = np.array([1.,0,0])
    exz = unit(e_z + e_x)
    M = Rmatrix(e_z, exz)
    exz_dash = np.dot(M, e_z)
    print(exz)
    print(exz_dash)

#check_rmatrix()
