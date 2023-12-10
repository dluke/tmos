
import time

import tmos
from tmos.base import Vector3d
import tmos.surface as surface
import tmos.pili as pili
import ctfp3d

import tfputils

vdummy = Vector3d()

# The problem is that python data holding objects are not being kept alive correctly
def test_wlc():

    cls = tfputils.read_configuration_npy()
    wlcg = pili.WLCgenerator(1.,0.1,len(cls))
    for i, cdata in list(enumerate(cls)):
        #cdata = np.array(cdata, dtype='float64')
        #wlcg.setnth(i, cdata)
        wlcg.append(cdata)
        #print i+1
        #print wlcg.sample(i+1, vdummy)

    return wlcg

def sample_forever(wlcg):
    while True:
        print(wlcg.sample(1, vdummy))
        time.sleep(1)

def step_forever(cell):
    print('step forever')
    while True:
        print(time.time())
        pili.kmcstep(cell, 0.)
        time.sleep(1)


def test_return_wlc():
    import scipy.constants
    kb = scipy.constants.Boltzmann
    T = 273.0 + 30.0 # absolute temperature
    kbT = (kb * T) * 10**18. 
    Lp = 5. # persistence length
    ka = Lp * kbT
    a = 1.0
    nmax = 20
    kpgen = pili.KPgenerator(ka, a, nmax, T)

    pl = surface.Plane()
    cap = surface.Capsule(Vector3d(0,0,0.495), Vector3d(1,0,0), 0.5, 2.)
    cell = pili.CellWLC3d(0, 0, cap, pl, kpgen)
    cell.pilivar = 4.0
    cell.common_init()
    pilus = cell.spawn_pilus()
    cell.add_pilus(cell.spawn_pilus())
    return cell

def test_return_cell():

    cell = test_return_wlc()
    print(cell)
    #wlcg = cell.get_wlc()

    kmc = pili.Kmc(0.)
    kmc.kmcstep(cell, 0.)


def main():
    #wlcg = ctfp3d.initialise_pili_sampler()
    #print 'out scope sample'
    #print wlcg.sample(0)

    cls = tfputils.read_configuration_npy()
    wlcg = pili.WLCgenerator(1.,0.1,len(cls))
    for i, cdata in list(enumerate(cls))[:5]:
        #cdata = np.array(cdata, dtype='float64')
        #wlcg.setnth(i, cdata)
        wlcg.append(cdata)
    ch = wlcg.sample(1, Vector3d(1.,0.,0.))
    print(ch)

    pl = surface.Plane()
    cap = surface.Capsule(Vector3d(), Vector3d(1,0,0), 0.5, 2.)
    cell = pili.CellWLC3d(0, 8, cap, pl, wlcg)
    #print cell


if __name__=='__main__':
    #main()

    #test_wlc()
    #test_return_wlc()
    test_return_cell()
    
