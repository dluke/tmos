

import sys
import numpy as np
import tmos 
import tmos.base
Vector3d = tmos.base.Vector3d
import ctfp3d
import parameters

norm = np.linalg.norm
sin = np.sin
cos = np.cos
pi = np.pi
R = 0.5

# Vector3d to array
def toarray(v):
    return np.array([v.x,v.y,v.z])

def setup_cell(pilivar=R * np.pi/2., polarisation=0.):
    args = parameters.read().apply()
    args.cell.length = 2
    args.cell.initaxis = (0,0,-1)
    args.cell.initcenter = (0,0,1.5)

    args.surface.shape = 'plane'
    args.Cell3d.repulsive_only = True
    args.ACell.pilivar = pilivar
    args.ACell.polarisation = polarisation
    cell = ctfp3d.setup_cell(args)
    return cell

def run_minimise(cell):
    # python wrapper on 
    def fun(state):
        # return double
        return cell.state_energy(list(state))
    def jac(state):
        # return array
        return np.array(cell.state_gradient(list(state)))

    alg, jac = 'L-BFGS-B', jac
    res = scipy.optimize.minimize(
            fun,
            cell.get_state(),
            jac=jac,
            method=alg,
            options={'ftol':1.}
            )
    print(res)
    if res.success == False:
        print("minimisation algorithm failed")
    cell.set_state(list(res.x))


def test_set_state():
    cell = setup_cell()
    delta = 0.1
    #newstate = [0.,0.,0.,0.,0.,0.]
    newstate = [0.,0.,1.5-delta,0,-1.10*np.pi,0]
    #newstate = [0.,0.,1.5,0,+*np.pi,0]

    #newstate = [2.13013196 , 0.47833595,  1.82216012, -0.00536568,  1.05884445, -0.20632684]
    #newstate = [-0.413332373061, 1.13978207129, -2.07969498611, 0.486412256704, -0.0171063618926, -0.687785625558]
    print("newstate ", newstate)

    cell.set_state(newstate)
    #print "surface energy", cell.energy_surface()
    #print cell.get_headpt()
    #print cell.get_trail()

    print('cell 1')
    #print 'axis', cell.get_axis()
    #print 'head/tail', cell.get_headpt(), cell.get_trail()
    print('grad', cell.grad())
    state = cell.get_state()
    print('get_state', state)
    cell.set_state(state)
    print('cell again')
    #print 'axis', cell.get_axis()
    #print 'head/tail', cell.get_headpt(), cell.get_trail()
    print('grad', cell.grad())

#test_set_state()

def test_surface_grad():
    cell = setup_cell()
    delta = 0.001
    sstate = np.array([0.,0.,0.5,0.,-pi/2.+delta,0.])
    cell.set_state(sstate)
    print('energy', cell.energy())
    print('sgrad', cell.surface_grad())
    assert(np.array_equal(cell.surface_grad(),  cell.grad()))
    assert(np.array_equal(cell.pili_grad(), np.zeros(6)))

#test_surface_grad()


def spherical(theta, phi):
    return Vector3d(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta))

# settings
theta = pi/4.
length = 2
R = 0.5

def setup_pilus(cell):
    pilus = cell.pili[0]
    lab_anchor = cell.get_lab_anchor(pilus)
    sph = Vector3d(0,0,length/2.) + 0.5 * (0.46/0.5) * spherical(theta, 0)
    pilus.anchor = sph
    pilus.axisEq = spherical(theta, 0)
    attach = Vector3d(lab_anchor.x+0.5, 0., 0.)
    pilus.attachment = attach 
    pilus.isbound = True
    pilus.pl = (pilus.attachment - lab_anchor).len()
    pilus.leq = pilus.pl - tmos.pili.Pili.d_free
    return pilus

sstate = [-1.,0.,.5,0.,pi/2.,0.]
def test_pili_grad():
    cell = setup_cell()
    # state with r_head sitting on the surface at (x,z) = (0.,0.5)
    cell.set_state(sstate)
    tol = 1e-10
    assert(cell.energy() < tol)
    pilus = setup_pilus(cell)
    lab_anchor = cell.get_lab_anchor(pilus)

    axis = pilus.attachment - lab_anchor
    print(cell.energy())
    print(cell.pili_grad())

def main_test_gradient():
    #test_set_state()
    print('surface grad')
    test_surface_grad()
    print('pili grad')
    test_pili_grad()

#main_test_gradient()
test_set_state()

#################################################################################
# 
# d_a
"""
cell = setup_cell()
pilus = setup_pilus(cell, lab_anchor)
cell.set_state(sstate)
d_a0 = toarray(pilus.anchor)
lab_anchor = cell.get_lab_anchor(pilus)
r_b = toarray(pilus.attachment)

# results
sgrad = [0.0, 0.0, 29.411980461275434, 0.0, -29.411965755286428, 0.0]
senergy = 0.0146
pgrad = [-71.83591337967758, 0.0, 25.103904814858108, 0.0, -9.903425543718349, 0.0]
penergy = 0.15219

import algorithm.axis_angle as axa
dot = np.dot

# reconstruct the calculation
def Rkforp(p):
    theta = norm(p)
    assert(theta != 0) # tmp
    ptilde = axa.tilde(p)
    p1, p2, p3 = p
    ptilde1 = axa.skew([1-p1**2/theta**2, -p1*p2/theta**2, -p1*p3/theta**2 ])/theta
    ptilde2 = axa.skew([ -p2*p1/theta**2, 1-p2**2/theta**2, -p2*p3/theta**2 ])/theta
    ptilde3 = axa.skew([ -p3*p1/theta**2, -p3*p2/theta**2, 1 - p3**2/theta**2])/theta
    ptildek = [ptilde1, ptilde2, ptilde3]
    def Rk(k):
        assert k in [0,1,2]
        return (p[k]*sin(theta)/theta * np.dot(ptilde,ptilde) 
                + (1-cos(theta))*(dot(ptildek[k],ptilde) + dot(ptilde,ptildek[k]))
                + p[k]*cos(theta)/theta * ptilde + sin(theta) * ptildek[k]
                )
    return Rk

ks = 1e4
def pgrad(rab, leq):
    if leq > rab:
        return 0.
    else:
        return (rab - leq)/leq

def dudpk():
    Rk = Rkforp(np.array([0,pi/2.,0],dtype=float))


"""
