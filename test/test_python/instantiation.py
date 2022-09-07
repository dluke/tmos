
# python script to test initialisation and __str__ function for 
# all my classes.

# can define a custom target to run this and add it to cmake

# We would like to visually inspect the __str__ output here 
# python unittests will be in another file and target.

import tmos
import tmos.base as base
import tmos.surface as surface 
import tmos.pili as pili
import tmos.mdynamics as md
import tmos.vtkwriter as vtkwriter

# now name all the classes
from tmos.base import Vector3d,Frame
from tmos.surface import Plane, Capsule

#...


# test
def test_init():

    print(Vector3d())
    print(Frame())
    pl = Plane()
    print(pl)
    cap = Capsule(Vector3d(), Vector3d(1,0,0), 0.5, 2.)
    print(cap)
    cell = pili.Cell3d(0, 4, cap, pl)
    cell.common_init()
    print(cell)





if __name__=='__main__':

    test_init()

