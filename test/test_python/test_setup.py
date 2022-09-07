
import runtfp_global
import ctfp3d 

import tmos
from tmos.base import Vector3d


def test_setup_cell():
    args = runtfp_global.setup()
    #print args
    idx = 0
    cell = ctfp3d.setup_cell(args, idx)
    print(cell.pili[0].n)
    return cell



def test_setup_pilus():
    args = runtfp_global.setup()
    kpgen = ctfp3d.setup_kpgen(args)
    pilus = tmos.pili.PiliWLC3d(0, 1, 1-0.004, Vector3d(), Vector3d(1,0,0), 0, 0, kpgen)
    print(pilus.force_l())

# cell = test_setup_cell()

test_setup_pilus()


