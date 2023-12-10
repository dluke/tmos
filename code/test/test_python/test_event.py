

import runtfp_global
import ctfp3d as tfp

import tmos.pili as pili

def test_null_event():

    #args = runtfp_global.setup()
    #idx = 0
    #cell = tfp.setup_cell(args, idx)
    kmc = pili.Kmc(0.)
    for ev in kmc.events:
        print(ev)

test_null_event()

