
import unittest
import numpy as np


import tmos
from tmos.base import Vector3d
import tmos.pili 

import parameters
import ctfp3d as tfp

plist = ["extension", "retraction", "ret_on", "ret_off", "ext_on", "ext_off", "resample", "detach", "dissolve", "spawn"]


FLOAT_TOL = 1e-6
def vectorequals(u, v):
    return (u-v).len() < FLOAT_TOL

def to_degrees(rad):
    return rad*180/np.pi

class TestPilus(unittest.TestCase):

    def setUp(self):
        self.args = parameters.read()
        self.args.apply()

        cell = tfp.setup_cell(self.args)
        ret_motor, ext_motor = 0, 1
        anchor = Vector3d(0,0,1.5)
        axis = Vector3d(0,0,1)
        pilus = tmos.pili.PiliWLC3d(0, 1.0, 1.0,
            anchor, axis, ret_motor, ext_motor)
        pilus.isbound = True
        cell.add_pilus(pilus)
        cell.update_anchors()

        self.cell = cell
        self.pilus = pilus

        self._set_attachment()
        pilus.isbound = True

    def _set_attachment(self):
        x = np.sqrt(self.pilus.leq**2 - self.pilus.lab_anchor.z**2)
        self.pilus.attachment = Vector3d(self.pilus.lab_anchor.x + x, 0, 0)
        pv = (self.pilus.attachment - self.pilus.lab_anchor).unit()
        self.assertAlmostEqual(pv.len() , 1.0)
        self.assertAlmostEqual(self.pilus.lab_axisEq.len() , 1.0)
        self.anchor_angle = np.arccos(pv.dot(self.pilus.lab_axisEq))

    def _set_angle(self, angle):
        self.pilus.leq = self.pilus.lab_anchor.z / np.cos(np.pi/2 - angle)
        self._set_attachment()
        self.assertAlmostEqual(angle, self.anchor_angle)
        
    def test_update_anchors(self):
        pilus = self.pilus
        cell = self.cell
        self.assertTrue(vectorequals(cell.body.frame.to_lab_rt(pilus.axisEq), pilus.lab_axisEq))
        self.assertTrue(vectorequals(cell.body.frame.to_lab(pilus.anchor), pilus.lab_anchor))

    def test_attachment(self):
        self.assertTrue(self.pilus.attachment.x > self.pilus.lab_anchor.x)
        self.assertTrue(self.pilus.attachment.z == 0)
        
    def test_anchor_angle(self):
        self.pilus.leq = self.pilus.lab_anchor.z
        self._set_attachment()
        self.assertAlmostEqual(self.anchor_angle, np.pi/2)
        self.assertAlmostEqual(self.anchor_angle, self.pilus.get_anchor_angle())
        self.pilus.leq = 1e38
        self._set_attachment()
        self.assertAlmostEqual(self.anchor_angle, 0)
        self.assertAlmostEqual(self.anchor_angle, self.pilus.get_anchor_angle())
        self.pilus.leq = 2
        self.pilus.lab_anchor = Vector3d(0,0,np.sqrt(3))
        self._set_attachment()
        self.assertAlmostEqual(self.anchor_angle, np.pi/3)
        self.assertAlmostEqual(self.anchor_angle, self.pilus.get_anchor_angle())

    def test_bound_retract_rate(self):
        pilus = self.pilus
        kb_sh = self.args.Pili.kb_sh
        self._set_angle(np.pi/4)
        rate = pilus.get_bound_retract_rate()
        self.assertEqual(rate, kb_sh)
        # set the threshold to pi/6
        self.args.Pili.anchor_angle_threshold = np.pi/6
        self.assertEqual(pilus.get_bound_retract_rate(), 0.0)
        self._set_angle(np.pi/12)
        self.assertEqual(pilus.get_bound_retract_rate(), kb_sh)

       
class TestRates(unittest.TestCase):

    def setUp(self):
        self.args = parameters.read()

        ret_motor, ext_motor = 0, 1
        taut_pilus = tmos.pili.PiliWLC3d(0, 1.0, 1.0 - 0.004,
            Vector3d(), Vector3d(1,0,0), ret_motor, ext_motor)
        taut_pilus.isbound = True
        self.pilus = taut_pilus

    def test_surface_sensing_talaI(self):
        # talaI
        pilus = self.pilus
        tmos.pili.Pili.allow_bound_ext_motor = False
        tmos.pili.Pili.force_bound_retraction = True
        pilus._attach_upkeep(0)
        self.assertTrue(pilus.ret_motor == 1)
        self.assertTrue(pilus.ext_motor == 0)
        rates = pilus.get_rates()
        self.assertTrue( rates[plist.index("ext_on")] == 0)
        self.assertTrue( rates[plist.index("ret_on")] == 0)

    def test_surface_sensing_talaII(self):
        # talaII
        pilus = self.pilus
        tmos.pili.Pili.allow_bound_ext_motor = False
        tmos.pili.Pili.force_bound_retraction = False
        pilus._attach_upkeep(0)
        self.assertTrue(pilus.ret_motor == 0)
        self.assertTrue(pilus.ext_motor == 0)
        rates = pilus.get_rates()
        self.assertTrue( rates[plist.index("ext_on")] == 0)
        self.assertTrue( rates[plist.index("ret_on")] > 0)

    def test_surface_sensing_koch(self):
        # Koch
        pilus = self.pilus
        tmos.pili.Pili.allow_bound_ext_motor = True 
        tmos.pili.Pili.force_bound_retraction = False
        pilus._attach_upkeep(0)
        self.assertTrue(pilus.ret_motor == 0)
        self.assertTrue(pilus.ext_motor == 1)
        # set no motor
        pilus.ext_motor = 0
        # check that either motor can bind
        rates = pilus.get_rates()
        self.assertTrue( rates[plist.index("ext_on")] > 0)
        self.assertTrue( rates[plist.index("ret_on")] > 0)


    def test_allow_bound_extension(self):
        ret_motor, ext_motor = 0, 1
        taut_pilus = tmos.pili.PiliWLC3d(0, 1.0, 1.0 - 0.004,
            Vector3d(), Vector3d(), ret_motor, ext_motor)
        bent_pilus = tmos.pili.PiliWLC3d(0, 1.0, 1.0 + 0.004,
            Vector3d(), Vector3d(), ret_motor, ext_motor)
        taut_pilus.isbound = True
        bent_pilus.isbound = True

        tmos.pili.Pili.allow_bound_extension = True
        self.assertTrue(self.args.Pili.allow_bound_extension == True)
        rates = taut_pilus.get_rates()
        ext_rate, ret_rate = rates[0], rates[1]
        self.assertTrue(taut_pilus.istaut() == True)
        self.assertTrue(ext_rate > 0)
        self.assertTrue(ret_rate == 0)
        
        rates = bent_pilus.get_rates()
        ext_rate, ret_rate = rates[0], rates[1]
        self.assertTrue(bent_pilus.istaut() == False)
        self.assertTrue(ext_rate == 0)

        tmos.pili.Pili.allow_bound_extension = False
        self.assertTrue(self.args.Pili.allow_bound_extension == False)

        # taut
        rates = taut_pilus.get_rates()
        ext_rate, ret_rate = rates[0], rates[1]
        self.assertTrue(ext_rate == 0)

        # not taut
        rates = bent_pilus.get_rates()
        ext_rate, ret_rate = rates[0], rates[1]
        self.assertTrue(ext_rate == 0)





if __name__ == '__main__':
    unittest.main()

