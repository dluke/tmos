
import numpy as np
import parameters
import ctfp3d as tfp

import tmos
from tmos.base import Vector3d
import tmos.pili 

import matplotlib.pyplot as plt

plt.rcParams.update({
    'text.usetex': True,
    'axes.labelsize': 30,
    'xtick.labelsize': 30,
    'ytick.labelsize': 30
    })


class TestAnchor():

    def setUp(self):
        self.args = parameters.read()
        self.args.apply()

        cell = tfp.setup_cell(self.args)
        ret_motor, ext_motor = 1, 0
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
        self.anchor_angle = np.arccos(pv.dot(self.pilus.lab_axisEq))

    def _set_angle(self, angle):
        self.pilus.leq = self.pilus.lab_anchor.z / np.cos(np.pi/2 - angle)
        self._set_attachment()

    def test(self):
        self.args.Pili.anchor_angle_threshold = np.pi/6
        self.args.Pili.anchor_angle_smoothing_fraction = 1.0
        self.args.Pili.enforce_stalling = False

        def plot_rate(out='retraction_rate.png'):
            fig, ax = plt.subplots(figsize=(5,5))
            aspace = np.linspace(0,np.pi/2,180)
            ret = np.empty_like(aspace)
            for i, angle in enumerate(aspace):
                self._set_angle(angle)
                rates = self.pilus.get_rates() 
                retrate = rates[1]
                ret[i] = retrate * 0.004 * 1.5
            ax.plot(aspace, ret, linewidth=6)
            ax.set_xlabel(r"$\theta_a$", fontsize=30)
            ax.set_ylabel(r"$v_{\mathrm{ret}}$")
            ax.set_xticks([0, 0.6*np.pi/2, np.pi/2])
            ax.set_xticklabels([0, r"$\alpha$", r"$\pi/2$"])
            print('save ', out)
            plt.tight_layout()
            plt.savefig(out, transparent=True)
        plot_rate('retraction_rate_00.svg')

        self.args.Pili.anchor_angle_smoothing_fraction = 0.2
        plot_rate('retraction_rate_01.svg')

        self.args.Pili.anchor_angle_smoothing_fraction = 0.6
        self.args.Pili.anchor_angle_threshold = 0.0
        plot_rate('retraction_rate_02.svg')

if __name__=='__main__':

    test = TestAnchor()
    test.setUp()
    test.test()

 