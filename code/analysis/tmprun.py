
import numpy as np

import shapeplot
import analysefj
import readyr
import command

from matplotlib import pyplot as plt

#yr = readyr.init()
#yr.compute_ax()
#analysefj.plot_corr_orient(yr)


@command.defaultsave()
def fjcapsdraw(fj):
    ax = plt.gca()
    shapeplot.exp_capsdraw(ax, fj)
#

def test_corr():
    fj = analysefj.init(first=False,alldirs=False)
    fj.filter_short(1000.)
    fj.override_orientation_using_lead_trail()
    fj.compute_ax()
    analysefj.plot_corr_orient(fj)
    fjcapsdraw(fj)

def test_orient():
    fj = analysefj.init(first=False,alldirs=False)
    fj.filter_short(1000.)
    fj.override_orientation_using_lead_trail()
    thetas = fj.get_whole_col('orientation')
    print(np.max(thetas))
    print(np.min(thetas))
    fj.compute_ax()
    analysefj.plot_frmodes(fj, 95., 1.)


if __name__=='__main__':

    test_corr()
    #test_orient()

