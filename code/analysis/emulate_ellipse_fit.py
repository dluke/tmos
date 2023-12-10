# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.6.0
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# we need to be careful in our definition  of leading and trailing poles since the 
# experimental ellipse fitting method is totally different from the perfect microscope we have in simulation
# also the simulated poles are generally at the centers of the spherocylindrical ends, not at the boundary of the shape

# %% 
import sys, os
import copy
join = lambda *x: os.path.abspath(os.path.join(*x))
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import seaborn as sns
import scipy.stats

import command
import readtrack
import pili
import parameters
import _fj
import fjanalysis
import twanalyse
import rtw
import twutils
import shapeplot

# %% 
target_dir = twanalyse.get_walking_target_dir()
trs = twanalyse.load_walking_target()
# %% 
# perform the operation
old_trs = trs
new_trs = [tr.extend_projected_axis_by_radius() for tr in trs]

# %% 

lin_old_trs = [_fj.linearize(tr) for tr in old_trs]
lin_new_trs = [_fj.linearize(tr) for tr in new_trs]

fig, ax = plt.subplots(figsize=(5,5))

ex = 0
shapeplot.longtracks(ax, [old_trs[ex]], linekw={})
shapeplot.longtracks(ax, [new_trs[ex]], linekw={})
# ^^ just checking it did something

# %% 
# compute summary statistics 
with command.chdir(target_dir):
    oldld = twanalyse._summary(old_trs)
    newld = twanalyse._summary(new_trs)

# %% 
# what is the difference?

_interest = ['lvel.mean', 'deviation.var', 'qhat.estimate', 'ahat.estimate']
def describe(ld, _interest):
    return {key: twutils.make_get(key)(ld) for key in _interest}
print('old')
print(describe(oldld, _interest))
print('new')
print(describe(newld, _interest))

# %% [markdown]
# the difference is extremely significant 
# we are going to have to bifurcate our analysis again the recompute all the summary statistics with this new procedure



