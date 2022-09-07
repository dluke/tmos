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
# We played around with the concept of "actions" before but we considered them synonymous with steps
# in the linearised trajectory
# Here we will redefine actions and examine their properties across fanjin and simulated data

# %% 
import warnings
import sys, os
import copy
join = lambda *x: os.path.abspath(os.path.join(*x))
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import pandas as pd

import parameters
import stats
import command
import twutils
import eventanalyse
import readtrack
import collections
import _fj
norm = np.linalg.norm

# %% 
# load subsets
subsets = _fj.load_subsets()
subset = subsets["top"]

# %% 
# what is the distribution of Angle(v_i, v_{i+1})
threshold = 10 # degrees

angles = 180/np.pi * np.concatenate([tr.step_angle() for tr in subset])
ax = sns.histplot(angles)
ax.set_xlabel("angle (degrees)")
ax.axvline(threshold, c='k', linestyle='--')

# %% 
# arbitrary threshold
arb = threshold * np.pi/180
print('arbitrary  threshold %f radians' % arb)
# take linearised  tracking data and return an array with action indices seperating actions

# df = ltr.decompose_actions(arbt=arb)
df = pd.concat([ltr.decompose_actions(arbt=arb, track_idx=i) for i, ltr in enumerate(subset)])
# %%

def describe(df):
    print("median velocities for actions")
    medv = df["velocity"].median()
    medv_filter = df[df["large"] == True]["velocity"].median()
    print("(all, large) = ({}, {})".format(medv, medv_filter))
    #
    print("mean velocities for actions")
    mv = df["velocity"].mean()
    mv_filter = df[df["large"] == True]["velocity"].mean()
    print("(all, large) = ({}, {})".format(mv, mv_filter))
describe(df)

# g = sns.displot(data=df, x="velocity", hue="large")
# g.ax.set_xlim((0, 1))
# fig, ax = plt.subplots(figsize=(5,5))
# sns.histplot(df["velocity"], ax=ax)
# sns.histplot(df[df["large"]==True]["velocity"], palette='red', ax=ax)

# %%
candidate = subsets["candidate"][0]
cdf = candidate.decompose_actions(arbt=arb)
describe(cdf)

# %%
candidate = subsets["candidate"][0]
vdf = candidate.decompose_actions_by_velocity()
vdf
vdf[vdf["fast"] == True]["angle"].mean()
vdf[vdf["fast"] == False]["angle"].mean()


