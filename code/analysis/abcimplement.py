

# %%
import warnings
import copy
import collections
import numpy as np
import pandas as pd
import sobol
import scipy.stats
import scipy.spatial
# import matplotlib.pyplot as plt

from matplotlib.ticker import FormatStrFormatter
import sigfig as sf
import matplotlib as mpl

# implementation  and utilities for approximate bayesian computation

# uniform sampling for problem
def mcsampling(problem, N):
    # sample a flat prior in the space defined by problem
    d = problem["num_vars"]
    rv = np.random.rand(N,d)
    for i, bounds in enumerate(problem["bounds"]):
        xn, xm = bounds
        rv[:,i] = rv[:,i] * (xm-xn) + xn
    return rv

def filter_missing(lookup, nan_idx):
    new = copy.deepcopy(lookup)
    uids = [new[0][idx] for idx in nan_idx]
    for uid in uids:
        del new[1][uid]
    for idx in sorted(nan_idx, reverse=True):
        del new[0][idx]
    return new 

def load_problem_simulation(psd, local_stats_path='local.json'): 
    # problem simulation data
    _lookup = sobol.read_lookup(psd["simdir"])
    _problem = sobol.read_problem(psd["simdir"])
    psd["M"] = len(_lookup[0])
    # check failed
    types= {"failed_condition": object}
    Y = sobol.collect(["failed", "failed_condition"], targetdir=psd["simdir"], tdct=types, local_stats_path=local_stats_path)
    Y["failed"][np.isnan(Y["failed"])] = False
    failed = Y["failed"].astype(bool)
    # assert( all(Y["failed_condition"][failed] == 'step_condition') )
    if 'objectives' in psd: 
        # second call to sobol.collect is inefficient
        Ym = sobol.collect(psd['objectives'], targetdir=psd["simdir"], local_stats_path=local_stats_path)
        nans = {}
        for name in Ym.keys():
            nan = np.isnan(Ym[name])
            if any(nan):
                print("nan found in {}. filtering {} samples".format(name, int(np.sum(nan))))
            nans[name] = nan
            failed[nan] = True
        # nan = np.logical_and.reduce(nans.values())
    failed_count = collections.Counter(Y["failed_condition"][failed])
    print("failed:", failed_count)
    nan_idx = np.argwhere(failed).ravel()
    N = len(_lookup[0])
    _lookup = filter_missing(_lookup, nan_idx)
    print("filtered out {}/{} samples".format(N - len(_lookup[0]), N))
    psd["n_failed"] = len(nan_idx)
    psd["lookup"] = _lookup
    psd["problem"] = _problem
    _lduid = sobol.collect_lduid(psd["simdir"], local_stats_path)
    psd["lduid"] = _lduid # backwards compatible. store lduid separate from MCData object.
    psd["data"] = sobol.MCData(_problem, _lookup, _lduid)
    print("loaded data from ", psd["simdir"])
    return psd

# brute force mc sampling of arbitrary function
def force_sample_target(target_density, bounds, size=100000):
    # use with scipy.interpolate.RectBivariateSpline(ydata,xdata,zdata)
    xlim, ylim, zlim = bounds
    # sample uniformly in 3d box
    xt = np.random.uniform(*xlim, size=size)
    yt = np.random.uniform(*ylim, size=size)
    zt = np.random.uniform(*zlim, size=size)
    # evaluate pdf at xt,yt 
    teval = np.array([target_density(yt[i], xt[i]) for i in range(size)])
    # accept if zt < target_density(xt,yt)
    accept = np.less(zt, teval.ravel())
    return xt[accept], yt[accept]

# utilities for computing kernel density

# https://stats.stackexchange.com/questions/450570/kernel-density-estimation-and-boundary-bias
def mirror(arr, bounds):
    pmin, pmax = bounds
    left = 2*pmin - arr 
    right = 2*pmax - arr
    return np.concatenate([left,arr,right])

# mirror 2d points in y boundaries
def mirrorpts(data, lims):
    data[1] = mirror(data[1], lims)
    data[0] = np.concatenate([data[0],data[0],data[0]])
    return data

protected = ['pdks']
def _is_simscore(name):
    return (name.startswith('fanjin') or name.startswith('ks_') or name.startswith('sim')
        or name in protected)

# make this resistant to nan
def regularise_stats(statdf, reference, objectives):
    newdf = statdf.copy()
    newref = reference.copy()
    std = {}
    for objective in objectives:
        stat = newdf[objective].to_numpy()
        if _is_simscore(objective): 
            # create a score starting at 0.0 by subtracting minimum

            stat = stat - np.nanmin(stat)
            _std = np.nanstd(stat)
            newdf[objective] = stat/_std
        else:

            _std = np.nanstd(stat)
            newref[objective] = reference[objective]/_std
            newdf[objective] = stat/_std
        std[objective] =  _std
    newdf.attrs['std']= std
    return newdf, newref

# simple rejection ABC helps us test more advanced algorithms
# generalise rejection_abc to any number of metrics
# handle similarity scores more naturally?

def rejection_abc(statdf, objectives, reference, N, ss_component=False, norm=2):
    # statdf is dataframe with parameters + only the objectives that we are using for abc 
    refpt = np.zeros(len(objectives))
    for i, objective in enumerate(objectives):
        # to make l2 norm work for similarity score we set reference value to 0.
        refpt[i] = 0. if _is_simscore(objective) else reference[objective]
    stats = np.stack([statdf[objective] for objective in objectives],axis=1)
    stats_diff = stats - refpt

    # score = np.linalg.norm(stats_diff, axis=1)
    score = np.linalg.norm(stats_diff, ord=norm, axis=1)

    statdf["score"] = score
    if ss_component:
        for i in range(len(objectives)):
            statdf["s"+str(i)] = stats_diff[:,i]
    sortdf = statdf.sort_values("score")
    accepted = sortdf.iloc[:N]
    delta = accepted["score"].to_numpy()[-1]
    # print(objectives)
    # print("N = {}, delta = {}, target = {}".format(N, delta, refpt))
    accepted.attrs['objectives'] = objectives
    return accepted

# compute
def mean_accepted(problem, accepted):
    _names = problem['names']
    
    epsilon = accepted['score'].iloc[-1]
    kernel = new_epanechnikov(epsilon)
    scores = accepted['score'].to_numpy()
    ws = kernel(scores)
    # print('weights {} ... {}'.format(ws[:3], ws[-3:]))

    # weighted mean
    mean = {name : np.average(accepted[name], weights=ws) for name in _names}
    # weighted variance
    variance = {name : np.average((accepted[name]-mean[name])**2, weights=ws) for name in _names}
    return mean, variance


def maximum_likelihood(problem, accepted):
    # * even though we have an explicit approximate posterior distribution
    # * finding the maximum is challenging because of the, high dimensionality
    # * few samples and the need for a smoothing kernel
    # * we have so far avoided constructing the 4d distribution and instead used 2d projections
    pass # see tabulate_inference


# ----------------------------------------------------------------
# implementing regression ABC

def new_epanechnikov(delta):
    # should be like numpy ufunc
    def epanechnikov(x):
        return max(3*(1 - (x/delta)**2)/(4*delta), 0)
        # return max(3*(1 - (x/delta)**2)/4, 0)
    return np.vectorize(epanechnikov)


inv = np.linalg.inv
def llregression(_acc, objectives, statref, projection):
    """
    _acc = accepted samples
    """
    print("using local linear regression")
    delta = _acc["score"].iloc[-1]
    print("max accepted", delta)
    _acc_rho = _acc["score"].to_numpy()
    epanechnikov = new_epanechnikov(delta)
    weight = np.array([epanechnikov(x) for x in _acc_rho])
    N = len(_acc)
    par = np.column_stack([_acc[p] for p in projection])

    xsls = []
    for name in objectives:
        if _is_simscore(name):
            xsls.append(_acc[name])
        else:
            xsls.append(_acc[name] - statref[name][0])
    Xs = np.column_stack(xsls)
    # Xs = np.column_stack([_acc[name] - statref[name] for name in objectives if not _is_simscore(name) else _acc[name]])
    X = np.column_stack([np.ones(N), Xs])
    W = np.diag(weight)
    coef = inv(X.T @ W @ X) @ X.T @ W @ par
    beta = coef[1:,:]
    llpar = par - X[:,1:] @ beta
    llacc = _acc.copy()
    for i, p in enumerate(projection):
        llacc[p] = llpar[:,i]
    llacc["weight"] = weight
    return llacc


# ----------------------------------------------------------------
# implementing regression ABC
# plotting utilities
import itertools
import matplotlib.pyplot as plt
import seaborn as sns

# dep this method
def plot_accepted_projection(problem, accepted):
    fig, axes = plt.subplots(1,3, figsize=(12,4))
    _names = problem["names"]
    _bounds = problem["bounds"]
    _cidx = itertools.combinations(range(3),2)
    for ax_i, idx in enumerate(_cidx):
        ax = axes[ax_i]
        _p1, _p2 = [_names[_i] for _i in idx]
        print(idx, _p1, _p2)
        _data = accepted[[_p1, _p2, "score"]]
        sns.scatterplot(data=_data, x=_p1, y=_p2, hue="score",  ax=ax, legend=False)
        # sns.kdeplot(data=_data, x=_p1,  y=_p2, fill=True, ax=ax)
        ax.set_xlabel(_p1)
        ax.set_ylabel(_p2)
        xlim, ylim = [_bounds[_i] for _i in idx]
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
    return fig, axes


def problemplot3d(problem, df3d, hue=None, snskw={}):
    fig, axes = plt.subplots(1,3, figsize=(12,4))
    _names = problem["names"]
    _bounds = problem["bounds"]
    _cidx = itertools.combinations(range(3),2)
    for ax_i, idx in enumerate(_cidx):
        ax = axes[ax_i]
        _p1, _p2 = [_names[_i] for _i in idx]
        print(idx, _p1, _p2)
        sns.scatterplot(data=df3d, x=_p1, y=_p2, hue=hue, ax=ax,  **snskw)
        ax.set_xlabel(_p1)
        ax.set_ylabel(_p2)
        xlim, ylim = [_bounds[_i] for _i in idx]
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
    plt.tight_layout()
    return fig, axes

 
def problemplot4d(problem, accepted, objective=None):
    # print(objective)
    _names = problem["names"]
    _bounds = problem["bounds"]
    n = len(_names)
    hue = None
    is_lvel = objective.startswith("fanjin")
    if is_lvel:
        data  = accepted[_names+[objective]] 
        hue = objective
    else:
        data  = accepted[_names+["score"]] 
        hue = "score" 
    pkws = {"alpha": 0.8}
    g = sns.pairplot(data, hue=hue, plot_kws=pkws)
    g.fig.suptitle(objective, y=1.08) # y= some height>1
    for i in range(n):
        for j in range(n):
            if i==j:
                continue
            _xlim = _bounds[j]
            _ylim = _bounds[i]
            g.axes[i,j].set_xlim(_xlim)
            g.axes[i,j].set_ylim(_ylim)


import statsmodels.api as sm

def _statsdens(data, x, y):
    _x = data[x].to_numpy()
    _y = data[y].to_numpy()
    dens = sm.nonparametric.KDEMultivariate(data=[_x,_y], var_type='cc', 
        bw='cv_ls')
    print('crossvalidated bandwidth', x, y, dens.bw)
    return dens

def _gkde(data, x, y, weights):
    # * no wei
    _x = data[x].to_numpy()
    _y = data[y].to_numpy()
    kde = scipy.stats.gaussian_kde(np.stack([_x,_y]), weights=weights)
    kde.set_bandwidth(kde.factor*0.8)
    return kde

def _heatmap_ax(ax, _bounds, res):
    xlim, ylim = _bounds
    ax.invert_yaxis()
    # simplify this tick locator business
    # ridx = [5,res//2,res-5-1]
    ridx = [0,res//2,res]
    ax.set_xticks(ridx)
    ax.set_yticks(ridx)
    xt= []
    yt= []
    for n in ridx:
        a = float(n)/res
        data_xa = a * (xlim[1]-xlim[0]) + xlim[0]
        xt.append( data_xa)
        data_ya = a * (ylim[1]-ylim[0]) + ylim[0]
        yt.append( data_ya)
    ax.set_xticklabels([str(round(_xt,2)) for _xt in xt])
    ax.set_yticklabels([str(round(_xt,2)) for _xt in yt])

# use the density object to  evaluate the heatmap image
def _statskdeplot(dens, bounds, ax, res=100, cmap=None):
    xlim, ylim = bounds
    xspace =  np.linspace(*xlim, num=res)
    yspace =  np.linspace(*ylim, num=res)
    X1, X2 = np.meshgrid(xspace, yspace)
    result = dens.pdf([X1.ravel(), X2.ravel()])
    _grid = result.reshape([res,res])

    if not cmap:

        # cmap = sns.color_palette("light:g", as_cmap=True)
        cmap = mpl.cm.get_cmap("Greens")
        clist = np.array([cmap(a) for a in np.linspace(0,0.8,256,True)])
        # clist[:10] = (1.0,1.0,1.0,1.0)
        cmap = mpl.colors.ListedColormap(clist)

    sns.heatmap(_grid, ax=ax, cbar=False, cmap=cmap)
    ax.invert_yaxis()
    # _heatmap_ax(ax, [xlim, ylim], res)



# pretty latex parameter  names
lpar = [r'$\tau_{\mathrm{dwell}}$', r'$\kappa$', r'$\alpha$', r'$k_{\mathrm{spawn}}$']

config = {}
config['fastdensity'] = True


def hide_ticks(ax):
    ax.xaxis.set_ticks([])
    ax.yaxis.set_ticks([])
def show_xonly(ax):
    for name in ['left', 'top', 'right']:
        ax.spines[name].set_visible(False) 
    ax.yaxis.set_ticks([])
def show_yonly(ax):
    for name in ['top', 'right', 'bottom']:
        ax.spines[name].set_visible(False) 
    ax.xaxis.set_ticks([])

def sfround(value, n):
    if value < 1:
        return sf.round(str(value), decimals=n)
    else:
        return sf.round(str(value), sigfigs=n)

def perfectplot4d(problem, accepted, simpar=None, mpar=None, lpar=lpar, _config={}):
    res = 100
    lconfig = config.copy()
    lconfig.update(_config)

    epsilon = accepted['score'].iloc[-1]
    kernel = new_epanechnikov(epsilon)
    scores = accepted['score'].to_numpy()
    weights = kernel(scores)


    # setup colormap
    marker_purple = '#9B5F9A'
    white  = '#ffffff'
    dmarkerkw = {'c':marker_purple, 'edgecolor':white, 'linewidths':2}
    xmarkerkw = {'c':marker_purple, 'edgecolor':white, 'linewidths':2}

    cmap = mpl.cm.get_cmap('Greens_r')
    clist = np.array([cmap(a) for a in np.linspace(0,0.8,256,True)])
    # clist[-10:] = (1.0,1.0,1.0,1.0)
    cmap = mpl.colors.ListedColormap(clist)

    # set context
    _context = {'text.usetex': True, 'axes.labelsize': 32, 'font.size' : 20, 'axes.facecolor':white}
    tpad = 20
    _context.update({'xtick.major.pad':tpad, 'ytick.major.pad':tpad})
    darkedge = '#546354'
    lightedge = '#C3D0C3'
    _context['axes.edgecolor'] = lightedge
    # unpack problem
    _names = problem["names"]
    _bounds = problem["bounds"]
    n = len(_names)

    # apply context
    with plt.rc_context(_context):
        #  construct the axes and handle the plotting myself
        r1, rc = 1, 32
        grkw = {'width_ratios':[r1,rc,rc,rc,rc], 'height_ratios':[rc,rc,rc,rc,r1]}
        fig, axes = plt.subplots(5,5, figsize=(20,20), gridspec_kw=grkw)
        xaxes = axes[-1,1:]
        yaxes = axes[:-1,0]
        axes[-1,0].axis('off')
        for i in range(n):
            lim = _bounds[i]
            label = lpar[i]
            xax = xaxes[i]
            yax = yaxes[i]
            xax.spines['bottom'].set_edgecolor(darkedge)
            yax.spines['left'].set_edgecolor(darkedge)
            xax.set_xlim(lim)
            yax.set_ylim(lim)
            xax.set_xlabel(label, fontsize=40, labelpad=10)
            yax.set_ylabel(label, fontsize=40, labelpad=10)
            show_xonly(xax)
            show_yonly(yax)
            _ticks = np.linspace(1.0*lim[0], 0.85*lim[1], 3, True)
            xax.set_xticks(_ticks)
            yax.set_yticks(_ticks)
            xax.set_xticklabels([sfround(_xt,2) for _xt in _ticks])
            yax.set_yticklabels([sfround(_yt,2) for _yt in _ticks])
            xax.xaxis.set_tick_params(width=4., length=10., color=darkedge)
            yax.yaxis.set_tick_params(width=4., length=10., color=darkedge)

        for y in range(n):
            for x in range(n):
                name1,  name2 = _names[y], _names[x]
                ax = axes[x,y+1]
                ax.grid(False)
                _xlim = _bounds[y]
                _ylim = _bounds[x]
                if x==y:
                    # ax.axis('off')
                    ax.set_xticks([])
                    ax.set_yticks([])
                    name = name1 
                    ins = ax.inset_axes([0.02, 0.02, 0.96, 0.96])

                    hist, bin_edges = np.histogram(accepted[name].to_numpy(), bins='auto')
                    nbins = 2*len(hist)
                    # ! note this weighting starts at 0, but the best score may be quite large
                    sns.histplot(x=accepted[name], bins=nbins, binrange=_xlim, ax=ins, color='#88CB89', weights=weights, 
                            alpha=0.6, kde=False,  line_kws={'linewidth':4})
                    # sns.histplot(x=accepted[name], binrange=_xlim, ax=ins, color='#88CB89')

                    ins.set_xlabel('')
                    ins.set_ylabel('')
                    ins.set_yticks([])
                    ins.set_xticks([])
                    ins.grid(False)
                    for spine in ax.spines.values():
                        spine.set_visible(True)
                        spine.set_edgecolor(lightedge)
                    for spine in ins.spines.values():
                        spine.set_visible(True)
                        spine.set_edgecolor(lightedge)
                if x > y:
                    if lconfig['fastdensity']:
                        dens = _gkde(accepted, x=_names[y], y=_names[x], weights=weights)
                    else:
                        dens = _statsdens(accepted, x=_names[y], y=_names[x])
                    bounds = [_xlim, _ylim]
                    _statskdeplot(dens, bounds, ax=ax)
                    def to_image_coords(dat, lim):
                        return (dat-lim[0])/(lim[1]-lim[0]) * res
                    if simpar:
                        ptx,  pty = simpar[name1], simpar[name2]
                        ptx, pty = to_image_coords(ptx, _xlim), to_image_coords(pty, _ylim)
                        ax.scatter(ptx, pty, s=300, marker='X', **xmarkerkw)
                    if mpar:
                        ptx,  pty = mpar[name1], mpar[name2]
                        ptx, pty = to_image_coords(ptx, _xlim), to_image_coords(pty, _ylim)
                        ax.scatter(ptx, pty, s=150, marker='D', **dmarkerkw)
                    hide_ticks(ax)
                    for spine in ax.spines.values():
                        spine.set_visible(True)
                        spine.set_edgecolor(lightedge)
                if x < y:
                    r_accepted = accepted.iloc[::-1]
                    scat = ax.scatter(r_accepted[_names[y]],r_accepted[_names[x]], c=r_accepted["score"], cmap=cmap, s=80, alpha=0.5)
                    ax.set_xlim(_xlim)
                    ax.set_ylim(_ylim)
                    if simpar:
                        ptx,  pty = simpar[name1], simpar[name2]
                        ax.scatter(ptx, pty, s=200, marker='X', **xmarkerkw)
                    if mpar:
                        ptx,  pty = mpar[name1], mpar[name2]
                        # ax.scatter(ptx, pty, s=160, marker='D', c=marker_purple)
                        ax.scatter(ptx, pty, s=160, marker='D', **dmarkerkw)
                    hide_ticks(ax)
    #
    fig.subplots_adjust(left=0.1,
                    bottom=0.1, 
                    right=0.8, 
                    top=0.8, 
                    wspace=0.08, 
                    hspace=0.08)
    #
    _ax = axes[0,0]
    _ax.set_ylabel(_ax.get_ylabel(), labelpad=20)

    #
    cb_ax = fig.add_axes([0.83, 0.21, 0.02, 0.5])
    score = accepted["score"]
    cbspace = np.linspace(score.min(), score.max(), 3)
    # cbar.ax.set_yticks(cbspace)
    form = lambda x, _: "{:.2f}".format(x)
    cbar = fig.colorbar(scat, cax=cb_ax, ticks=cbspace, format=form)
    cbar.ax.invert_yaxis()
    cbar.ax.tick_params(labelsize=20) 

    # spine coloring doesn't work
    for spine in cb_ax.spines.values():
        spine.set_visible(True)
        spine.set_color(lightedge)

    return fig, axes


# 
def plot_accepted_stats(params, accepted, objectives, simref=None):
    #
    if 'lims' not in params.attrs:
        params.attrs['lims'] = {obs : (params[obs].min(), params[obs].max()) for obs in objectives}
    # untransformed accepted samples
    unt_accepted = params.iloc[accepted.index]
    n = len(objectives)
    fig, axes = plt.subplots(n,2,figsize = (2*5,n*5))

    for i, ax in enumerate(axes[:,0]):
        obs = objectives[i]
        obs_lims = params.attrs['lims'][obs]
        sns.histplot(unt_accepted[obs], binrange=obs_lims, ax=ax)
        if simref:
            _refv = simref[obs]
            ax.axvline(_refv, linestyle='--', c='k')
        ax.set_xlabel(obs)

    for i, ax in enumerate(axes[:,1]):
        obs = objectives[i]
        obs_lims = params.attrs['lims'][obs]
        sns.histplot(params[obs], binrange=obs_lims, ax=ax)
        ax.set_xlabel(obs)

# ----------------------------------------------------------------
# transform annoying anchor parameter between tmos and paper

def transform_anchor_parameter_data(psd):
    problem = psd["problem"]
    _bounds = problem["bounds"]
    _names = problem["names"]
    anchor_name = 'anchor_angle_smoothing_fraction'
    psd["params"][anchor_name] *= np.pi/2
    _i = _names.index(anchor_name)
    mn, mx = _bounds[_i]
    _bounds[_i] = (mn*np.pi/2, mx*np.pi/2)
    problem["bounds"] = _bounds

def transform_vret_parameter_data(psd):
    # kb_sh is defined as a rate of 0.004 micron retractions -> we want the velocity
    problem = psd["problem"]
    _bounds = problem["bounds"]
    _names = problem["names"]
    vret_name = 'kb_sh'
    contract = 0.004
    psd["params"][vret_name] *= contract
    _i = _names.index(vret_name)
    mn, mx = _bounds[_i]
    _bounds[_i] = (mn*contract, mx*contract)
    problem["bounds"] = _bounds


# ----------------------------------------------------------------
# dump relevant info

def describe_abc(psd, _accepted):
    print('--- describing accepted samples')
    N = len(_accepted)
    M = len(psd["lookup"][0])
    print('accepted {}/{}'.format(N,M))
    # print epsilon
    epsilon = _accepted['score'].iloc[-1]
    print('epsilon',  epsilon)
    # 
    m_par, v_par = mean_accepted(psd["problem"], _accepted)
    print('weighted means', m_par)
    # print('pilivar', m_par['pilivar'])
    # print('anchor * pi/2', m_par['anchor_angle_smoothing_fraction'] * np.pi/2)
    std_par = {name : np.sqrt(v_par[name]) for name in v_par}
    print('accepted standard deviation', std_par)
    # rescale this as well
    # print('accepted anchor std', np.pi/2 * std_par['anchor_angle_smoothing_fraction'])
    # print('---')

# 
def _statsuniformkde(data, x):
    _x = data[x].to_numpy()
    dens = sm.nonparametric.KDEMultivariate(data=(_x,), var_type='c',
        bw='cv_ls')
    return dens
    
Cinterval = collections.namedtuple('Cinterval', ["low", "high"])

from scipy.stats import bootstrap
# parameter inference
# cols [parameter, reference data, expected value, std, _min, _max, 95% confidence]
def tabulate_inference(problem, accepted, objectives, simpar=None, extra=False, use_mle=True):
    # print("using maximum likelihood estimate ", use_mle)
    _names = problem["names"]
    _bounds = problem["bounds"]
    col_simpar = np.full(len(_names), np.nan)
    if simpar:
        col_simpar = [simpar[name] for name in _names]
    mpar, vpar = mean_accepted(problem, accepted)

    # quantiles
    qmin, qmax = 0.05, 0.95
    col_qmin = [np.quantile(accepted[name], qmin) for name in _names]
    col_qmax = [np.quantile(accepted[name], qmax) for name in _names]
    epsilon = accepted['score'].iloc[-1]
    kernel = new_epanechnikov(epsilon)
    scores = accepted['score'].to_numpy()
    weights = kernel(scores)
    # col_ref = [reference_used for _ in range(len(_names))]

    # integrate Posterior projection
    confidence_data = {}
    mle = {}
    for i, name in enumerate(_names):
        low, high = _bounds[i]
        space = np.linspace(low, high, 1000)

        kernel = 'gau'
        kde1 = sm.nonparametric.KDEUnivariate(accepted[name])
        if kernel == 'gau':
            # automatic bandwidth selection
            dens = _statsuniformkde(accepted, name)
            # this object accepts weights
            kde1.fit(weights=weights, 
                    bw=dens.bw[0],
                    fft=False)


        pdf = kde1.evaluate(space)
        mle[name] = space[np.argmax(pdf)]
        cdf = np.cumsum(pdf) / np.sum(pdf)
        cinterval = Cinterval(space[np.searchsorted(cdf, qmin)], space[np.searchsorted(cdf, qmax)])
        confidence_data[name] = cinterval

    df = pd.DataFrame({
        "parameter": _names,
        # "reference data": col_ref,
        "reference value": col_simpar,
        "weighted mean": list(mpar.values()),
        # "weighted std": [np.sqrt(v) for v in vpar.values()],
        "MLE" : [mle[name] for name in _names],
        "confidence(0.05)": [confidence_data[name].low for name in _names],
        "confidence(0.95)": [confidence_data[name].high for name in _names]
        # "quantile(0.05)": col_qmin,
        # "quantile(0.95)": col_qmax,
    })


    if extra:
        # also compute a row for a table that compares the choice of statistics
        # convert the accepted parameters into the [0,1]
        unit = accepted.copy()
        error = {}
        interval = {}
        unit_error = {}
        for i, name in enumerate(_names):
            low, high = _bounds[i]
            def transform(x):
                return (x - low)/(high - low)
            par = unit[name].to_numpy()
            unit[name] = transform(par)
            if use_mle:
                mle = df['MLE'].iloc[i]
            else:
                mle = df['weighted mean'].iloc[i]

            unit_mle = transform(mle)
            if simpar:
                unit_sim = transform(simpar[name])
                unit_error[name] = unit_mle - unit_sim
                error[name] = (mle - simpar[name])/simpar[name]
            else:
                unit_error[name] = np.nan
                error[name] = np.nan


            lc, hc = transform(confidence_data[name].low), transform(confidence_data[name].high)
            interval[name] = hc - lc

        uerrv = np.array([unit_error[name] for name in _names])
        unit_error = np.linalg.norm(uerrv)
        S4 = np.prod(list(interval.values()))

        points = np.column_stack([unit[name] for name in _names])
        hull = scipy.spatial.ConvexHull(points)
        errv = list(error.values())

        n1, n2, n3, n4 = _names
        sd = pd.DataFrame({
            'objectives': ' '.join(objectives),
            '\epsilon' : epsilon,
            '\Delta '+n1 : error[n1],
            '\Delta '+n2 : error[n2],
            '\Delta '+n3 : error[n3],
            '\Delta '+n4 : error[n4],
            '\I'+n1 : interval[n1],
            '\I'+n2 : interval[n2],
            '\I'+n3 : interval[n3],
            '\I'+n4 : interval[n4],
            'mean error' : np.mean(np.abs(errv)),
            'max error' : np.max(np.abs(errv)),
            '|unit_error|' : unit_error,
            'max(unit_error)' : np.max(np.abs(uerrv)),
            'S^4': S4,
            'cxhvolume': hull.volume
            },
            index=[0])

        return df, sd

    else:
        return df


# ----------------------------------------------------------------
# testing

def test_kernel_function(delta):
    epanechnikov = new_epanechnikov(delta)
    rho_space = np.linspace(-2*delta,2*delta,1000)
    fig, ax = plt.subplots(figsize=(4,2))
    ax.plot(rho_space, [epanechnikov(x) for x in rho_space])
    ax.set_title("kernel for scores")
    plt.show()
