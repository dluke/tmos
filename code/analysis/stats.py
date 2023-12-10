#!/usr/bin/env python3

"""
Handle local save/load of json file for computed data
"""

# anyway to protect this file from being opened by multiple threads?
# USAGE: use stats.keep decorator and return a dictionary -- key names are not checked for uniqueness

import pickle
import os
import sys
import json

import numpy as np
import scipy.stats

from functools import wraps

DEBUG = False
localjson = 'local.json'


def load(json_path=localjson):
    if not os.path.exists(json_path):
        print('local.json does not exist at {} .\n Creating it...'.format(os.path.abspath(json_path)))
        dump({})
    with open(json_path, 'r') as f:
        if DEBUG:
            print('Reading from {}'.format(json_path))
        return json.load(f)


def dump(sd, json_path=localjson):
    with open(json_path, 'w') as f:
        if DEBUG:
            print('Writing to {}'.format(json_path))
        json.dump(sd, f, indent='\t')


def extend(sd, json_path=localjson):
    jd = load(json_path)
    jd.update(sd)
    dump(jd, json_path)


statsdir = 'stats/'


def namedkeep(json_name, statsdir=statsdir):
    if not json_name.endswith('.json'):
        json_name += '.json'
    if not os.path.exists(statsdir):
        os.mkdir(statsdir)
    json_path = os.path.join(statsdir, json_name)

    def keep(f):
        @wraps(f)
        def analysis_method(*args, **kw):
            sd = f(*args, **kw)
            if not isinstance(sd, dict):
                raise RuntimeError(
                    'Method {} expected to return dictionary.'.format(f.__name__))
            extend(sd, json_path)
            return sd
        return analysis_method
    return keep


def keep(f):
    @wraps(f)
    def analysis_method(*args, **kw):
        sd = f(*args, **kw)
        if not isinstance(sd, dict):
            raise RuntimeError(
                'Method {} expected to return dictionary.'.format(f.__name__))
        extend(sd)
        return sd
    return analysis_method

# TODO why doesn't this work?
# def keep(f):
#     return namedkeep(localjson, statsdir='.')

# import inspect
# print(keep.__name__)
# print(inspect.signature(keep))


def load_or_compute(named_var, compute):
    local = load()
    if named_var in local:
        return local[named_var]
    else:
        sd = compute()
        if not (isinstance(sd, dict)):
            raise RuntimeError('{} is not dict'.format(sd))
        # don't use compute() return value, instead load again to make sure we indeed saved the data
        local = load()
        if named_var not in local:
            raise RuntimeError('function {} expected to return quantity {}'.format(
                compute.__name__, named_var))
        return local[named_var]

# helper functions for calculating stats


def col_stats(trs, param_f, weight_f=None, each=True):
    # param_f a function that operates on a Track to obtain a parameter
    colsd = {}
    param = [param_f(tr) for tr in trs]
    if weight_f:
        weight =[weight_f(tr) for tr in trs]
        each_mean = [np.sum(par*w)/np.sum(w) for par, w in zip(param, weight)]
    else:
        each_mean = [np.mean(par) for par in param]
    col_all = np.concatenate(param)
    colsd['mean'] = np.mean(each_mean)
    colsd['min'] = float(col_all.min())
    colsd['max'] = float(col_all.max())
    colsd['median'] = np.median(col_all)
    colsd['std'] = np.std(col_all)
    colsd['std_error'] = scipy.stats.sem(col_all)
    #
    if each:
        colsd['each_mean'] = each_mean
        colsd['each_median'] = [np.median(par) for par in param]
    return colsd


def stats(arr):
    sdct = {}
    sdct['mean'] = np.mean(arr)
    sdct['median'] = np.median(arr)
    sdct['std'] = np.std(arr)
    sdct['std_error'] = scipy.stats.sem(arr)
    # each ...
    return sdct


# same idea but using pickle -- Can use this for more varied data (?)
# currently not used

# create meta data
metafile = 'local.meta'


def loadmeta():
    if os.path.exists(metafile):
        with open(metafile, 'r') as fp:
            meta = pickle.load(fp)
        return meta
    else:
        # create an empty meta data dictionary
        return {}


def dumpmeta(meta):
    with open(metafile, 'w') as fp:
        pickle.dump(meta, fp)


###################################################################################
# command line tool

def clean(var):
    jd = load()
    if var not in jd:
        raise ValueError(
            "Request to delete field {} does not exist.".format(var))
    print('Deleting field {}...'.format(var))
    del jd[var]
    dump(jd)


def print_keys():
    jd = load()
    for k in jd.keys():
        print(k)


def pretty_list(var):
    jd = load()
    try:
        var = jd[var]
    except KeyError:
        print("Requested field {} does not exist.".format(var))
        raise
    if isinstance(var, dict):
        return json.dumps(var, indent='\t')
    else:
        return str(var)


if __name__ == '__main__':
    import argparse

    parse = argparse.ArgumentParser()
    subparsers = parse.add_subparsers(dest='which')
    #
    parse_clean = subparsers.add_parser('clean')
    parse_clean.add_argument('clean')
    # TODO add documentation for clean command
    #
    parse_list = subparsers.add_parser('list')
    parse_list.add_argument('--keys', action='store_true')
    parse_list.add_argument('--var', action='store', default=None)
    args = parse.parse_args()

    if args.which == 'clean':
        clean(args.clean)
    elif args.which == 'list':
        if args.keys:
            print_keys()
        elif args.var:
            print(pretty_list(args.var))
        else:
            print(json.dumps(load(), indent='\t'))
