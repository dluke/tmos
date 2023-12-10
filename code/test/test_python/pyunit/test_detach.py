
import os
import sys
import unittest
import numpy as np


import tmos
from tmos.base import Vector3d
import tmos.base as base
import tmos.surface as surface
import tmos.pili 

import ctfp3d
import tfputils
import parameters
import pilush # for reconstructing frames

args = parameters.read()
R = 0.5
length = 2.0

# ---------------------------------------------------------------------------------
# setup test data 

import json
import txtdata
import pili
testcasepath = os.path.join(pili.root, "src/pyunit/detach_data.json")
if not os.path.exists(testcasepath):
    print('{} does not exist. Creating it.'.format(testcasepath))
    with open(testcasepath, 'w') as f:
        json.dump([],f,indent=1)

from collections import OrderedDict


# ---------------------------------------------------------------------------------
# construct cells

def get_attach(pr):
    return Vector3d(pr['attach_x'], pr['attach_y'], pr['attach_z'])
def get_anchor(pr):
    return Vector3d(pr['anchor_x'], pr['anchor_y'], pr['anchor_z'])
def get_paxis(pr):
    return Vector3d(pr['paxis_x'], pr['paxis_y'], pr['paxis_z'])

def construct(case):
    # IMPORTANT: setup up static variables using testcase config
    paramlist = parameters.ParameterList.parse_str(case['args'])
    paramlist.apply()
    # trl: one row from the tracking where the failed release occurs
    trl = case['track']
    evr = case['event']
    assert(evr['process'] == 'release')
    assert(evr['trigger'] == 'no_release')
    # create the pilus
    kpgen = ctfp3d.setup_kpgen(args)
    attachment, anchor, paxis = get_attach(evr), get_anchor(evr), get_paxis(evr)
    pilus = tmos.pili.PiliWLC3d(evr['pidx'], evr['plength'], evr['pleq'], 
        anchor, paxis, evr['ret_motor'], evr['ext_motor'], kpgen)
    pilus.isbound = evr['isbound']
    pilus.attachment = attachment
    # guessing
    if pilus.istaut():
        pilus.lastistaut = True
    pilus.cycles = 1
    pilus.new_cycle = False
    pilus.is_fully_retracted = False

    plane = surface.Plane()

    frame = pilush._construct_frame(trl)
    caps = tmos.surface.Capsule(frame, 0.5, 2.0)
    cell = tmos.pili.CellWLC3d(0, 0, caps, plane, kpgen)
    cell.add_pilus(pilus)
    
    assert(evr['nseg'] == 1) # we dont reconstruct the chain details because we don't know it

    return cell

# ---------------------------------------------------------------------------------
# construct test data

def row_to_dict(row):
    # convert row of numpy structured array to dictionary
    return OrderedDict([( name, txtdata.txtform[name][0](row[name]) ) for name in row.dtype.names])

def create_test_data(trl, evr, args, caseid=None):
    with open(testcasepath, 'r') as f:
        casedata = json.load(f)
    case = {'track':row_to_dict(trl), 'event':row_to_dict(evr)}
    case['args'] = str(args)
    # 
    if caseid is None or len(casedata) == caseid:
        caseid = len(casedata)
        casedata.append(case)
    else:
        assert (caseid < len(casedata))
        casedata[caseid] = case
    with open(testcasepath, 'w') as f:
        print('writing case {} to {}'.format(caseid, testcasepath))
        f.write('[\n')
        for case in casedata:
            f.write(json.dumps(case))
            f.write('\n')
        f.write(']')

def load_data():
    with open(testcasepath, 'r') as f:
        casedata = json.load(f)
    return casedata


# ---------------------------------------------------------------------------------
# 

class TestDetach(unittest.TestCase):

    def setUp(self):
        kpgen = ctfp3d.setup_kpgen(args)
        pl = surface.Plane()
        bodyctr = Vector3d(0,0,0.495)
        bodyaxis = Vector3d(1,0,0)
        e1 = Vector3d(0,1,0) # rotation axis
        e2 = bodyaxis.cross(e1)
        frame = base.Frame(bodyctr, e1, e2, bodyaxis)
        cap = surface.Capsule(frame, 0.5, 2.)
        cell = tmos.pili.CellWLC3d(0, 0, cap, pl, kpgen)
        cell.pilivar = 4.0
        cell.common_init()
        self.cell = cell
        self.casedata = load_data()
        
    def test_spawn_pilus(self):
        pilus = self.cell.spawn_pilus()

    def test_case_00(self):
        case = self.casedata[0]
        cell = construct(case)
        self.assertEqual(tmos.pili.Pili.allow_bound_ext_motor, False)
        self.assertEqual(tmos.pili.Pili.allow_bound_extension, True)
        self.assertEqual(tmos.pili.Pili.force_bound_retraction, True)
        # -------------------------------------------------------------------------
        pilus = cell.pili[0]
        print('pilus attached at', pilus.attachment)
        chain = pilus.transformed_chain_instance(cell.body)
        bound_pilus = pilus.attachment - chain.source
        print('bound_pilus', bound_pilus, bound_pilus.len())
        print('bound_pilus axis', bound_pilus * (1/bound_pilus.len()))
        lab_axis = cell.body.frame.to_lab_rt(chain.axis)
        print('pilus lab axis', lab_axis)

        print('chain first segment', chain.lines[0], 'length', chain.lines[0].len())
        print('target ', chain.targets[0])
        print('axis', chain.axis)
        print('source', chain.source)
        self.assertEqual(pilus.detach_grace_length, 0.016)
        print('lasta', pilus.lasta)
        pre_detach_leq = pilus.leq

        self.assertEqual(pilus.isbound, 1)
        pilus.detach(0.)
        print()
        
        # -------------------------------------------------------------------------
        # after detach
        chain = pilus.transformed_chain_instance(cell.body)
        print('chain first segment', chain.lines[0], 'length', chain.lines[0].len())
        print('target ', chain.targets[0])
        print('axis', chain.axis)
        print('source', chain.source)
        print(pilus)
        self.assertEqual(pilus.isbound, 0)
        self.assertEqual(pilus.leq, pre_detach_leq - pilus.detach_grace_length)

        # -------------------------------------------------------------------------
        # attempt attach
        _err_at = pilus.attach(0., cell.body, cell.surface)
        # self.assertEqual(pilus.isbound, 0)

    def test_shorten_to_detach(self):
        case = self.casedata[0]
        cell = construct(case)
        pilus = cell.pili[0]
        shortening = pilus.shorten_to_detach(cell.body, cell.surface)
        self.assertTrue(shortening > 0.)
        pilus.shrink_by(shortening + 0.004)
        ztarget = pilus.transformed_chain_instance(cell.body).targets[-1].z
        self.assertTrue( ztarget > 0 )

    def test_detach(self):
        case = self.casedata[0]
        cell = construct(case)
        tmos.pili.Pili.force_detachment = True
        pilus = cell.pili[0]
        self.assertTrue(pilus.last_shortening == 0)
        cell.detach(pilus, 0.)
        self.assertTrue(pilus.last_shortening > 0)


if __name__ == '__main__':
    unittest.main()
