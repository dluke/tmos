

import tmos.base as base
from tmos.base import Vector3d


def test_frame():
    origin = Vector3d(0,0,0)
    ax = Vector3d(1,0,0)
    rax = Vector3d(0,1,0)
    frame = base.Frame(origin, ax, rax)
    print(frame)
    frame = base.Frame(origin, ax, rax, Vector3d(0,0,1))
    print(frame)
test_frame()


