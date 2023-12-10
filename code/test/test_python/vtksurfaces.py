import os

import tmos.base as base
import tmos.surface as surface 
import tmos.vtkwriter as vtkwriter

surfacedir = './surfaces/'


def write_hgrid(): 
    print("writing hexsurface representation to file.")
    R = 0.5
    hgrid = surface.HexGrid(R)
    size = 20
    res = 20
    vtkwriter.write_hgrid(hgrid, size, res)

    original_name = "HexGrid.vtk"
    outname = "Hgrid_R_{:06.3f}.vtk".format(R)
    os.rename(original_name, outname)


def write_infsteps(): 
    print("writing infsteps representation to file.")
    height = 5.
    smallr = 0.25
    sep = 10.
    infsteps = surface.InfSteps(height, sep, smallr)
    vtkwriter.write_infsteps(infsteps)

def write_segplane():
    print("construct segplane")
    splane = surface.SegPlane()
    print("writing segplane representation to file.")
    original_name = "SegPlane.vtk"
    vtkwriter.write_segplane(splane)
    outname = os.path.join(surfacedir, original_name)
    os.rename(original_name, outname)

def write_sineplane():
    print("construct sineplane")
    frame = base.Frame()
    splane = surface.SinePlane(frame, 1., 1., 1)
    print("writing sineplane representation to file.")
    #print splane.A, splane.B, splane.invB, splane.xperiod
    #print splane.form(1)
    original_name = "SinePlane.vtk"
    vtkwriter.write_sineplane(splane)
    outname = os.path.join(surfacedir, original_name)
    os.rename(original_name, outname)

#write_segplane()
#write_sineplane()
write_hgrid()

