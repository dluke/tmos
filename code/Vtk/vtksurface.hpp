
#ifndef __VTKSURFACE_HPP__
#define __VTKSURFACE_HPP__

#include <memory>
#include <vector>

#include <vtkSmartPointer.h>
#include <vtkLinearExtrusionFilter.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPolyLine.h>
#include <vtkPolyData.h>
#include <vtkXMLPolyDataWriter.h>

#include <vtkAppendPolyData.h>
#include <vtkSphereSource.h>


//
#include "vector3d.hpp"
#include "cylinder.hpp"
#include "hexsurface.hpp"
#include "curvestep.hpp"
#include "periodic.hpp"

int write_infsteps(InfSteps& step);

int write_hgrid(HexGrid& Hgrid, int range = 5, int res = 20);

int write_segplane(SegPlane& sp);

int write_sineplane(SinePlane& sp);

#endif
