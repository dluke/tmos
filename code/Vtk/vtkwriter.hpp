
#ifndef __VTK_HPP__
#define __VTK_HPP__

#include <cmath>
#include <memory>
#include <vector>
#include <utility>

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkUnsignedCharArray.h>
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkLine.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkPolyData.h>
#include <vtkFieldData.h>
#include <vtkPolyDataMapper.h>
//#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkXMLPolyDataWriter.h>
    
#include <vtkCleanPolyData.h>
#include <vtkAppendPolyData.h>
#include <vtkSphereSource.h>
#include <vtkCylinderSource.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkTransform.h>

#include "lvec.hpp"
#include "chain.hpp"
#include "pili.hpp"
#include "cell.hpp"
#include "cell3d.hpp"

using std::vector;

//can pybind11 implicitely convert a list of tuples of Vector3d to this type?
typedef std::vector<std::pair<Vector3d, Vector3d>> lsegl;

int write_pili3d(Cell3d& cell, char* fout);
int write_cell3d(Cell3d& cell, char* fout);

int write_cell3d_from_body(Capsule& body, char* fout);
int write_pili_vectors(Capsule& body, lsegl ppts, char* fout);


#endif

