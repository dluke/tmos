
#include <array>
#include <vector>
#include <string>
#include <exception>

#include <iostream>

#include "vtkwriter.hpp"
#include "vector3d.hpp"
#include "matrix3d.hpp"

using std::cout;
using std::endl;

using std::vector;
using std::string;

// Haven't figured out how to set colors properly in paraview 5
// Can rework set single value scalar colors using a broad color spectrum

// Create two colors - one for each line
//const unsigned char white[3] = { 255, 255, 255 };
const unsigned char white[3] = { 255, 255, 255 };

const unsigned char red[3] = { 200, 0, 0 };
const unsigned char yellow[3] = { 170, 230, 0};
//const unsigned char redflash[3] = { 170, 80, 30 };
const unsigned char blue[3] = { 0, 0, 255 };
const unsigned char blueflash[3] = { 50, 100, 255 };

/*const unsigned char red[3] = { 255, 0, 0 };*/
//const unsigned char redflash[3] = { 255, 80, 30 };
//const unsigned char blue[3] = { 0, 0, 255 };
//const unsigned char blueflash[3] = { 100, 100, 255 };

//const unsigned char grey[3] = { 102, 102, 102};
/*const unsigned char purple[3] = { 204, 0, 204 };*/

const double blue_scale = 50;
const double red_scale = 255;
const double yellow_scale = 200;


typedef std::array<double, 3> array;

std::array<double, 3> toarray(Vector3d v) {
  std::array<double,  3> arr = {v.x, v.y, v.z};
  return arr;
}

void InsertNextTupleValue(vtkSmartPointer<vtkUnsignedCharArray> car, const unsigned char col[]) {
  car->InsertNextTuple3(col[0], col[1], col[2]);
}

int add_pilus_line(
    vtkSmartPointer<vtkPoints> pts, 
    vtkSmartPointer<vtkCellArray> lines,
    vtkSmartPointer<vtkDoubleArray> colors,
    Cell3d& cell,
    PiliRod3d& pilus,
    int ptid)
{
  // add anchor point
  Vector3d anc = cell.get_lab_anchor(pilus);
  double pt[3] = {anc.x, anc.y, anc.z};
  // add target point
  Vector3d tgt;
  if (pilus.isbound) {
    tgt = pilus.attachment;
  }
  else {
    tgt = anc + pilus.leq * cell.get_lab_axisEq(pilus);
  }
  double topt[3] = {tgt.x, tgt.y, tgt.z};

  pts->InsertNextPoint(pt);
  pts->InsertNextPoint(topt);
  vtkSmartPointer<vtkLine> line =
    vtkSmartPointer<vtkLine>::New();
  line->GetPointIds()->SetId(0, ptid); // anchor
  line->GetPointIds()->SetId(1, ptid+1); // pili point
  ptid += 2;
  lines->InsertNextCell(line);

  // colours
  if (pilus.isbound) {
    if (pilus.istaut()) {
      colors->InsertNextValue(red_scale);
    }
    else {
      colors->InsertNextValue(yellow_scale);
    }
  }
  else {
    colors->InsertNextValue(blue_scale);
  }

  return ptid;
}

int add_pilus_chain(
    vtkSmartPointer<vtkPoints> pts, 
    vtkSmartPointer<vtkCellArray> lines,
    vtkSmartPointer<vtkDoubleArray> colors,
    Cell3d& cell,
    PiliWLC3d& pilus,
    int ptid)
{

  Vector3d anc = cell.get_lab_anchor(pilus);
  if (pilus.isbound) {
    // add anchor point
    double source[3] = {anc.x, anc.y, anc.z};
    pts->InsertNextPoint(source);

    //
    //Vector3d tgt = anc + pilus.leq * cell.get_lab_axisEq(pilus);
    Vector3d tgt = pilus.get_attachment();
    double topt[3] = {tgt.x, tgt.y, tgt.z};
    pts->InsertNextPoint(topt);
    vtkSmartPointer<vtkLine> line =
      vtkSmartPointer<vtkLine>::New();
    line->GetPointIds()->SetId(0, ptid); // anchor
    line->GetPointIds()->SetId(1, ptid+1); // pili point
    lines->InsertNextCell(line);
    ptid += 2;

    if (pilus.istaut()) {
      colors->InsertNextValue(red_scale);
    }
    else {
      colors->InsertNextValue(yellow_scale);
    }
    
  }
  else {
    // if the pili is not bound then we add all the line segments
    std::shared_ptr<Chain> chain = pilus.new_chain_instance();
    Chain& pchain = *chain;
    pchain.set_lasta(pilus.lasta);
    // Chain pchain{chain}; // copy

    if (true)
    {
      pchain.rotate(cell.body.frame.get_rmatrix());
      pchain.translate(anc);
      pchain.compute_targets();

      Vector3d run = pchain.source;
      double r[3] = {run.x, run.y, run.z};
      pts->InsertNextPoint(r);
      ptid += 1;

      // get line segments
      for (Vector3d target : pchain.targets)
      {
        double s[3] = {target.x, target.y, target.z};
        pts->InsertNextPoint(s);
        vtkSmartPointer<vtkLine> line =
            vtkSmartPointer<vtkLine>::New();
        line->GetPointIds()->SetId(0, ptid-1); // source
        line->GetPointIds()->SetId(1, ptid); // target
        lines->InsertNextCell(line);
        ptid += 1;

        colors->InsertNextValue(blue_scale);
      }
    }
  }
    

  return ptid;
}

    
// just the attached pili as reconstructed from the text output
// Need to reconstruct the rotational state of the body from the e1 and axis vectors
int write_pili_vectors(Capsule& body, lsegl ppts, char* fout)
{
}

int write_pili3d(Cell3d& cell, char* fout)
{
  vtkSmartPointer<vtkPolyData> allpili =
    vtkSmartPointer<vtkPolyData>::New();
 
  // storage for points
  vtkSmartPointer<vtkPoints> pts =
    vtkSmartPointer<vtkPoints>::New();
  // storage for lines
  vtkSmartPointer<vtkCellArray> lines =
    vtkSmartPointer<vtkCellArray>::New();

  // Create a vtkUnsignedCharArray container and store the colors in it
  vtkSmartPointer<vtkDoubleArray> colors =
    vtkSmartPointer<vtkDoubleArray>::New();
  colors->SetNumberOfComponents(1);
  colors->SetName("Color");
 
  // calculate bound pili positions
  //
  int ptid = 0;
  for (int i = 0 ; i < cell.pili.size() ; i++) {
    Pili* pilus = cell.pili[i].get();

    if (pilus->get_typestr() == "rod")
    {
      PiliRod3d *pp = dynamic_cast<PiliRod3d *>(pilus);
      ptid = add_pilus_line(pts, lines, colors, cell, *pp, ptid);
    }
    else if (pilus->get_typestr() == "wlc")
    {
      PiliWLC3d *pp = dynamic_cast<PiliWLC3d *>(pilus);
      ptid = add_pilus_chain(pts, lines, colors, cell, *pp, ptid);
    }
    else { throw std::runtime_error("Failed to identify pilus object type."); }
  }

  // Add the points to the polydata container
  allpili->SetPoints(pts);
  // Set the points and vertices we created as the geometry and topology of the polydata
  allpili->SetLines(lines);

  allpili->GetCellData()->SetScalars(colors);

  // writeout
  vtkSmartPointer<vtkXMLPolyDataWriter> writer 
    = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName(fout);
  writer->SetInputData(allpili);
 
  // Optional - set the mode. The default is binary.
  //writer->SetDataModeToBinary();
  writer->SetDataModeToAscii();
 
  writer->Write();

  return EXIT_SUCCESS;
}

// the constant that controls the number of facets in the polydata sphere and cylinder
const int res = 20;

// can get this data from text file
int write_cell3d_from_body(Capsule& body, char* fout)
{
  Vector3d oopt = body.get_headpt();
  Vector3d tailpt = body.get_endpt();
  double R = body.R;
  double length = body.length;

  Vector3d axis = (oopt - tailpt).unit();
  Vector3d cpt = 0.5 * (oopt + tailpt);
  //cout << "writing to vtk file" << endl;
  //cout << "axis " << axis << endl;
  //cout << "headpt " << oopt << endl;
  //cout << "tailpt " << tailpt << endl;

  // PolyData 
  vtkSmartPointer<vtkPolyData> spherepd
    = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkPolyData> cylinderpd
    = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkPolyData> sphere2pd
    = vtkSmartPointer<vtkPolyData>::New();

  // Create a sphere
  vtkSmartPointer<vtkSphereSource> sphere = 
    vtkSmartPointer<vtkSphereSource>::New();
  sphere->SetCenter(oopt.x, oopt.y, oopt.z);
  sphere->SetRadius(R);
  // Make the surface smooth.
  sphere->SetPhiResolution(res);
  sphere->SetThetaResolution(res);
  sphere->Update();
  spherepd->ShallowCopy(sphere->GetOutput());

  // Create a cylinder 
  vtkSmartPointer<vtkCylinderSource> cylinder = 
    vtkSmartPointer<vtkCylinderSource>::New();
  cylinder->SetCenter(cpt.x, cpt.y, cpt.z);
  cylinder->SetRadius(R);
  cylinder->SetHeight(length);
  // Make the surface smooth.
  cylinder->SetResolution(res);
  cylinder->Update();
  cylinderpd->ShallowCopy(cylinder->GetOutput());

  // Cylinder default is with rotation axis in e_y direction
  double theta = (180/M_PI) * std::acos(axis.dot(e_y));
  Vector3d omega = -axis.cross(e_y);

  vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
  transform->PostMultiply();
  transform->Translate(-cpt.x, -cpt.y, -cpt.z);
  transform->RotateWXYZ(theta, omega.x, omega.y, omega.z);
  transform->Translate(cpt.x, cpt.y, cpt.z);
  
  vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter = 
      vtkSmartPointer<vtkTransformPolyDataFilter>::New();
  transformFilter->SetTransform(transform);
  transformFilter->SetInputData(cylinderpd);
  transformFilter->Update();
  //transformFilter->GetOutput(); 
  // do I need to do a shallowcopy operation like last time?
  cylinderpd->ShallowCopy(transformFilter->GetOutput());

  // Create the other sphere
  vtkSmartPointer<vtkSphereSource> sphere2 = 
    vtkSmartPointer<vtkSphereSource>::New();
  sphere2->SetCenter(tailpt.x, tailpt.y, tailpt.z);
  sphere2->SetRadius(R);
  // Make the surface smooth.
  sphere2->SetPhiResolution(100);
  sphere2->SetThetaResolution(100);
  sphere2->Update();
  //vtkSmartPointer<vtkPolyData> sphere2pd = sphere2->GetOutput();
  sphere2pd->ShallowCopy(sphere2->GetOutput());

  vtkSmartPointer<vtkAppendPolyData> appendFilter =
    vtkSmartPointer<vtkAppendPolyData>::New();
  appendFilter->AddInputData(cylinderpd);
  appendFilter->AddInputData(spherepd);
  appendFilter->AddInputData(sphere2pd);
  appendFilter->Update();

  vtkSmartPointer<vtkXMLPolyDataWriter> writer = 
    vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName(fout);
  writer->SetInputConnection(appendFilter->GetOutputPort());
  writer->SetDataModeToAscii();
  writer->Write();

  return EXIT_SUCCESS;
}

// used by simulation code
int write_cell3d(Cell3d& cell, char* fout)
{
  return write_cell3d_from_body(cell.body, fout);
}

