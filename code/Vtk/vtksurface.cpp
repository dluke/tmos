
#include "vtksurface.hpp"
#include <iostream>

using std::vector;


int write_infsteps(InfSteps& step) {
  // const , make variable
  int Nsteps = 5; // +/minus 5 steps in both directions
  int npts = 21;
  double ywidth = 50.;

  double sep = step.get_sep();
  double xlim = 5 * sep;
  double smallr = step.get_smallr();
  double chunk = sep/2. - smallr;

  double delta = (2*smallr)/(npts-1);
  vector<double> xspace;


  double x;
  for (int i = -Nsteps; i < Nsteps; i++) {
    xspace.push_back(i * sep);
    double xc = i * sep + chunk;
    for (int j = 0; j < npts; j++) {
      x = xc + j * delta;
      xspace.push_back(x);
    }
  }

  // last element
  xspace.push_back(xlim);

  // use space to a VTK polyline
  vtkSmartPointer<vtkPoints> pts =
    vtkSmartPointer<vtkPoints>::New();

  vtkSmartPointer<vtkPolyLine> polyLine = 
    vtkSmartPointer<vtkPolyLine>::New();
  polyLine->GetPointIds()->SetNumberOfIds(npts*2*Nsteps + 2*Nsteps + 1); // needed?

  int ptid = 0;
  for (double x : xspace) {
    double pt[3] = {x, -ywidth, step(x)};
    pts->InsertNextPoint(pt);
    polyLine->GetPointIds()->SetId(ptid, ptid);
    ptid++;
  }

  vtkSmartPointer<vtkCellArray> cells = 
    vtkSmartPointer<vtkCellArray>::New();
  cells->InsertNextCell(polyLine);
 
  // Create a polydata to store everything in
  vtkSmartPointer<vtkPolyData> polyData = 
    vtkSmartPointer<vtkPolyData>::New();
 
  // Add the points to the dataset
  polyData->SetPoints(pts);
 
  // Add the lines to the dataset
  polyData->SetLines(cells);

  // Linear Extrusion
  vtkSmartPointer<vtkLinearExtrusionFilter> extrude = 
    vtkSmartPointer<vtkLinearExtrusionFilter>::New();
  extrude->SetInputData(polyData);
  extrude->SetExtrusionTypeToNormalExtrusion();
  extrude->SetVector(0, 1, 0 );
  extrude->SetScaleFactor(2*ywidth);
  extrude->Update();

  // final surface object
  vtkSmartPointer<vtkPolyData> surface = 
    vtkSmartPointer<vtkPolyData>::New();
 
  surface->ShallowCopy(extrude->GetOutput());
  
  vtkSmartPointer<vtkXMLPolyDataWriter> writer = 
    vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  cout << "Writing infsteps to " << "InfSteps.vtk" << endl;
  writer->SetFileName("InfSteps.vtk");
  writer->SetInputData(surface);
  writer->SetDataModeToAscii();
  writer->Write();

  return EXIT_SUCCESS;
}


/*int write_step(Step& step)*/
//{
  //cout << "preparing to write surface step" << endl;
  //double delta = 0.1;

  //vtkSmartPointer<vtkPoints> pts =
    //vtkSmartPointer<vtkPoints>::New();

  //vtkSmartPointer<vtkPolyLine> polyLine = 
    //vtkSmartPointer<vtkPolyLine>::New();
  //polyLine->GetPointIds()->SetNumberOfIds(101); // 10/delta + 1

  //int ptid = 0;
  //for (double x = -5; x < 5+delta/2.; x += delta)
  //{
    //double pt[3] = {x, -3, step(x)};
    //pts->InsertNextPoint(pt);
    //polyLine->GetPointIds()->SetId(ptid, ptid);
    //ptid++;
  //}

  //vtkSmartPointer<vtkCellArray> cells = 
    //vtkSmartPointer<vtkCellArray>::New();
  //cells->InsertNextCell(polyLine);
 
  //// Create a polydata to store everything in
  //vtkSmartPointer<vtkPolyData> polyData = 
    //vtkSmartPointer<vtkPolyData>::New();
 
  //// Add the points to the dataset
  //polyData->SetPoints(pts);
 
  //// Add the lines to the dataset
  //polyData->SetLines(cells);

  //// Linear Extrusion
  //vtkSmartPointer<vtkLinearExtrusionFilter> extrude = 
    //vtkSmartPointer<vtkLinearExtrusionFilter>::New();
  //extrude->SetInputData(polyData);
  //extrude->SetExtrusionTypeToNormalExtrusion();
  //extrude->SetVector(0, 1, 0 );
  //extrude->SetScaleFactor(6);
  //extrude->Update();

  //// final surface object
  //vtkSmartPointer<vtkPolyData> surface = 
    //vtkSmartPointer<vtkPolyData>::New();
 
  //surface->ShallowCopy(extrude->GetOutput());
  
  //vtkSmartPointer<vtkXMLPolyDataWriter> writer = 
    //vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  //cout << "Writing step to " << "Step.vtk" << endl;
  //writer->SetFileName("Step.vtk");
  //writer->SetInputData(surface);
  //writer->SetDataModeToAscii();
  //writer->Write();

  //return EXIT_SUCCESS;
//}


int write_hgrid(HexGrid& Hgrid, int range, int res)
{

  // going to append spheres to this append filter
  vtkSmartPointer<vtkAppendPolyData> appendFilter =
    vtkSmartPointer<vtkAppendPolyData>::New();

  // get a list of hexes
  vector<Hexc> hexes = Hgrid.coordinate_range(HexGrid::origin, range);

  // iterate through the hexes and append a sphere size Hgrid.R at that point
  for (Hexc hx : hexes) {
    Vector3d target = Hgrid.get_xyz(hx);
    // create a polydata sphere source
    vtkSmartPointer<vtkPolyData> spherepd
        = vtkSmartPointer<vtkPolyData>::New();
    // create a sphere source
    vtkSmartPointer<vtkSphereSource> sphere = 
        vtkSmartPointer<vtkSphereSource>::New();
    // set sphere specific values
    sphere->SetCenter(target.x, target.y, target.z);
    sphere->SetRadius(Hgrid.R);
    sphere->SetPhiResolution(res);
    sphere->SetThetaResolution(res);
    sphere->Update();
// at this step we create a mesh object so the spehre representation becomes data intensive
// can display spheres more efficiently?
    spherepd->ShallowCopy(sphere->GetOutput());

    appendFilter->AddInputData(spherepd);
    // need to call every loop?
    appendFilter->Update();
  }

  // 
  vtkSmartPointer<vtkXMLPolyDataWriter> writer = 
    vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName("HexGrid.vtk");
  writer->SetInputConnection(appendFilter->GetOutputPort());
  writer->SetDataModeToAscii();
  writer->Write();

}

vector<Vector3d> get_polyline(SegPlane& sp)
{
  // return the vertices of the Segmented plane 
  vector<Vector3d> vline;
  for (int i = -sp.xnwidth; i < sp.xnwidth; ++i) {
    // iterating through repeating elements 
    vector<Vector3d> verts = sp.get_element_at(i)->verts;
    vline.insert(vline.end(), verts.begin(), verts.end());
  }
  // finally add the last point
  vline.push_back(sp.get_element_at(sp.xnwidth-1)->get_lmpt());
  return vline;
}


vtkSmartPointer<vtkPolyData> extrude_surface(vector<Vector3d> vline, double yextent)
{

  // use space to a VTK polyline
  vtkSmartPointer<vtkPoints> pts =
    vtkSmartPointer<vtkPoints>::New();

  vtkSmartPointer<vtkPolyLine> polyLine = 
    vtkSmartPointer<vtkPolyLine>::New();
  polyLine->GetPointIds()->SetNumberOfIds(vline.size()); 

  int ptid = 0;
  for (Vector3d v : vline) {
    double pt[3] = {v.x, -yextent, v.z};
    pts->InsertNextPoint(pt);
    polyLine->GetPointIds()->SetId(ptid, ptid);
    ptid++;
  }

  vtkSmartPointer<vtkCellArray> cells = 
    vtkSmartPointer<vtkCellArray>::New();
  cells->InsertNextCell(polyLine);
 
  // Create a polydata to store everything in
  vtkSmartPointer<vtkPolyData> polyData = 
    vtkSmartPointer<vtkPolyData>::New();
 
  // Add the points to the dataset
  polyData->SetPoints(pts);
 
  // Add the lines to the dataset
  polyData->SetLines(cells);

  // Linear Extrusion
  vtkSmartPointer<vtkLinearExtrusionFilter> extrude = 
    vtkSmartPointer<vtkLinearExtrusionFilter>::New();
  extrude->SetInputData(polyData);
  extrude->SetExtrusionTypeToNormalExtrusion();
  extrude->SetVector(0, 1, 0 );
  extrude->SetScaleFactor(2*yextent);
  extrude->Update();

  // final surface object
  vtkSmartPointer<vtkPolyData> surface = 
    vtkSmartPointer<vtkPolyData>::New();
 
  surface->ShallowCopy(extrude->GetOutput());
  return surface;
}


int write_segplane(SegPlane& sp)
{
  double yextent = 50; // extent in infinite direction 

  // get the polyline has a vector of Vector3d

  vector<Vector3d> vline = get_polyline(sp);
  vtkSmartPointer<vtkPolyData> surface = extrude_surface(vline, yextent);

  // write this 
  vtkSmartPointer<vtkXMLPolyDataWriter> writer = 
    vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  //cout << "Writing infsteps to " << "InfSteps.vtk" << endl;
  writer->SetFileName("SegPlane.vtk");
  writer->SetInputData(surface);
  writer->SetDataModeToAscii();
  writer->Write();

  return EXIT_SUCCESS;
}

vector<Vector3d> get_polyline(SinePlane& sp, int nperiod)
{
  double period, step, x;
  int density = 100; // 100 points per period
  period = sp.B * 2 * M_PI;
  step = period/density;
  vector<Vector3d> vline;
  int nmax = density * 2 * nperiod + 1;
  vline.resize(nmax-1);
  x = -nperiod*period;
  int i;
  for (i = 0; i < nmax-1; i++) {
    vline[i] = sp.sp_form(x);
    x += step;
    //cout << x << " " << sp.form(x) << endl;
  }
  //vline[i] = sp.sp_form(nperiod*period);
  return vline;
}


int write_sineplane(SinePlane& sp)
{
  double yextent = 50; // extent in infinite direction 
  int nperiod = 5;

  // get the polyline has a vector of Vector3d

  vector<Vector3d> vline = get_polyline(sp, nperiod);
  vtkSmartPointer<vtkPolyData> surface = extrude_surface(vline, yextent);

  // write this 
  vtkSmartPointer<vtkXMLPolyDataWriter> writer = 
    vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  //cout << "Writing infsteps to " << "InfSteps.vtk" << endl;
  writer->SetFileName("SinePlane.vtk");
  writer->SetInputData(surface);
  writer->SetDataModeToAscii();
  writer->Write();

  return EXIT_SUCCESS;
}
