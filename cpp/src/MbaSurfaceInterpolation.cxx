/* MbaSurfaceInterpolation.cxx
 */

#ifdef _MSC_VER
#pragma warning ( disable : 4786 )
#endif

// C++ functions
#include <iostream>

// MBA libary
#include <MBA.h>
#include <UCButils.h>
#include <PointAccessUtils.h>

// entry point for the program
int main(int argc, char** argv)
{
  
  typedef std::vector<double> dVec;
  boost::shared_ptr<dVec> x_arr(new std::vector<double>);
  boost::shared_ptr<dVec> y_arr(new std::vector<double>);
  boost::shared_ptr<dVec> z_arr(new std::vector<double>);
  UCBspl::readScatteredData("/tmp/foo.csv", *x_arr, *y_arr, *z_arr);
  
  MBA mba(x_arr, y_arr, z_arr);
  
  // Create spline surface.
  mba.MBAalg(1,1.0240,7);
  
  UCBspl::SplineSurface surf = mba.getSplineSurface();
  
  std::cout << "umin = " << surf.umin() << std::endl;
  std::cout << "umax = " << surf.umax() << std::endl;

  // Find height and normal vector of surface in (x,y).
  double x = 1.52e-2, y = 1.17e-2;
  // double x = 9.64e-3, y = 6.40e-3;
  double z = surf.f(x,y);         
  // double nx, ny, nz;
  // surf.normalVector(x, y, nx, ny, nz);
  std::cout << "z-value in (x,y) = " << z << std::endl;
  // std::cout << "Normal in (x,y) = (" << nx << "," << ny << "," << nz << ")"
  //           << std::endl;
  
  // Sample surface and print to VRML file.
  UCBspl::printCSVgrid("/tmp/grid.csv", surf, 50, 50);
  UCBspl::printVRMLgrid("/tmp/vrmlgrid.wrl", surf, 50, 50);
  UCBspl::printGNUgrid("/tmp/gnugrid.txt", surf, 50, 50);
  
  return 0;
    // return EXIT_SUCCESS; 
}
