#include <stdio.h>
#include <math.h>
#define M_PI        3.14159265358979323846264338327950288  

struct GeometryParams {
  double angle;
  double radius; 
};

/** 00. Default function, sdf geometry and gradient for a horizontal fault*/
void horizontal_sdf(double coor[], struct GeometryParams GeoParamList, double *phi)
{
  *phi = coor[1];
}
void horizontal_grad_sdf(double coor[], struct GeometryParams GeoParamList, double grad[])
{
  grad[0] = 0;
  grad[1] = 1.0;
}

/** 01. Counterclock-wise Tilted Function: sdf geometry and gradient*/
void tilted_sdf(double coor[], struct GeometryParams GeoParamList, double *phi)
{
    *phi = -sin(GeoParamList.angle* M_PI/180.0) * coor[0]+ cos(GeoParamList.angle* M_PI/180.0) * coor[1];
}

void tilted_grad_sdf(double coor[], struct GeometryParams GeoParamList, double grad[])
{
  grad[0] = -sin(GeoParamList.angle* M_PI/180.0);
  grad[1] = cos(GeoParamList.angle* M_PI/180.0);
}


/** Definition of function to pointer for SDF */
void (*sdf_func[])(double coor[], struct GeometryParams GeoParamList, double *phi) =
  {horizontal_sdf, tilted_sdf};
void (*sdf_grad_func[])(double coor[], struct GeometryParams GeoParamList, double grad[]) =
  {horizontal_grad_sdf, tilted_grad_sdf};



/** ----------------
 *  Functions for Unit testing
 *  -----------------*/
void PrintCoordSDFvalues(double coor[], double grad[], double phi)
{
  printf("-- Coords: %.2f, %.2f ", coor[0], coor[1]);
  printf("-- Grad: %.2f, %.2f ", grad[0], grad[1]);
  printf("-- SDF Phi: %.2f\n\n", phi);
}

void UnitTest_SDFgeometry()
{
  double FakeCoord[2];
  double FakeGrad[2];
  double FakePhi;
  struct GeometryParams FakeGeoParam = {0}; //Initialized to null

  
  FakeCoord[0] = -2.0;
  FakeCoord[1] =  2.0;

  printf("Test 1 - Horizontal geometry\n");
  (*sdf_func[0])(FakeCoord, FakeGeoParam,  &FakePhi);
  (*sdf_grad_func[0])(FakeCoord, FakeGeoParam,  FakeGrad);
  PrintCoordSDFvalues(FakeCoord, FakeGrad, FakePhi);

  printf("Test 2 - Tilted (00°) geometry\n");
  FakeGeoParam.angle = 00;
  (*sdf_func[1])(FakeCoord, FakeGeoParam,  &FakePhi);
  (*sdf_grad_func[1])(FakeCoord, FakeGeoParam,  FakeGrad);
  PrintCoordSDFvalues(FakeCoord, FakeGrad, FakePhi);

  printf("Test 3 - Tilted (90°) geometry\n");
  FakeGeoParam.angle = 90 ;
  (*sdf_func[1])(FakeCoord, FakeGeoParam,  &FakePhi);
  (*sdf_grad_func[1])(FakeCoord, FakeGeoParam,  FakeGrad);
  PrintCoordSDFvalues(FakeCoord, FakeGrad, FakePhi);

  printf("Test 4 - Tilted (45°) geometry\n");
  FakeGeoParam.angle = 45;
  (*sdf_func[1])(FakeCoord, FakeGeoParam,  &FakePhi);
  (*sdf_grad_func[1])(FakeCoord, FakeGeoParam,  FakeGrad);
  PrintCoordSDFvalues(FakeCoord, FakeGrad, FakePhi);

  FakeCoord[0] =  2.0;
  FakeCoord[1] =  2.0;

  printf("Test 3 - Tilted (90°) geometry and negative x component of the coord\n");
  FakeGeoParam.angle = 90 ;
  (*sdf_func[1])(FakeCoord, FakeGeoParam,  &FakePhi);
  (*sdf_grad_func[1])(FakeCoord, FakeGeoParam,  FakeGrad);
  PrintCoordSDFvalues(FakeCoord, FakeGrad, FakePhi);


  printf("Test 4 - Tilted (45°) geometry and negative x component of the coord\n");
  FakeGeoParam.angle = 45;
  (*sdf_func[1])(FakeCoord, FakeGeoParam,  &FakePhi);
  (*sdf_grad_func[1])(FakeCoord, FakeGeoParam,  FakeGrad);
  PrintCoordSDFvalues(FakeCoord, FakeGrad, FakePhi);
}
