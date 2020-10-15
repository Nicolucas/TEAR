#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include "FuzzyFault.h"

void demo1(double loc[],double delta)
{
  double Phi;
  bool  InFault;
  double GradPhi[2];
  double ProjectedLoc[2];
  double TwinPoint[2];
  double TanDir[2];
  
  double Velocity[2];
  double tVel[2];
  double SlipRate[2];
  
  double PostLocation[2];
  double PostTwin[2];
  double PhiPL;
  double PhiTwin;
  double VelocityOne[2], VelTanOne[2], TanDirOne[2];
  double VelocityTwo[2], VelTanTwo[2], TanDirTwo[2];

  printf("=== %s ===\n",__FUNCTION__);
  
  printf("Test Start \n");
  printf("Location given is : %f, %f\n", loc[0], loc[1]);
  
  EvalPhi(loc, &Phi);
  printf("Signed distance (Phi) is : %f\n", Phi);
  printf("Distance (|Phi|) is : %f\n", fabs(Phi));
  
  PointInFaut(Phi, delta, &InFault);
  printf("Is point in the fault? : %s\n", InFault ? "true": "false");
  
  NablaPhi(loc,GradPhi);
  printf("GradPhi direction is : %f, %f\n", GradPhi[0], GradPhi[1]);
  
  GetProj(loc, GradPhi, Phi, ProjectedLoc);
  printf("Projected location is : %f, %f\n", ProjectedLoc[0], ProjectedLoc[1]);
  
  GetTwins(loc, Phi, GradPhi, TwinPoint);
  printf("Twin point location is : %f, %f\n", TwinPoint[0], TwinPoint[1]);
  
  TangentVect(GradPhi, TanDir);
  printf("Tangent direction is : %f, %f\n", TanDir[0], TanDir[1]);
  
  EvalVelocity(loc, Velocity);
  TangentVel(Velocity, TanDir, tVel);
  printf("Tangent velocity is : %f, %f\n", tVel[0], tVel[1]);
  
  
  printf("-----------------------------------------\n");
  LocateInFault(loc, InFault, GradPhi, Phi, delta, PostLocation);
  EvalPhi(PostLocation, &PhiPL);
  GetTwins(PostLocation, PhiPL, GradPhi, PostTwin);
  EvalPhi(PostTwin, &PhiTwin);
  printf("Projected fault coordinate (post) : %f, %f\n", PostLocation[0], PostLocation[1]);
  printf("Projected fault coordinate (twin) : %f, %f\n", PostTwin[0], PostTwin[1]);
  
  EvalVelocity(PostLocation, VelocityOne);
  TangentVect(GradPhi, TanDirOne);
  TangentVel(VelocityOne, TanDirOne, VelTanOne);
  
  EvalVelocity(PostTwin, VelocityTwo);
  TangentVect(GradPhi, TanDirTwo);
  TangentVel(VelocityTwo, TanDirTwo, VelTanTwo);
  
  CalcSlipRate(VelTanOne, VelTanTwo, PhiPL, PhiTwin, InFault, SlipRate);
  printf("Slip rate is : %f, %f\n", SlipRate[0], SlipRate[1]);
}

void demo2(double loc[],double delta)
{
  double sdot;
  
  printf("=== %s ===\n",__FUNCTION__);
  EvaluateSlipRateAtPoint(loc,delta,&sdot);
  printf("Slip rate is : %f \n",sdot);
}

int main(int nargs, char *args[])
{
  double loc[2];
  double delta;

  loc[0] = 4.0;
  loc[1] = 5.0;
  delta = 6.0;

  demo1(loc,delta);

  demo2(loc,delta);

  return(0);
}
