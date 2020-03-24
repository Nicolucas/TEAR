#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void ExpDisplacementFunc(double x, double y, double t, double DispVect[], double Grad[], double VelocityVect[])
{
  DispVect[0] = exp(-x - y) - exp(y)*pow(t,2.0);
  DispVect[1] = exp(-y - x);

  VelocityVect[0] = - exp(y)*t*2.0;
  VelocityVect[1] = 0.0;  

  Grad[0] = -exp(-x - y); //DU_1/Dx
  Grad[1] = -exp(-x - y); //DU_2/Dy
  Grad[2] = -exp(-x - y) - exp(y)*pow(t,2.0); //DU_1/Dy
  Grad[3] = -exp(-x - y); //DU_2/Dx
}

void LinearDisplacementFunc(double x, double y, double t, double DispVect[], double Grad[], double VelocityVect[])
{
  DispVect[0] = exp(-x - y) - exp(y)*t;
  DispVect[1] = exp(-y - x);

  VelocityVect[0] = - exp(y);
  VelocityVect[1] = 0.0;  

  Grad[0] = -exp(-x - y); //DU_1/Dx
  Grad[1] = -exp(-x - y); //DU_2/Dy
  Grad[2] = -exp(-x - y) - exp(y)*t; //DU_1/Dy
  Grad[3] = -exp(-x - y); //DU_2/Dx
}


void RectDisplacementFunc(double x, double y, double t, double DispVect[], double Grad[], double VelocityVect[])
{

  if (t > 0.4 && t < 0.6){
  DispVect[0] = exp(-x - y) - exp(y)*t*2.0;
  DispVect[1] = exp(-y - x);

  VelocityVect[0] = - exp(y)*2.0;
  VelocityVect[1] = 0.0;  

  Grad[0] = -exp(-x - y); //DU_1/Dx
  Grad[1] = -exp(-x - y); //DU_2/Dy
  Grad[2] = -exp(-x - y) - exp(y)*t*2.0; //DU_1/Dy
  Grad[3] = -exp(-x - y); //DU_2/Dx
  } else {
  DispVect[0] = exp(-x - y) - exp(y)*t;
  DispVect[1] = exp(-y - x);

  VelocityVect[0] = - exp(y);
  VelocityVect[1] = 0.0;  

  Grad[0] = -exp(-x - y); //DU_1/Dx
  Grad[1] = -exp(-x - y); //DU_2/Dy
  Grad[2] = -exp(-x - y) - exp(y)*t; //DU_1/Dy
  Grad[3] = -exp(-x - y); //DU_2/Dx
  }
}