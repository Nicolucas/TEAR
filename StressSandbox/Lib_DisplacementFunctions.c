#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void ExpDisplacementFunc(double x, double y, double t, double DispVect[], double Grad[], double VelocityVect[])
{
  double Amp = 2., DiagAmp=0.2;

  DispVect[0] = - Amp*exp(x + y) - DiagAmp*exp(y)*pow(t,2.0);
  DispVect[1] = - Amp*exp(y + x);

  VelocityVect[0] = - DiagAmp*exp(y)*t*2.0;
  VelocityVect[1] = 0.0;  

  Grad[0] = - Amp*exp(x + y); //DU_1/Dx
  Grad[1] = - Amp*exp(x + y); //DU_2/Dy
  Grad[2] = - Amp*exp(x + y) - DiagAmp*exp(y)*pow(t,2.0); //DU_1/Dy
  Grad[3] = - Amp*exp(x + y); //DU_2/Dx
}

void LinearDisplacementFunc(double x, double y, double t, double DispVect[], double Grad[], double VelocityVect[])
{
  double Amp = 2., DiagAmp=0.2;

  DispVect[0] = - Amp*exp(x + y) - DiagAmp*exp(y)*t;
  DispVect[1] = - Amp*exp(y + x);

  VelocityVect[0] = - DiagAmp*exp(y);
  VelocityVect[1] = 0.0;  

  Grad[0] = - Amp*exp(x + y); //DU_1/Dx
  Grad[1] = - Amp*exp(x + y); //DU_2/Dy
  Grad[2] = - Amp*exp(x + y) - DiagAmp*exp(y)*t; //DU_1/Dy
  Grad[3] = - Amp*exp(x + y); //DU_2/Dx
}


void RectDisplacementFunc(double x, double y, double t, double DispVect[], double Grad[], double VelocityVect[])
{
  double Amp = 2., DiagAmp=0.2;


  if (t > 0.4 && t < 0.6){
  DispVect[0] = -Amp*exp(x + y) - DiagAmp*exp(y)*t;
  DispVect[1] = -Amp*exp(y + x);

  VelocityVect[0] = - DiagAmp*exp(y);
  VelocityVect[1] = 0.0;  

  Grad[0] = -Amp*exp(x + y); //DU_1/Dx
  Grad[1] = -Amp*exp(x + y); //DU_2/Dy
  Grad[2] = -Amp*exp(x + y) - DiagAmp*exp(y)*t; //DU_1/Dy
  Grad[3] = -Amp*exp(x + y); //DU_2/Dx
  } else {
  Amp = 0.1, DiagAmp=0.002;

  DispVect[0] = -Amp*exp(x + y) - DiagAmp*exp(y)*t;
  DispVect[1] = -Amp*exp(y + x);

  VelocityVect[0] = - DiagAmp*exp(y);
  VelocityVect[1] = 0.0;  

  Grad[0] = -Amp*exp(x + y); //DU_1/Dx
  Grad[1] = -Amp*exp(x + y); //DU_2/Dy
  Grad[2] = -Amp*exp(x + y) - DiagAmp*exp(y)*t; //DU_1/Dy
  Grad[3] = -Amp*exp(x + y); //DU_2/Dx
  }
}