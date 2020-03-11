#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include "FuzzyFault.h"

/* Calculates Phi given a point x,y. */
void EvalPhi(double loc[], double* phi)
{
    phi[0] = loc[1] - 10.0; /* y - 10 */
}

/* Checks if the location is within the fault */
void PointInFaut(double distLoc, double delta, bool *InFault)
{
    if (fabs(distLoc) <= delta) {
        *InFault = true;
    } else {
        *InFault = false;
    }
}

/* Calculates GradPhi and returns an array with direction, normalized by def. */
void NablaPhi(double loc[],double GradPhi[])
{
    GradPhi[0] = 0.0; /* as in example */
    GradPhi[1] = 1.0;
}

/* Calculates the projected point given a point */
void GetProj(double loc[], double GradPhi[], double PhiEval, double Projected[])
{
    Projected[0] = loc[0] - GradPhi[0] * PhiEval;
    Projected[1] = loc[1] - GradPhi[1] * PhiEval;
}

/* Calculates the projected point on the +/-delta surface given a point */
void GetProjDelta(double loc[], double GradPhi[], double PhiEval, double delta,double Projected[])
{
    GetProj(loc, GradPhi, PhiEval, Projected);

    Projected[0] = Projected[0] + GradPhi[0] * delta;
    Projected[1] = Projected[1] + GradPhi[1] * delta;
}

/* Calculates the twin point given a point */
void GetTwins(double loc[], double PhiEval, double GradPhi[], double Twin[])
{
    Twin[0] = loc[0] - 2.0 * GradPhi[0] * PhiEval;
    Twin[1] = loc[1] - 2.0 * GradPhi[1] * PhiEval;
}

/* Calculate tangent vector */
void TangentVect(double GradPhi[], double TanDir[])
{
    TanDir[0] =  GradPhi[1];
    TanDir[1] = -GradPhi[0];
}

/* Calculate tangent velocity */
void TangentVel(double Velocity[], double TanVect[], double TanVel[])
{
    double VecDotVel;

    VecDotVel = Velocity[0] * TanVect[0] + Velocity[1] * TanVect[1];

    TanVel[0] = VecDotVel * TanVect[0];
    TanVel[1] = VecDotVel * TanVect[1];
}

/* Evaluate velocity */
void EvalVelocity(double loc[], double Velocity[])
{
    double PhiEval;
  
    EvalPhi(loc, &PhiEval);
    if (PhiEval < 0) {
        Velocity[0] = 3.0;
        Velocity[1] = 1.0;
    } else if (PhiEval > 0) {
        Velocity[0] = -3.0;
        Velocity[1] =  1.0;
    } else {
        Velocity[0] = 0.0;
        Velocity[1] = 1.0;
    }
}

/* Calculates the twin point given a point */
void LocateInFault(double loc[], bool LocInFault, double GradPhi[], double PhiEval, double delta, double PostLocation[])
{
    if (LocInFault) {
        GetProjDelta(loc, GradPhi, PhiEval, delta, PostLocation);
    } else {
        PostLocation[0] = loc[0];
        PostLocation[1] = loc[1];
    }
}

/* Calculate slip rate for a location */
void CalcSlipRate(double VelPlus[], double VelMinus[], double PhiPlus, double PhiMinus, bool LocInFault, double SlipRate[])
{
    SlipRate[0] = 0.0;
    SlipRate[1] = 0.0;
    if (LocInFault) {
        if (PhiPlus * PhiMinus > 0) {
            printf("[Error] CalcSlipRate(): phi(+) and phi(-) must be of opposite sign!\n");
            printf("[Error] CalcSlipRate(): phi(+) = %+1.4e, phi(-) = %+1.4e\n",PhiPlus,PhiMinus);
            return;
        }
      
        if (PhiPlus > 0) {
            SlipRate[0] = VelPlus[0] - VelMinus[0];
            SlipRate[1] = VelPlus[1] - VelMinus[1];
        } else {
            /* 
             Note: I believe this else case is redundent.
             One should compute slip = V+ - V-. 
             If you always project from any coordinate first to phi = 0,
             and then to phi = delta, we should never have to enter here.
             Something to think about.
             ~ DAM
            */
            SlipRate[0] = VelMinus[0] - VelPlus[0];
            SlipRate[1] = VelMinus[1] - VelPlus[1];
        }
    }
}


void NormalVecGetTangentVec(double n[], double t[])
{
  t[0] =  n[1];
  t[1] = -n[0];
}

void VecDot(double x[], double y[], double *d)
{
  *d = x[0]*y[0] + x[1]*y[1];
}

void VecMag(double x[], double *m)
{
  *m = sqrt(x[0]*x[0] + x[1]*x[1]);
}


void EvaluateSlipRateAtPoint(double coor[],double delta,double *sdot)
{
  double coor_plus[] = {0,0};
  double coor_minus[] = {0,0};
  double vel_plus[] = {0,0};
  double vel_minus[] = {0,0};
  double normal[] = {0,0},tangent[] = {0,0};
  double gradphi[2],phi_p,f_phi_p;
  double vdott_plus,vdott_minus;
  
  sdot[0] = sdot[1] = 0.0;
  
  
  EvalPhi(coor,&phi_p);
  f_phi_p = fabs(phi_p);
  
  /* early return if not inside the fault */
  if (f_phi_p > delta) { return; }
  
  NablaPhi(coor,gradphi);
  {
    double mag;
    
    VecMag(gradphi,&mag);
    normal[0] = gradphi[0] / mag;
    normal[1] = gradphi[1] / mag;
  }
  NormalVecGetTangentVec(normal,tangent);
  

  /* project to phi = 0 */
  coor_plus[0] = coor[0] - phi_p * normal[0];
  coor_plus[1] = coor[1] - phi_p * normal[1];
  
  /* project to phi = delta */
  coor_plus[0] = coor_plus[0] + delta * normal[0];
  coor_plus[1] = coor_plus[1] + delta * normal[1];

  /* project to phi = -delta */
  coor_minus[0] = coor_plus[0] - 2.0 * delta * normal[0];
  coor_minus[1] = coor_plus[1] - 2.0 * delta * normal[1];

  EvalVelocity(coor_plus,vel_plus);
  EvalVelocity(coor_minus,vel_minus);
  
  VecDot(vel_plus,tangent,&vdott_plus);
  VecDot(vel_minus,tangent,&vdott_minus);
  
  *sdot = vdott_plus - vdott_minus;
}
