
#ifndef __FuzzyFault_h__
#define __FuzzyFault_h__

/* function prototypes */
void EvalPhi(double loc[],double* phi);
void PointInFaut(double distLoc,double delta,bool *InFault);
void NablaPhi(double loc[],double GradPhi[]);
void GetProj(double loc[],double GradPhi[], double PhiEval, double Projected[]);
void GetProjDelta(double loc[],double GradPhi[], double PhiEval, double delta,double Projected[]);
void GetTwins(double loc[],double PhiEval,double GradPhi[],double Twin[]);
void TangentVect(double GradPhi[],double TanDir[]);
void TangentVel(double Velocity[],double TanVect[],double TanVel[]);
void EvalVelocity(double loc[],double Velocity[]);
void CalcSlipRate(double VelPlus[], double VelMinus[],double PhiPlus,double PhiMinus, bool LocInFault, double SlipRate[]);
void LocateInFault(double loc[], bool LocInFault,double GradPhi[], double PhiEval,double delta, double PostLocation[]);
void EvaluateSlipRateAtPoint(double coor[],double delta,double *sdot);

void NormalVecGetTangentVec(double n[], double t[]);
void VecDot(double x[], double y[], double *d);
void VecMag(double x[], double *m);



#endif
