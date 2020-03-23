
#ifndef __Lib_NewmarkTS_h__
#define __Lib_NewmarkTS_h__

#include <stdbool.h>

void PartialUpScalar(double ScalarIn, double *ScalarHalf, double TimeStep, double ScalarDot);
void PartialUpVector(double VectIn[], double VectHalf[], double TimeStep, double VectDot[]);
void CalcSigmaComponent(double Sigma[],double n_i[], double n_j[], double *SigmaScalar);
void CompTauCritic(double Sigma[], double Sdot, double Theta, double ListOfParameters[], double n[], double *TauC, double *Friction);
void GetFaultTraction(double Sigma[],double n_T[], double n[], double TauC, double *Traction, bool *UpStress);
void GetSlipFromTraction(double delta, double G, bool UpStress, double Traction, double TauC, double OldSlip, double *NewSlip);
void GetSlipRate(double OldSlip, double Slip, double TimeStep, double *SlipRate);

#endif
