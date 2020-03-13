#include <stdio.h>
#include <math.h>
#include "Lib_SetOfFrictionLaws.h"
#include "Lib_NewmarkTS.h"


void PartialUpScalar(double ScalarIn, double *ScalarHalf, double TimeStep, double ScalarDot)
{
    ScalarHalf[0] = ScalarIn + 0.5*TimeStep*ScalarDot;
}

void PartialUpVector(double VectIn[], double VectHalf[], double TimeStep, double VectDot[])
{
    VectHalf[0] = VectIn[0] + 0.5*TimeStep*VectDot[0];
    VectHalf[1] = VectIn[1] + 0.5*TimeStep*VectDot[1];
}

void CalcSigmaComponent(double Sigma[],double n_i[], double n_j[], double *SigmaScalar)
{
    SigmaScalar[0] = Sigma[0]*n_i[0]*n_j[0] + Sigma[1]*n_i[1]*n_j[1] + Sigma[2]*n_i[1]*n_j[0] + Sigma[2]*n_i[0]*n_j[1];
}

void CompTauCritic(double Sigma[], double Sdot, double Theta, double ListOfParameters[], double n[], double *TauC)
{
    double Fric;
    double SigmaN;

    FricRS(&Fric, Sdot, Theta, ListOfParameters);
    CalcSigmaComponent(Sigma, n, n, &SigmaN);
    TauC[0] = SigmaN * Fric;
    //printf("%f - %f\n",SigmaN,Fric);
}

void GetFaultTraction(double Sigma[],double n_T[], double n[], double TauC, double *Traction, bool *UpStress)
{
    CalcSigmaComponent(Sigma, n_T, n, &Traction[0]);
    UpStress[0] = false; 

    if (fabs(Traction[0]) > fabs(TauC))
    {
        UpStress[0] = true;
    }  
}

void GetSlipFromTraction(double delta, double G, bool UpStress, double Traction, double TauC, double OldSlip, double *NewSlip)
{

    if (UpStress)
    {
        NewSlip[0] = OldSlip + (Traction - TauC) * delta / G;
    } 
    else
    {
        NewSlip[0] = OldSlip;
    } 
} 

void GetSlipRate(double OldSlip, double Slip, double TimeStep, double *SlipRate)
{
    SlipRate[0] = 2.0 * (Slip - OldSlip) /  TimeStep;
}
