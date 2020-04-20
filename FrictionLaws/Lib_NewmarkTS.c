#include <stdio.h>
#include <math.h>
#include "Lib_SetOfFrictionLaws.h"
#include "Lib_NewmarkTS.h"

/**PartialUpScalar
 * Partial update for half timestep of a Scalar variable
*/
void PartialUpScalar(double ScalarIn, double *ScalarHalf, double TimeStep, double ScalarDot)
{
    ScalarHalf[0] = ScalarIn + 0.5*TimeStep*ScalarDot;
}

/**PartialUpVector
 * Partial update for half timestep of a 2D vector variable
*/
void PartialUpVector(double VectIn[], double VectHalf[], double TimeStep, double VectDot[])
{
    VectHalf[0] = VectIn[0] + 0.5*TimeStep*VectDot[0];
    VectHalf[1] = VectIn[1] + 0.5*TimeStep*VectDot[1];
}

/**CalcSigmaComponent
 * Calculate the 2D sigma component:
 * Sigma_n1n2 = \sigma_ij * n1_ij * n2_ij 
*/
void CalcSigmaComponent(double Sigma[],double n_i[], double n_j[], double *SigmaScalar)
{
    SigmaScalar[0] = Sigma[0]*n_i[0]*n_j[0] + Sigma[1]*n_i[1]*n_j[1] + Sigma[2]*n_i[1]*n_j[0] + Sigma[2]*n_i[0]*n_j[1];
}

/**CompTauCritic
 * Calculate the Critical Shear Stress from the normal component of the stress \sigma_n
 * - Checks for positive SlipRate, Slip, and Theta, otherwise it throws an error
 * - Checks for negative Sigma_N, if positive it leaves Tau_c = 0
*/
void GetFricValue(double Slip, double SlipDot, double Theta,\
                  double ListOfParameters[], int FricFuncNo, double *Friction)
{
    if (FricFuncNo = 0) // 0 -> LSW 
    {
        FricSW(Friction, ListOfParameters[5], ListOfParameters[6], ListOfParameters[7],  Slip);
    } else if (FricFuncNo = 1) { // 1 -> VW 
        printf("VW Not Implemented\n");
    } else { // 2 -> RSF 
        FricRS(Friction, SlipDot, Theta, ListOfParameters);
    }
}




void CompTauCriticRS(double Sigma[], double Sdot, double Theta, double ListOfParameters[], double n[], double *TauC, double *Friction)
{
    double SigmaN;

    FricRS(Friction, Sdot, Theta, ListOfParameters);
    CalcSigmaComponent(Sigma, n, n, &SigmaN);

    TauC[0] = 0.0;

    if(SigmaN < 0.0)
    {
        TauC[0] = - SigmaN * Friction[0];
    }
}

/**GetFaultTraction
 * Calculate the Fault traction component T and compare with the critical shear traction Tau_c. 
 * The comparison sets a boolean variable true in case of higher T than Tau_c.
*/
void GetFaultTraction(double Sigma[],double n_T[], double n[], double TauC, double *Traction, bool *UpStress)
{
    CalcSigmaComponent(Sigma, n_T, n, &Traction[0]);
    UpStress[0] = false; 

    if (fabs(Traction[0]) > fabs(TauC))
    {
        UpStress[0] = true;
    }  
}

/**GetSlipFromTraction
 * Calculates a new Slip if a boolean is set to True
*/
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

/** GetSlipRate
 * Calculates a slip rate from a change in slip in a half timestep
*/
void GetSlipRate(double OldSlip, double Slip, double TimeStep, double *SlipRate)
{
    SlipRate[0] = 2.0 * (Slip - OldSlip) /  TimeStep;
}
