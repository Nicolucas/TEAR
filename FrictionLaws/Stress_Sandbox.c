#include <stdio.h>
#include "SetOfFrictionLaws.c"


void PartialUpScalar(float ScalarIn, float ScalarHalf, float TimeStep, float ScalarDot)
{
    ScalarHalf = ScalarIn + 0.5*TimeStep*ScalarDot;
}

void PartialUpVector(float VectIn[], float VectHalf[], float TimeStep, float VectDot[])
{
    VectHalf[0] = VectIn[0] + 0.5*TimeStep*VectDot[0];
    VectHalf[1] = VectIn[1] + 0.5*TimeStep*VectDot[1];
}

void CalcSigmaComponent(float Sigma[],float n_i[], float n_j[], float *SigmaScalar)
{
    SigmaScalar[0] = Sigma[0]*n_i[0]*n_j[0] + Sigma[1]*n_i[1]*n_j[1] + Sigma[2]*n_i[1]*n_j[0] + Sigma[2]*n_i[0]*n_j[1];
}

void CompTauCritic(float Sigma[], float Sdot, float Theta, float ListOfParameters[], float n[], float *TauC)
{
    float Fric;
    float SigmaN;

    FricRS(&Fric, Sdot, Theta, ListOfParameters);
    CalcSigmaComponent(Sigma, n, n, &SigmaN);

    TauC[0] = SigmaN * Fric;

}

void GetFaultTraction(float Sigma[],float n_T[], float n[], float TauC, float *Traction, bool *UpStress)
{
    CalcSigmaComponent(Sigma, n_T, n, &Traction);
    UpStress[0] = false; 

    if (Traction[0] > TauC)
    {
        UpStress[0] = true;
    }  
}

void GetSlipFromTraction(float delta, float G, bool UpStress, float Traction, float TauC, float OldSlip, float *NewSlip)
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
}



int main () {
    FILE *fp;

    float Tau[1000], Slip[1000], SlipRate[1000];
    float mu_s = 0.6, mu_d = 0.3, D_c = 0.01;
    float sigma[3];
    int i,j;

    float DeltaTime = 0.1;
    float DeltaSlip = 0.0002;

    float ListOfParameters[5] = {0.011, 0.001, 0.2, DeltaSlip/(2.0*DeltaTime), D_c};
    
    float Theta, theta_oo, theta_o;
    
    theta_o = D_c/(ListOfParameters[3]);
    
    sigma[0] = 10, sigma[1] = 5, sigma[2] = 2;

    Slip[0]=0;
    SlipRate[0]=0;



}