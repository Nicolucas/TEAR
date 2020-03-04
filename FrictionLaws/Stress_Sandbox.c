#include <stdio.h>
#include "SetOfFrictionLaws.c"

/**
 * Rotate the stress tensor:
 * \sigma' = R \sigma R^T
 * And rotate back
 * \sigma = R^T \sigma' R
 * 
 * R = (n_T n)
 * where n_T and n are the orthonormal basis of the local coordinates of the fault
 * Thus we get the local stress tensor and think of the 
 * parallel and normal vector as (1 0) and (0 1)
 * 
 * Remember to remember that the stress tensor is expressed in Voigt Notation
*/

void CanToLocRot(float n_T[], float n[]) {
    float CanToLoc[2][2]={n_T[0], n[0], n_T[1], n[1]};
    float CanToLoc_T[2][2]={n_T[0], n_T[1], n[0], n[1]};
}

void CalcSigmaSG(float Sigma[], float Sdot, float Theta, float ListOfParameters[]){
    float Fric;

    FricRS(&Fric, Sdot, Theta, ListOfParameters);
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