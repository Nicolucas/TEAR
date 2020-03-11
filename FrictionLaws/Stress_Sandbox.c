
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Lib_NewmarkTS.h"


int main () {
    FILE *fp;

    double Tau[1000], Slip[1000], SlipRate[1000];
    double mu_s = 0.6, mu_d = 0.3, D_c = 0.01;
    double sigma[3];
    int i,j;

    double DeltaTime = 0.1;
    double DeltaSlip = 0.0002;

    double ListOfParameters[5] = {0.011, 0.001, 0.2, DeltaSlip/(2.0*DeltaTime), D_c};
    
    double Theta, theta_oo, theta_o;
    
    theta_o = D_c/(ListOfParameters[3]);
    
    sigma[0] = 10, sigma[1] = 5, sigma[2] = 2;

    Slip[0]=0;
    SlipRate[0]=0;



}
