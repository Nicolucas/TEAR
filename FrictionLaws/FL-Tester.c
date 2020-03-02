#include <stdio.h>
#include "SetOfFrictionLaws.c"


int main () {
    FILE *fp;
    float Tau[1000], Slip[1000], SlipRate[1000];
    float mu_s = 0.6, mu_d = 0.3, D_c = 0.01;
    float sigma_n[3];
    int i,j;

    float DeltaTime = 0.1;
    float DeltaSlip = 0.0002;

    float ListOfParameters[5] = {0.011, 0.001, 0.2, DeltaSlip/(2.0*DeltaTime), D_c};
    
    float Theta, theta_oo, theta_o;
    
    theta_o = D_c/(ListOfParameters[3]);
    
    sigma_n[0] = 10, sigma_n[1] = 5, sigma_n[2] = 2;
    Slip[0]=0;
    SlipRate[0]=0;

    { // SW Test
        fp = fopen("./Output/TestSW.txt","w+");
        for (i=1; i<1000; i++)
        {
            Slip[i] = Slip[i-1] + DeltaSlip;
            EvalSlipWeakening(&Tau[i],sigma_n, mu_s, mu_d, D_c, &Slip[i]);
            fprintf(fp, "%f ; %f\n", Tau[i], Slip[i]);
        }
        fclose(fp);
    }
    

    theta_oo=theta_o;
    { // SR Test
        fp = fopen("./Output/TestSR.txt","w+");
        for (i=1; i<1000; i++)
        {
            SlipRate[i] = (Slip[i] - Slip[i-1])/DeltaTime;

            State_AgingLaw(theta_o, SlipRate[i] , ListOfParameters, DeltaTime, &Theta);
            theta_o = Theta;
            
            EvalRateStateFriction(&Tau[i], sigma_n, SlipRate[i], Theta, ListOfParameters);
            fprintf(fp, "%f ; %f\n", Tau[i], Slip[i]);
        }
        fclose(fp);
    }
    theta_o=theta_oo;
    { // ModSR Test
        fp = fopen("./Output/TestModSR.txt","w+");
        for (i=1; i<1000; i++)
        {

            State_AgingLaw(theta_o, SlipRate[i] , ListOfParameters, DeltaTime, &Theta);
            theta_o = Theta;
            
            EvalModRateStateFriction(&Tau[i], sigma_n, SlipRate[i], Theta, ListOfParameters);
            fprintf(fp, "%f ; %f\n", Tau[i], Slip[i]);
        }
        fclose(fp);
    }
}