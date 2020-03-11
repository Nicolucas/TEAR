
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Lib_SetOfFrictionLaws.h"


void RectFunction(double Y[], double y_o, double x, double x1, double x2, double Amplitude)
{
    if (x > x1 && x < x2)
    {
        Y[0] = y_o + Amplitude;
    } else {
        Y[0] = y_o;
    }
}


int main () {
    FILE *fp;
    double Tau[1000], Slip[1000], SlipRate[1000], Fric[1000],time[1000];
    double mu_s = 0.6, mu_d = 0.3, D_c = 0.01;
    double sigma_n[3];
    int i,j;

    double DeltaTime = 0.1;
    double DeltaSlip = 0.0002;
     

    double ListOfParameters[5] = {0.011, 0.016, 0.2, DeltaSlip/(2.0*DeltaTime), D_c};
    
    double Theta, theta_oo, theta_o;
    
    theta_o = D_c/(ListOfParameters[3]);
    
    sigma_n[0] = 10, sigma_n[1] = 5, sigma_n[2] = 2;
    Slip[0]=0.001;
    SlipRate[0]=ListOfParameters[3];
    time[0] = 0.0;

    { // SW Test
        fp = fopen("./Output/TestSW.txt","w+");
        for (i=1; i<1000; i++)
        {
            Slip[i] = Slip[i-1] + SlipRate[0]*DeltaTime;
            time[i] = DeltaTime + time[i-1];
            SlipRate[i] = (Slip[i] - Slip[i-1])/DeltaTime;

            FricSW(&Fric[i], mu_s, mu_d, D_c, &Slip[i]);
            EvalSlipWeakening(&Tau[i],sigma_n, mu_s, mu_d, D_c, &Slip[i]);
            fprintf(fp, "%f ; %f ; %f ; %f  ; %f\n", time[i], Tau[i], Slip[i], SlipRate[i], Fric[i]);
        }
        fclose(fp);
    }
    

    theta_oo=theta_o;
    { // SR Test
        fp = fopen("./Output/TestSR.txt","w+");
        for (i=1; i<1000; i++)
        {
            RectFunction(&SlipRate[i], ListOfParameters[3], i, 200, 300, 0.002);
            Slip[i] = Slip[i-1] + SlipRate[i-1]*DeltaTime;
            

            State_AgingLaw(theta_o, SlipRate[i] , ListOfParameters, DeltaTime, &Theta);
            theta_o = Theta;
            
            FricRS(&Fric[i], SlipRate[i], Theta, ListOfParameters);
            EvalRateStateFriction(&Tau[i], sigma_n, SlipRate[i], Theta, ListOfParameters);
            fprintf(fp, "%f ; %f ; %f ; %f ; %f\n", time[i], Tau[i], Slip[i], SlipRate[i], Fric[i]);
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
            
            FricModRS(&Fric[i], SlipRate[i], Theta, ListOfParameters);
            EvalModRateStateFriction(&Tau[i], sigma_n, SlipRate[i], Theta, ListOfParameters);
            fprintf(fp, "%f ; %f ; %f ; %f ; %f\n", time[i], Tau[i], Slip[i], SlipRate[i], Fric[i]);
        }
        fclose(fp);
    }
}
