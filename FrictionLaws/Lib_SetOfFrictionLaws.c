#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <stdlib.h>
#include "Lib_SetOfFrictionLaws.h"


/**
 * Checking if the value is positive. Otherwise, the program comes to a halt.
*/
void PositiveValueControl(double Variable, char VarName[])
{
  if (Variable < 0.0) 
  {
    printf("Error[%s]: Variable [%s] is not positive. Terminating Program\n",__FUNCTION__,VarName);
    exit(1);
  }
}


/**
 * Different sets of friction laws
 */

/**
 * Static/Dynamic Friction:
 * The Shear stress \tau = \sigma\mu is proportional to the normal stress \sigma.
 * The proportionality constant is \mu. 
 * Unphysical instantaneous stress change from peak value to the sliding value.
*/
void FricSD(double* Fric, double mu)
{
    Fric[0] = mu;
}
void EvalStaticDynamic(double Tau[], double sigma[], double mu)
{
    double Fric;
    FricSD(&Fric, mu);
    Tau[0]=Fric*sigma[0];
}

/**
 * Slip-Weakening (SW)
 * 
 * Shear stress is a decreasing function of slip S up to some distance d_c, 
 * beyond which a constant stress is prescribed.
 * Linear Slip Weakening in the form introduced by Andrews (1976):
 * if s < D_c: Tau = \sigma_n * (\mu_s - (\mu_s - \mu_d) * s / D_c) 
 * Otherwise:  Tau = \sigma_n * \mu_d 
*/

void FricSW(double *Fric, double mu_s, double mu_d, double D_c, double Slip)
{
    PositiveValueControl(Slip,"Slip");
    if (Slip < D_c)
    {
        Fric[0] = (mu_s - (mu_s - mu_d) * Slip / D_c);
    } else {
        Fric[0] = mu_d;
    }

}
void EvalSlipWeakening(double Tau[], double sigma_n[], double mu_s, double mu_d, double D_c,double Slip)
{
    double Fric;
    
    PositiveValueControl(Slip,"Slip");
    FricSW(&Fric, mu_s, mu_d, D_c, Slip);
    Tau[0] = - sigma_n[0] * Fric;
}

/**
 * Rate and State friction (Dieterich-Ruina)
 * 
 * Captures steady state velocity dependence and transient slip and time dependence.
 * \tau = \sigma_n * (\mu_o + a * ln(Sdot / V_o) + b * ln(V_o * \theta / D_c));
*/ 
void FricRS(double *Fric, double Sdot, double Theta, double ListOfParameters[])
{
    // Unpacking ListOfParameters = [a, b, mu_o, V_o, D_c]
    double a = ListOfParameters[0]; // Direct effect coefficient
    double b = ListOfParameters[1]; // Evolution effect coefficient
    double mu_o = ListOfParameters[2]; // Ref. friction coefficient
    double V_o = ListOfParameters[3]; // Ref. slip Velocity
    double D_c = ListOfParameters[4]; // Length scale

    PositiveValueControl(Sdot, "SlipRate");

    Fric[0] = mu_o + a * logf(Sdot / V_o) + b * log((V_o * Theta) / D_c);
}
void EvalRateStateFriction(double Tau[], double sigma_n[], double Sdot, double Theta, double ListOfParameters[])
{
    double Fric;

    PositiveValueControl(Sdot, "SlipRate");
    PositiveValueControl(Theta, "Theta");
    
    FricRS(&Fric, Sdot, Theta, ListOfParameters);
    Tau[0] = - sigma_n[0] * Fric;
    PositiveValueControl(Tau[0], "Tau");
}

/**
 * Modified Rate and State friction law
 * 
 * Captures steady state velocity dependence and transient slip and time dependence.
 * \tau = \sigma_n * a * asinh(( V / (2.0 * V_o)) * exp((\mu_o + b * ln(V_o * \theta / D_c)) / a))
*/ 
void FricModRS(double Fric[], double Sdot, double Theta, double ListOfParameters[])
{
    // Unpacking ListOfParameters = [a, b, mu_o, V_o, D_c]
    double a = ListOfParameters[0]; // Direct effect coefficient
    double b = ListOfParameters[1]; // Evolution effect coefficient
    double mu_o = ListOfParameters[2]; // Ref. friction coefficient
    double V_o = ListOfParameters[3]; // Ref. slip Velocity
    double D_c = ListOfParameters[4]; // Length scale

    PositiveValueControl(Sdot,"SlipRate");
    PositiveValueControl(Theta, "Theta");
    Fric[0] = a * asinhf(( Sdot / (2.0 * V_o)) * expf((mu_o + b * log(V_o * Theta / D_c)) / a));
}

void EvalModRateStateFriction(double Tau[], double sigma_n[], double Sdot, double Theta, double ListOfParameters[])
{
    double Fric;

    PositiveValueControl(Sdot,"SlipRate");
    PositiveValueControl(Theta, "Theta");
    FricModRS(&Fric, Sdot, Theta, ListOfParameters);
    Tau[0] = - sigma_n[0] * Fric;
}


/**
 * State Evolution
*/

/**
 * Aging law
 * \Dot{\theta} = 1 - (\Dot{S} / L) \theta 
*/
void DotState_AgingLaw(double ListOfParameters[], double Sdot, double Theta, double* ThetaDot)
{
    double D_c = ListOfParameters[4]; //Critical Slip weakening Distance

    PositiveValueControl(Sdot,"SlipRate");
    PositiveValueControl(Theta, "Theta");
    ThetaDot[0] = 1.0 - Theta * Sdot / D_c ; 
}

/**
 * Aging law
 * \theta(\theta_o, \Dot{S}, t) = C * exp(-\Dot{S} * t / L) + L / \Dot{S}
*/
void State_AgingLaw(double theta_o, double Sdot, double ListOfParameters[], double time,double* Theta)
{
    double D_c = ListOfParameters[4]; //Critical Slip weakening Distance
    double C;

    PositiveValueControl(Sdot,"SlipRate");
    
    C = theta_o - D_c / Sdot;
    Theta[0] = C * expf(-Sdot * time / D_c) + D_c / Sdot;  // Actually this only would be the case if Sdot was constant
    PositiveValueControl(Theta[0], "Theta");
} 


/**
 * Slip law
 * \Dot{\theta} = - (\Dot{S} / L) \theta * ln((\Dot{S} / L) \theta) 
*/
void DotState_SlipLaw(double ListOfParameters[], double Sdot, double Theta, double* ThetaDot)
{
    double D_c = ListOfParameters[4]; //Critical Slip weakening Distance

    PositiveValueControl(Sdot,"SlipRate");
    PositiveValueControl(Theta, "Theta");
    ThetaDot[0] =  - (Theta * Sdot / D_c) * log(Theta * Sdot / D_c); 
}

/**
 * PRZ law
 * \Dot{\theta} = 1 - ((\Dot{S} / L) \theta )**2
*/
void DotState_PerrinRiceZhengLaw(double ListOfParameters[], double Sdot, double Theta, double* ThetaDot)
{
    double D_c = ListOfParameters[4]; //Critical Slip weakening Distance

    PositiveValueControl(Sdot,"SlipRate");
    PositiveValueControl(Theta, "Theta");
    ThetaDot[0] = 1.0 - pow((Theta * Sdot /(2.0 * D_c)),2.0); 
}