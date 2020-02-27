/**
 * Different sets of friction laws UNTESTED
 */

#include <stdio.h>
#include <stdbool.h>
#include <math.h>



/**
 * Static/Dynamic Friction:
 * The Shear stress \tau = \sigma\mu is proportional to the normal stress \sigma.
 * The proportionality constant is \mu. 
 * Unphysical instantaneous stress change from peak value to the sliding value.
*/
void StaticDynamic(float Tau[], float sigma[], float mu)
{
    Tau[0]=mu*sigma[0];
}
/**
 * Slip-Weakening (SW)
 * Shear stress is a decreasing function of slip S up to some distance d_c, 
 * beyond which a constant stress is prescribed.
 * Max(\tau_d | \tau_p - d_c S)
*/
void SlipWeakening(float Tau[], float tau_d,float tau_p,float D_c,float Slip[])
{
    float DecreasingStress;

    DecreasingStress = tau_p - D_c * Slip[0]; 
    if (DecreasingStress > tau_d)
    {
        Tau[0]=DecreasingStress;
    } else {
        Tau[0]=tau_d;
    }

}

/**
 * Rate and State friction (Dieterich-Ruina)
 * Captures steady state velocity dependence and transient slip and time dependence.
 * \tau = \sigma_n * (\mu_o + a * ln(V / V_o) + b * ln(V_o * \theta / D_c));
*/ 
void RateStateFriction(float Tau[], float sigma_n[], float V, float Theta, float ListOfParameters[])
{
    // Unpacking ListOfParameters = [a, b, mu_o, V_o, D_c]
    float a = ListOfParameters[0];
    float b = ListOfParameters[1];
    float mu_o = ListOfParameters[2]; // Ref. friction coefficient
    float V_o = ListOfParameters[3]; //Ref. slip Velocity
    float D_c = ListOfParameters[4]; //Length scale

    Tau[0] = sigma_n[0] * (mu_o + a * logf(V / V_o) + b * log(V_o * Theta / D_c));
}

/**
 * Modified Rate and State friction law
 * Captures steady state velocity dependence and transient slip and time dependence.
 * \tau = \sigma_n * a * asinh(( V / (2.0 * V_o)) * exp((\mu_o + b * ln(V_o * \theta / D_c)) / a))
*/ 

void ModRateStateFriction(float Tau[], float sigma_n[], float V, float Theta, float ListOfParameters[])
{
    // Unpacking ListOfParameters = [a, b, mu_o, V_o, D_c]
    float a = ListOfParameters[0];
    float b = ListOfParameters[1];
    float mu_o = ListOfParameters[2]; // Ref. friction coefficient
    float V_o = ListOfParameters[3]; //Ref. slip Velocity
    float D_c = ListOfParameters[4]; //Length scale

    Tau[0] = sigma_n[0] * a * asinhf(( V / (2.0 * V_o)) * expf((mu_o + b * log(V_o * Theta / D_c)) / a));
}


/**
 * Evolution Effect
*/

/**
 * Aging law
 * \Dot{\theta} = 1 - (\Dot{S} / L) \theta 
*/
void DotState_AgingLaw(float ListOfParameters[], float Sdot, float* Theta, float* ThetaDot)
{
    float D_c = ListOfParameters[4]; //Length scale

    ThetaDot[0] = Theta[0] * Sdot / D_c ; 
}
