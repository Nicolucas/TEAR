/**
 * Different sets of friction laws
 */

#include <stdio.h>
#include <stdbool.h>
#include <math.h>


/**
 * Univariate Relationships
*/

/**
 * Static/Dynamic Friction:
 * The Shear stress \tau = \sigma\mu is proportional to the normal stress \sigma.
 * The proportionality constant is \mu. 
 * Unphysical instantaneous stress change from peak value to the sliding value.
*/
void StaticDynamic(float Tau[],float mu, float sigma[]){
    Tau[0]=mu*sigma[0];
}
/**
 * Slip-Weakening (SW)
 * Shear stress is a decreasing function of slip S up to some distance d_c, 
 * beyond which a constant stress is prescribed.
 * Max(\tau_d | \tau_p - d_c S)
*/
void SlipWeakening(float Tau[], float tau_d,float tau_p,float d_c,float Slip[]){
    float DecreasingStress;

    DecreasingStress = tau_p - d_c * Slip[0]; 
    if (DecreasingStress > tau_d){
        Tau[0]=DecreasingStress;
    }
    else
    {
        Tau[0]=tau_d;
    }
    

    
}

/**
 *  Rate and State Friction (Dieterich-Ruina)
 * Captures steady state velocity dependence and transient slip and time dependence.
 * Defined in terms of 
*/ 

void RateStateFriction(){

}