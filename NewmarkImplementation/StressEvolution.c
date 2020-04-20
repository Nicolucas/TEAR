#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../SDF-FaultRep/FuzzyFault.h"
#include "../FrictionLaws/Lib_NewmarkTS.h"
#include "../FrictionLaws/Lib_SetOfFrictionLaws.h"

/**
 * Initialization has to be done in a structure object with 
 * the size of all the quadrature points
 * It is updated per timestep
 * 
 * Each quadrature point structure is updated per timestep
 * Hence, it requires an initialization step somewhere
 */

/**
 * Fault Geometry
*/

void DotProductTangent2D(double Vect[], double Tangent[], double *VectComp)
{
  VectComp[0] = Vect[0]*Tangent[0] + Vect[1]*Tangent[1];
}

void GetSlipSlipRate(double delta, double loc[], double Normal[], double Tangent[],\
                     double *Slip, double *SlipDot)
{ 
  EvaluateSlipRateAtPoint(loc, delta, SlipDot);
}

// [End] Fault Geometry Definition

/**
 * Friction Law: Value Initialization and Selection
 * Options (0-2)-> LSW, VW, RSF
 * 
 * FricParameters = [a, b, mu_o, V_o, L, mu_s, mu_d, D_c]
*/ 

void FLinit_LSW(double FricParam[])
{
    FricParam[5] = 0.05 ; // mu_s
    FricParam[6] = 0.02 ; // mu_d
    FricParam[7] = 0.03 ; // D_c
}

void FLinit_VW(double FricParam[])
{
    FricParam[5] = 0.05 ; // mu_s
    FricParam[6] = 0.02 ; // mu_d
    FricParam[3] = 1.0 * pow(10.0,-5.0) ; // V_o
}

void FLinit_RSF(double FricParam[])
{
    FricParam[0] = 0.011; // a
    FricParam[1] = 0.016; // b
    FricParam[2] = 0.54; // mu_o
    FricParam[3] = 1.0 * pow(10.0,-5.0); // V_o
    FricParam[4] = 0.01; // L
}

void FrictionParametersInitialization(int FuncNo, double FricParam[])
{
    void (*FLFunc[])(double[]) =\
        {FLinit_LSW, FLinit_VW, FLinit_RSF};

    if (FuncNo > 2) exit(1);
    
    (*FLFunc[FuncNo])(FricParam);
}

// [End] Friction Law: Value Initialization and Selection

/**
 *  Newmark Time Stepping: Implementation function
*/




void StressCorrector(double deltaTime, double loc[], double delta, double ListOfParameters[], int FricFuncNo,\
                     double LameDelta, double LameG, double *Traction,\
                     double *VarFric, double *Theta, double *Slip, double *SlipDot, double sigma[])
{
    double Normal[2];
    double Tangent[2];

    //FricFuncNo: Options (0-2)-> LSW, VW, RSF
    double Friction, TauC;

    bool UpdateStress;
    double NewSlip = Slip[0];

    /** 
     * Preamble:
     * get the Normal and the Tangent directions based on the Fault geometry in case it has changed
     * 
    */
    NablaPhi(loc,Normal);
    NormalVecGetTangentVec(Normal,Tangent);

    /**
     * Step 1. 
     * Extracting the Slip given the fault representation
    */
    
    GetSlipSlipRate(delta, loc, Normal, Tangent, &Slip, &SlipDot);


    /**
     * Step 2. 
     * Calculate the friction coefficient given a friction law
     * 
    */
   
    if (FricFuncNo = 0) // 0 -> LSW 
    {
        GetFricValue(Slip[0], SlipDot[0], Theta[0],\
                 ListOfParameters, FricFuncNo, &Friction);   
    } else if (FricFuncNo = 1) { // 1 -> VW 
        printf("VW Not Implemented\n");
        exit(1);
    } else { // 2 -> RSF 
        GetFricValue(Slip[0], SlipDot[0], Theta[0],\
                 ListOfParameters, FricFuncNo, &Friction);
    }

    /**
     * Step 3.
     * Calculate the critical shear traction Tau_Critical
     */
    CompTauCritic(sigma, Normal, &TauC,  Friction);

    /**
     * Step 4.
     * Get the fault-plane traction Tau, and compare with Tau_Critical
     * if Tau > Tau_Critical, update the off-diagonal component of the stress 
     */
    GetFaultTraction(sigma, Tangent, Normal, TauC, Traction, &UpdateStress);
    
    if(UpdateStress)
    {
        sigma[2] = Traction[0]*(Normal[0]*Tangent[1]+Normal[1]*Tangent[0]);
        NewSlip = Slip[0] + (LameDelta/LameG)*(Traction[0]-TauC);
    } 
        

    /**
     * Step 5. 
     * Update on the Slip and the Slip rate
    */
    SlipDot[0] = (NewSlip-Slip[0])/deltaTime;
    Slip[0] = NewSlip;
}