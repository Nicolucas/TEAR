#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../SDF-FaultRep/FuzzyFault.h"
#include "../FrictionLaws/Lib_NewmarkTS.h"
#include "../FrictionLaws/Lib_SetOfFrictionLaws.h"
#include "../StressSandbox/Lib_MiniVoigt.h"
#include "../StressSandbox/Lib_DisplacementFunctions.h"

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

// [End] Fault Geometry Definition

/**
 * Analytical displacement function for testing 
*/
void Displacement(int FuncNo, double x, double y, double t, double DispVect[], double Grad[], double VelocityVect[])
{
  void (*PresFunArray[])(double, double, double, double[], double[], double[]) =\
    {RectDisplacementFunc, LinearDisplacementFunc, ExpDisplacementFunc};

    if (FuncNo > 2) exit(1);
    (*PresFunArray[FuncNo])(x, y, t, DispVect, Grad, VelocityVect);
}

/**
 * Friction Law: Value Initialization and Selection
 * Options (0-2)-> LSW, VW, RSF
 * 
 * FricParameters = [a, b, mu_o, V_o, L, mu_s, mu_d, D_c]
*/ 

void FLinit_LSW(double FricParam[])
{
    FricParam[5] = 0.6 ; // mu_s
    FricParam[6] = 0.2 ; // mu_d
    FricParam[7] = 2.1 ; // D_c
}

void FLinit_VW(double FricParam[])
{
    FricParam[5] = 0.6 ; // mu_s
    FricParam[6] = 0.2 ; // mu_d
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
                     double LameDelta, double LameG, double *Traction, double Normal[], double Tangent[],\
                     double *Friction, double *Theta, double *Slip, double *SlipDot, double sigma[])
{

    //FricFuncNo: Options (0-2)-> LSW, VW, RSF
    double TauC;

    bool UpdateStress;
    double NewSlip = Slip[0];

    double SigmaN;
    

    /**
     * Step 2. 
     * Calculate the friction coefficient given a friction law
     * 
    */
   
    if (FricFuncNo == 0) // 0 -> LSW 
    {
        GetFricValue(Slip[0], SlipDot[0], Theta[0],\
                 ListOfParameters, FricFuncNo, Friction);   
    } else if (FricFuncNo == 1) { // 1 -> VW 
        printf("VW Not Implemented\n");
        exit(1);
    } else { // 2 -> RSF 
        GetFricValue(Slip[0], SlipDot[0], Theta[0],\
                 ListOfParameters, FricFuncNo, Friction);
    }

    CalcSigmaComponent(sigma, Normal, Normal, &SigmaN);
    
    /**
     * Step 3.
     * Calculate the critical shear traction Tau_Critical
     */
    CompTauCritic(sigma, Normal, &TauC,  Friction[0]);

    /**
     * Step 4.
     * Get the fault-plane traction Tau, and compare with Tau_Critical
     * if Tau > Tau_Critical, update the off-diagonal component of the stress 
     */
    GetFaultTraction(sigma, Tangent, Normal, TauC, Traction, &UpdateStress);
    
    printf("Friction: %f - Tau_C: %f - Traction: %f \n", Friction[0], TauC, Traction[0]);

    if(UpdateStress)
    {
        
        sigma[2] = Traction[0]*(Normal[0]*Tangent[1]+Normal[1]*Tangent[0]);
        NewSlip = Slip[0] + (LameDelta/LameG)*(Traction[0]-TauC);
        Traction[0] = TauC;
    } 
    

    /**
     * Step 5. 
     * Update on the Slip and the Slip rate
    */
    SlipDot[0] = (NewSlip-Slip[0])/deltaTime;
    Slip[0] = NewSlip;
}



void demo1(int FuncNo, int DispFuncNo, char *fn)
{
    FILE *fp;

    
    double loc[2], Normal[2], Tangent[2];
    double displacement[2], velocity[2], Grad[4], e_vect[3], sigma[3];
    double delta = 1.0;

    double lambda = 30.0, G = 10.;

    double time = 0.0, deltaTime = 0.0001;
    
    double ListOfParameters[8];
    double Slip, SlipDot, Friction, Traction, Theta;

    int i;

    loc[0] = 0.001;
    loc[1] = 0.001;
    Theta = 1.0;

    fp = fopen(fn,"w+");
    
    FrictionParametersInitialization(FuncNo, ListOfParameters);

    Displacement(DispFuncNo, loc[0],  loc[1],  time,  displacement,  Grad,  velocity);
    NablaPhi(loc,Normal);
    NormalVecGetTangentVec(Normal,Tangent);
    DotProductTangent2D(displacement, Tangent, &Slip);


    for (i=1; i<100000; i++)
    {
        time += deltaTime;

        /** 
         * Preamble: Geometry
         * get the Normal and the Tangent directions based on the Fault geometry in case it has changed
        */
        NablaPhi(loc,Normal);
        NormalVecGetTangentVec(Normal,Tangent);

        /**
         * Stress components calculation
         * Use the prescribed displacement field to calculate loading of the stress tensor
        */
        Displacement(DispFuncNo, loc[0],  loc[1],  time + deltaTime,  displacement,  Grad,  velocity);

        e_vect[0] = Grad[0];
        e_vect[1] = Grad[1];
        e_vect[2] = (Grad[2] + Grad[3]);
        CalcStress(lambda, G, e_vect, sigma);
        
        // Slip predictor
        Slip = Slip + deltaTime*SlipDot;

        StressCorrector(deltaTime, loc, delta, ListOfParameters, FuncNo,\
                        lambda, G, &Traction, Normal, Tangent,\
                        &Friction, &Theta, &Slip, &SlipDot, sigma);
        
        fprintf(fp, "%f ; %f ; %f ; %f ; %f ; %f \n", time, Slip, SlipDot, Theta, Friction, Traction);
    }
    fclose(fp);
}


int main(int nargs,char *args[])
{
    printf("Running demo: Func 1 \n");
    demo1(0, 0, "./Output/Demo_01.txt");
    printf("Running demo: Func 2 \n");
    demo1(0, 1, "./Output/Demo_02.txt");
    printf("Running demo: Func 3 \n");
    demo1(0, 2, "./Output/Demo_03.txt");

    return(0);
}