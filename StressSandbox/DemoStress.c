#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../SDF-FaultRep/FuzzyFault.h"
#include "../FrictionLaws/Lib_NewmarkTS.h"
#include "../FrictionLaws/Lib_SetOfFrictionLaws.h"
#include "Lib_MiniVoigt.h"


void DoTheEvolution(double loc[],double delta, double displacement[],double *Slip, double deltaTime, double deltaSlip, double ListOfParameters[], double *Theta, double *Traction)
{
  double Normal[2];
  double Tangent[2];
  double velocity[2];
  double displacementHalf[2];
  double SlipHalf=0.0, ThetaHalf=0.0, ThetaDot;
  double SlipDot[2];

  double e_vect[3], sigma[3];
  double lambda = 30.0, G = 24.0;

  double TauC;
  bool UpdateStress;
  //double Traction;

  
  NablaPhi(loc,Normal);
  NormalVecGetTangentVec(Normal,Tangent);
  EvalVelocity(loc,velocity);

  
  /**
   * Step 1. Displacement update 
   *         --> Getting displacement at t+1
   *         --> We leave the velocity as it is (Constant velocity extracted from the field in Phi)
  */
  PartialUpVector(displacement,displacementHalf,2.0*deltaTime,velocity);
  displacement[0] = displacementHalf[0];
  displacement[1] = displacementHalf[1];

  /**
   * Step 2. Extraction of SlipDot from the SDF
   *         --> Getting Slip (Scalar) and State Theta t+1/2
  */
  EvaluateSlipRateAtPoint(loc, delta, SlipDot);
  

  PartialUpScalar(Slip[0], &SlipHalf,deltaTime,SlipDot[0]);

  DotState_AgingLaw(ListOfParameters, SlipDot[0], Theta[0], &ThetaDot);
  
  PartialUpScalar(Theta[0], &ThetaHalf, deltaTime,ThetaDot);

  /**
   * Step 3. Calculate Stress and get Tau
  */
  CalcStrainTensor(e_vect, deltaSlip, displacementHalf);
  CalcStress(lambda, G, e_vect, sigma);

  
  /**
   * Step 3.1 Calculate Tau_Critical
  */
  CompTauCritic(sigma, SlipDot[0], ThetaHalf, ListOfParameters, Normal, &TauC);

  /**
   * Step 3.2 Get tau and compare with 
  */
  GetFaultTraction(sigma, Tangent,Normal, TauC, &Traction[0], &UpdateStress);
  sigma[2]=Traction[0];
  printf("Tau_c: %f - Traction: %f\n",TauC , Traction[0]);
  if(UpdateStress)
  {
    printf("Its Happening!\n");
  }
  /**
   * Step 4 Update of the slip and the state variable
  */
  GetSlipFromTraction(delta, G, UpdateStress, Traction[0], TauC, SlipHalf, Slip);
  
  DotState_AgingLaw(ListOfParameters, SlipDot[0], ThetaHalf, &ThetaDot);
  PartialUpScalar(ThetaHalf, Theta, deltaTime, SlipDot[0]);
  

}



int main(int nargs,char *args[])
{
    FILE *fp;


    double loc[2], displacement[2], Slip;
    double delta;
    double Tau;

    double deltaTime, deltaSlip, Theta_o;
    double ListOfParameters[5];
    double D_c = 0.01;
    double time = 0.0;
    int i;

    deltaTime = 0.001; 
    deltaSlip = 0.00001;
    // Unpacking ListOfParameters = [a, b, mu_o, V_o, D_c]
    ListOfParameters[0] = 0.011 ;
    ListOfParameters[1] = 0.016;
    ListOfParameters[2] = 0.01;
    ListOfParameters[3] = 4.0 * pow(10.0,-9.0);//deltaSlip/(2.0*deltaTime); //
    ListOfParameters[4] = D_c;



    loc[0] = 4.0;
    loc[1] = 15.0;
    delta = 6.0;
    

    displacement[0] = 0.0;
    displacement[1] = 0.0;

    Theta_o=D_c/ListOfParameters[3] ;

    fp = fopen("./Output/Drama.txt","w+");


    printf("=== %s : START! ===\n",__FUNCTION__);


    for (i=1; i<1000; i++)
    {
      time += deltaTime;
      DoTheEvolution(loc, delta, displacement, &Slip, deltaTime, deltaSlip, ListOfParameters, &Theta_o, &Tau);

      fprintf(fp, "%f ; %f ; %f \n", time, Tau, Slip);
    }
    fclose(fp);
    return(0);
}
