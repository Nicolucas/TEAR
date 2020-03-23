#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../SDF-FaultRep/FuzzyFault.h"
#include "../FrictionLaws/Lib_NewmarkTS.h"
#include "../FrictionLaws/Lib_SetOfFrictionLaws.h"
#include "Lib_MiniVoigt.h"


void Displacement(double x, double y, double t, double DispVect[], double Grad[], double VelocityVect[])
{
  DispVect[0] = exp(-x - y) - exp(y)*pow(t,2.0)/2.0;
  DispVect[1] = exp(-y - x);

  VelocityVect[0] = - exp(y)*t;
  VelocityVect[1] = 0.0;  

  Grad[0] = -exp(-x - y); //DU_1/Dx
  Grad[1] = -exp(-x - y); //DU_2/Dy
  Grad[2] = -exp(-x - y) - exp(y)*pow(t,2.0)/2.0; //DU_1/Dy
  Grad[3] = -exp(-x - y); //DU_2/Dx
}

void SlipFun(double DispVect[], double Tangent[], double *Slip)
{
  Slip[0] = DispVect[0]*Tangent[0] + DispVect[1]*Tangent[1];
}
void SlipRate(double VelVect[], double Tangent[], double *SlipDot)
{
  SlipDot[0] = VelVect[0]*Tangent[0]+VelVect[1]*Tangent[1];
}


void DoTheEvolution(double loc[],double delta, double displacement[],double *Slip, double deltaTime, double time, double deltaSlip, double ListOfParameters[], double *Theta, double *Traction, double *VarFric)
{
  double Normal[2];
  double Tangent[2];
  double velocity[2];
  //double displacementHalf[2];
  double SlipHalf=0.0, ThetaHalf=0.0, ThetaDot;
  double SlipDot;
  


  double e_vect[3], sigma[3];
  double lambda = 30.0, G = 25.0;

  double TauC;
  bool UpdateStress;

  double Grad[4];

  
  NablaPhi(loc,Normal);
  NormalVecGetTangentVec(Normal,Tangent);

  //EvalVelocity(loc,velocity);

  Displacement( loc[0],  loc[1],  time + deltaTime/2.0,  displacement,  Grad,  velocity);
  SlipDot=velocity[0];

  SlipFun( displacement,  Tangent, &SlipHalf);
  SlipRate( velocity,  Tangent, &SlipDot);

  /**
   * Step 1. Displacement update 
   *         --> Getting displacement at t+1
   *         --> We leave the velocity as it is (Constant velocity extracted from the field in Phi)
  */
  //PartialUpVector(displacement,displacementHalf,2.0*deltaTime,velocity);
  //displacement[0] = displacementHalf[0];
  //displacement[1] = displacementHalf[1];

  /**
   * Step 2. Extraction of SlipDot from the SDF
   *         --> Getting Slip (Scalar) and State Theta t+1/2
  */
  //EvaluateSlipRateAtPoint(loc, delta, SlipDot);
  

  //PartialUpScalar(Slip[0], &SlipHalf,deltaTime,SlipDot[0]);

  DotState_AgingLaw(ListOfParameters, SlipDot, Theta[0], &ThetaDot);
  
  PartialUpScalar(Theta[0], &ThetaHalf, deltaTime,ThetaDot);

  /**
   * Step 3. Calculate Stress and get Tau
  */
  e_vect[0] = Grad[0];
  e_vect[1] = Grad[1];
  e_vect[2] = (Grad[2] + Grad[3]);
  CalcStress(lambda, G, e_vect, sigma);
  

  printf("Sigma_n : %f - ", sigma[0]);
  /**
   * Step 3.1 Calculate Tau_Critical
  */
  CompTauCritic(sigma, SlipDot, ThetaHalf, ListOfParameters, Normal, &TauC, VarFric);
  printf("Friction: %f \n",VarFric[0]);
  //VarFri[0]=VarFric;

  /**
   * Step 3.2 Get tau and compare with 
  */
  
  GetFaultTraction(sigma, Tangent,Normal, TauC, Traction, &UpdateStress);
  sigma[2] = Traction[0];
  printf("Time: %f - Tau_c: %f - Traction: %f\n",time ,TauC , Traction[0]);
  if(UpdateStress)
  {
    printf("Its Happening!\n");
  }
  /**
   * Step 4 Update of the slip and the state variable
  */
  GetSlipFromTraction(delta, G, UpdateStress, Traction[0], TauC, SlipHalf, Slip);
  
  DotState_AgingLaw(ListOfParameters, SlipDot, ThetaHalf, &ThetaDot);
  PartialUpScalar(ThetaHalf, &Theta, deltaTime, SlipDot);
}



int main(int nargs,char *args[])
{
    FILE *fp;


    double loc[2], displacement[2], Slip;
    double delta;
    double Tau;
    double VarFr;

    double deltaTime, deltaSlip, Theta_o;
    double ListOfParameters[5];
    double time = 0.001;
    int i;

    deltaTime = 0.0001; 
    deltaSlip = 0.00005;
    // Unpacking ListOfParameters = [a, b, mu_o, V_o, L]
    ListOfParameters[0] = 0.011 ;
    ListOfParameters[1] = 0.016;
    ListOfParameters[2] = 0.48;
    ListOfParameters[3] = deltaSlip/(2.0*deltaTime); //4.0 * pow(10.0,-9.0);//
    ListOfParameters[4] = deltaSlip;



    loc[0] = 0.00001;
    loc[1] = 0.00001;
    delta = 6.0;
    


    Theta_o=ListOfParameters[4]/ListOfParameters[3] ;

    fp = fopen("./Output/Drama.txt","w+");


    printf("=== %s : START! ===\n",__FUNCTION__);


    for (i=1; i<1000; i++)
    {
      time += deltaTime;
      DoTheEvolution(loc, delta, displacement, &Slip, deltaTime, time, deltaSlip, ListOfParameters, &Theta_o, &Tau, &VarFr);

      fprintf(fp, "%f ; %f ; %f ; %f ; %f \n", time, Tau, Slip, Theta_o, VarFr);
    }
    fclose(fp);
    return(0);
}
