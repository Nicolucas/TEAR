#include <stdio.h>
#include <stdbool.h>
#include <math.h>

/* function prototypes */
void EvalPhi(double loc[],double* phi);
void PointInFaut(double distLoc,double delta,bool *InFault);
void NablaPhi(double GradPhi[]);
void GetProj(double loc[],double GradPhi[], double PhiEval, double Projected[]);
void GetProjDelta(double loc[],double GradPhi[], double PhiEval, double delta,double Projected[]);
void GetTwins(double loc[],double PhiEval,double GradPhi[],double Twin[]);
void TangentVect(double GradPhi[],double TanDir[]);
void TangentVel(double Velocity[],double TanVect[],double TanVel[]);
void EvalVelocity(double PhiEval,double Velocity[]);
void CalcSlipRate(double VelPlus[], double VelMinus[],double PhiPlus,double PhiMinus, bool LocInFault, double SlipRate[]);
void LocateInFault(double loc[], bool LocInFault,double GradPhi[], double PhiEval,double delta, double PostLocation[]);


int main (int nargs, char *args[])
{
    double loc[2];
    double Phi;
    double delta;
    bool  InFault;
    double GradPhi[2];
    double ProjectedLoc[2];
    double TwinPoint[2];
    double TanDir[2];
    
    double Velocity[2];
    double tVel[2];
    double SlipRate[2];

    double PostLocation[2];
    double PostTwin[2];
    double PhiPL;
    double PhiTwin;
    double VelocityOne[2], VelTanOne[2], TanDirOne[2];
    double VelocityTwo[2], VelTanTwo[2], TanDirTwo[2];

  
    loc[0] = 4.0;
    loc[1] = 5.0;
    delta = 6.0;
    
    printf("Test Start \n");
    printf("Location given is : %f, %f\n", loc[0], loc[1]);

    EvalPhi(loc, &Phi);
    printf("Signed distance (Phi) is : %f\n", Phi);
    printf("Distance (|Phi|) is : %f\n", fabsf(Phi));

    PointInFaut(Phi, delta, &InFault);
    printf("Is point in the fault? : %s\n", InFault ? "true": "false");

    NablaPhi(GradPhi);
    printf("GradPhi direction is : %f, %f\n", GradPhi[0], GradPhi[1]);

    GetProj(loc, GradPhi, Phi, ProjectedLoc);
    printf("Projected location is : %f, %f\n", ProjectedLoc[0], ProjectedLoc[1]);

    GetTwins(loc, Phi, GradPhi, TwinPoint);
    printf("Twin point location is : %f, %f\n", TwinPoint[0], TwinPoint[1]);

    TangentVect(GradPhi, TanDir);
    printf("Tangent direction is : %f, %f\n", TanDir[0], TanDir[1]);

    EvalVelocity(Phi, Velocity);
    TangentVel(Velocity, TanDir, tVel);
    printf("Tangent velocity is : %f, %f\n", tVel[0], tVel[1]);


    printf("-----------------------------------------\n");
    LocateInFault(loc, InFault, GradPhi, Phi, delta, PostLocation);
    EvalPhi(PostLocation, &PhiPL);
    GetTwins(PostLocation, PhiPL, GradPhi, PostTwin);
    EvalPhi(PostTwin, &PhiTwin);
    printf("Projected fault coordinate (post) : %f, %f\n", PostLocation[0], PostLocation[1]);
    printf("Projected fault coordinate (twin) : %f, %f\n", PostTwin[0], PostTwin[1]);

    EvalVelocity(PhiPL, VelocityOne);
    TangentVect(GradPhi, TanDirOne);
    TangentVel(VelocityOne, TanDirOne, VelTanOne);

    EvalVelocity(PhiTwin, VelocityTwo);
    TangentVect(GradPhi, TanDirTwo);
    TangentVel(VelocityTwo, TanDirTwo, VelTanTwo);

    CalcSlipRate(VelTanOne, VelTanTwo, PhiPL, PhiTwin, InFault, SlipRate);
    printf("Slip rate is : %f, %f\n", SlipRate[0], SlipRate[1]);
    return(0);
}

/* Calculates Phi given a point x,y. */
void EvalPhi(double loc[], double* phi)
{
    phi[0] = loc[1] - 10.0; /* y - 10 */
}

/* Checks if the location is within the fault */
void PointInFaut(double distLoc, double delta, bool *InFault)
{
    if (fabsf(distLoc) <= delta) {
        *InFault = true;
    } else {
        *InFault = false;
    }
}

/* Calculates GradPhi and returns an array with direction, normalized by def. */
void NablaPhi(double GradPhi[])
{
    GradPhi[0] = 0.0; /* as in example */
    GradPhi[1] = 1.0;
}

/* Calculates the projected point given a point */
void GetProj(double loc[], double GradPhi[], double PhiEval, double Projected[])
{
    Projected[0] = loc[0] - GradPhi[0] * PhiEval;
    Projected[1] = loc[1] - GradPhi[1] * PhiEval;
}

/* Calculates the projected point on the +/-delta surface given a point */
void GetProjDelta(double loc[], double GradPhi[], double PhiEval, double delta,double Projected[])
{
    GetProj(loc, GradPhi, PhiEval, &Projected);

    Projected[0] = Projected[0] + GradPhi[0] * delta;
    Projected[1] = Projected[1] + GradPhi[1] * delta;
}

/* Calculates the twin point given a point */
void GetTwins(double loc[], double PhiEval, double GradPhi[], double Twin[])
{
    Twin[0] = loc[0] - 2.0 * GradPhi[0] * PhiEval;
    Twin[1] = loc[1] - 2.0 * GradPhi[1] * PhiEval;
}

/* Calculate tangent vector */
void TangentVect(double GradPhi[], double TanDir[])
{
    TanDir[0] =  GradPhi[1];
    TanDir[1] = -GradPhi[0];
}

/* Calculate tangent velocity */
void TangentVel(double Velocity[], double TanVect[], double TanVel[])
{
    double VecDotVel;

    VecDotVel = Velocity[0] * TanVect[0] + Velocity[1] * TanVect[1];

    TanVel[0] = VecDotVel * TanVect[0];
    TanVel[1] = VecDotVel * TanVect[1];
}

/* Evaluate velocity */
void EvalVelocity(double PhiEval, double Velocity[])
{
    if (PhiEval < 0) {
        Velocity[0] = 3.0;
        Velocity[1] = 1.0;
    } else if (PhiEval > 0) {
        Velocity[0] = -3.0;
        Velocity[1] =  1.0;
    } else {
        Velocity[0] = 0.0;
        Velocity[1] = 1.0;
    }
}

/* Calculates the twin point given a point */
void LocateInFault(double loc[], bool LocInFault, double GradPhi[], double PhiEval, double delta, double PostLocation[])
{
    if (LocInFault) {
        GetProjDelta(loc, GradPhi, PhiEval, delta, PostLocation);
    } else {
        PostLocation[0] = loc[0];
        PostLocation[1] = loc[1];
    }
}

/* Calculate slip rate for a location */
void CalcSlipRate(double VelPlus[], double VelMinus[], double PhiPlus, double PhiMinus, bool LocInFault, double SlipRate[])
{
    SlipRate[0] = 0.0;
    SlipRate[1] = 0.0;
    if (LocInFault) {
        if (PhiPlus > PhiMinus) {
            SlipRate[0] = VelPlus[0] - VelMinus[0];
            SlipRate[1] = VelPlus[1] - VelMinus[1];
        } else if (PhiPlus < PhiMinus) {
            SlipRate[0] = VelMinus[0] - VelPlus[0];
            SlipRate[1] = VelMinus[1] - VelPlus[1];
        }
    }
}

