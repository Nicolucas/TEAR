#include <stdio.h>
#include <stdbool.h>
#include <math.h>

/* function prototypes */
void EvalPhi(float loc[],float* phi);
void PointInFaut(float distLoc,float delta,bool *InFault);
void NablaPhi(float GradPhi[]);
void GetProj(float loc[],float GradPhi[], float PhiEval, float Projected[]);
void GetProjDelta(float loc[],float GradPhi[], float PhiEval, float delta,float Projected[]);
void GetTwins(float loc[],float PhiEval,float GradPhi[],float Twin[]);
void TangentVect(float GradPhi[],float TanDir[]);
void TangentVel(float Velocity[],float TanVect[],float TanVel[]);
void EvalVelocity(float PhiEval,float Velocity[]);
void CalcSlipRate(float VelPlus[], float VelMinus[],float PhiPlus,float PhiMinus, bool LocInFault, float SlipRate[]);
void LocateInFault(float loc[], bool LocInFault,float GradPhi[], float PhiEval,float delta, float PostLocation[]);


int main (int nargs, char *args[])
{
    float loc[2];
    float Phi;
    float delta;
    bool  InFault;
    float GradPhi[2];
    float ProjectedLoc[2];
    float TwinPoint[2];
    float TanDir[2];
    
    float Velocity[2];
    float tVel[2];
    float SlipRate[2];

    float PostLocation[2];
    float PostTwin[2];
    float PhiPL;
    float PhiTwin;
    float VelocityOne[2], VelTanOne[2], TanDirOne[2];
    float VelocityTwo[2], VelTanTwo[2], TanDirTwo[2];

  
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
void EvalPhi(float loc[], float* phi)
{
    phi[0] = loc[1] - 10.0; /* y - 10 */
}

/* Checks if the location is within the fault */
void PointInFaut(float distLoc, float delta, bool *InFault)
{
    if (fabsf(distLoc) <= delta) {
        *InFault = true;
    } else {
        *InFault = false;
    }
}

/* Calculates GradPhi and returns an array with direction, normalized by def. */
void NablaPhi(float GradPhi[])
{
    GradPhi[0] = 0.0; /* as in example */
    GradPhi[1] = 1.0;
}

/* Calculates the projected point given a point */
void GetProj(float loc[], float GradPhi[], float PhiEval, float Projected[])
{
    Projected[0] = loc[0] - GradPhi[0] * PhiEval;
    Projected[1] = loc[1] - GradPhi[1] * PhiEval;
}

/* Calculates the projected point on the +/-delta surface given a point */
void GetProjDelta(float loc[], float GradPhi[], float PhiEval, float delta,float Projected[])
{
    float PhiSign, DistToDelta;

    if (fabsf(PhiEval) < 1.0e-12) {
        DistToDelta = delta;
    } else {
        PhiSign = PhiEval/fabsf(PhiEval);
        DistToDelta = (delta * PhiSign - PhiEval);
    }
    
    Projected[0] = loc[0] + GradPhi[0] * DistToDelta;
    Projected[1] = loc[1] + GradPhi[1] * DistToDelta;
}

/* Calculates the twin point given a point */
void GetTwins(float loc[], float PhiEval, float GradPhi[], float Twin[])
{
    Twin[0] = loc[0] - 2.0 * GradPhi[0] * PhiEval;
    Twin[1] = loc[1] - 2.0 * GradPhi[1] * PhiEval;
}

/* Calculate tangent vector */
void TangentVect(float GradPhi[], float TanDir[])
{
    TanDir[0] =  GradPhi[1];
    TanDir[1] = -GradPhi[0];
}

/* Calculate tangent velocity */
void TangentVel(float Velocity[], float TanVect[], float TanVel[])
{
    float VecDotVel;

    VecDotVel = Velocity[0] * TanVect[0] + Velocity[1] * TanVect[1];

    TanVel[0] = VecDotVel * TanVect[0];
    TanVel[1] = VecDotVel * TanVect[1];
}

/* Evaluate velocity */
void EvalVelocity(float PhiEval, float Velocity[])
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
void LocateInFault(float loc[], bool LocInFault, float GradPhi[], float PhiEval, float delta, float PostLocation[])
{
    if (LocInFault) {
        GetProjDelta(loc, GradPhi, PhiEval, delta, PostLocation);
    } else {
        PostLocation[0] = loc[0];
        PostLocation[1] = loc[1];
    }
}

/* Calculate slip rate for a location */
void CalcSlipRate(float VelPlus[], float VelMinus[], float PhiPlus, float PhiMinus, bool LocInFault, float SlipRate[])
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

