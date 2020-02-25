#include <stdio.h>
#include <stdbool.h>
#include <math.h>

/* function declaration */
void EvalPhi(float loc[],float* phi);
void PointInFaut(float distLoc,float delta,bool *InFault);
void NablaPhi(float GradPhi[]);
void GetProj(float loc[],float GradPhi[], float PhiEval, float Projected[]);
void GetProjDelta(float loc[],float GradPhi[], float PhiEval, float delta,float Projected[]);
void GetTwins(float loc[],float PhiEval,float GradPhi[],float Twin[]);
void TangentVect(float GradPhi[],float TanDir[]);
void TangentVel(float Velocity[],float TanVect[],float TanVel[]);
void EvalVelocity(float PhiEval,float Velocity[]);
void CalcSlipRate(float VelPlus[], float VelMinus[], float SlipRate[]);
void LocateInFault(float loc[], bool LocInFault,float GradPhi[], float PhiEval,float delta, float PostLocation[]);


int main () {
    float loc[2];
    float Phi;
    float delta;
    bool InFault;
    float GradPhi[2];
    float ProjectedLoc[2];
    float TwinPoint[2];
    float TanDir[2];
    
    float Velocity[2];
    float tVel[2];
    float SlipRate[2];

    loc[0] = 4.0;
    loc[1] = 15.0;
    delta = 6.0;
    
    printf("Test Start \n");
    printf( "Location given is : %f, %f\n", loc[0],loc[1] );

    EvalPhi(loc, &Phi);
    printf( "Signed Distance (Phi) is: %f\n", Phi );
    printf( "Distance (|Phi|) is: %f\n", fabs(Phi) );

    PointInFaut(Phi, delta, &InFault);
    printf("Is point in the fault? %s\n",InFault ? "true": "false");

    NablaPhi(GradPhi);
    printf( "GradPhi direction is : %f, %f\n", GradPhi[0],GradPhi[1] );

    GetProj( loc, GradPhi, Phi, ProjectedLoc);
    printf( "Projected location is : %f, %f\n", ProjectedLoc[0],ProjectedLoc[1] );

    GetTwins( loc, Phi, GradPhi, TwinPoint);
    printf( "Twin Point Location is : %f, %f\n", TwinPoint[0],TwinPoint[1] );

    TangentVect( GradPhi, TanDir);
    printf( "Tangent direction is : %f, %f\n", TanDir[0],TanDir[1] );

    EvalVelocity( Phi, Velocity);
    TangentVel(Velocity, TanDir, tVel);
    printf( "Tangent Velocity is : %f, %f\n", tVel[0],tVel[1] );


    printf("-----------------------------------------\n");
    float PostLocation[2];
    float PostTwin[2];
    LocateInFault( loc,  InFault, GradPhi,  Phi, delta,  PostLocation);
    float PhiPL; 
    EvalPhi(PostLocation, &PhiPL);
    GetTwins( PostLocation, PhiPL, GradPhi, PostTwin);
    float PhiTwin;
    EvalPhi(PostTwin, &PhiTwin);
    printf( "PostLocation is %f, %f\nTwin Point Location is : %f, %f\n", PostLocation[0],PostLocation[1],PostTwin[0],PostTwin[1] );
    
    

    float VelocityOne[2], VelTanOne[2], TanDirOne[2];
    EvalVelocity( PhiPL, VelocityOne);
    TangentVect( GradPhi, TanDirOne);
    TangentVel(VelocityOne, TanDirOne, VelTanOne);

    float VelocityTwo[2], VelTanTwo[2], TanDirTwo[2];
    EvalVelocity( PhiTwin, VelocityTwo);
    TangentVect( GradPhi, TanDirTwo);
    TangentVel(VelocityTwo, TanDirTwo, VelTanTwo);

    CalcSlipRate(VelTanOne, VelTanTwo, SlipRate);
    printf( "Slip Rate is : %f, %f\n", SlipRate[0],SlipRate[1] );
    return 0;
}

/*Calculates Phi given a point x,y. Returns scalar*/
void EvalPhi(float loc[],float* phi) {
    phi[0]=loc[1]-10.0; /*y - 10*/
}

/*Checks if the location is within the fault*/
void PointInFaut(float distLoc,float delta,bool *InFault){
    if(fabs(distLoc) <= delta){
        InFault[0]=true;
    }
    else{
        InFault[0]=false;
    }
}

/*Calculates GradPhi and returns an array with direction, normalized by def.*/
void NablaPhi(float GradPhi[]) { 
    GradPhi[0]=0.0; /*as in example*/
    GradPhi[1]=1.0;
}

/*Calculates the Projected Point given a point*/
void GetProj(float loc[],float GradPhi[], float PhiEval, float Projected[]) {
    Projected[0]=loc[0]-GradPhi[0]*PhiEval; 
    Projected[1]=loc[1]-GradPhi[1]*PhiEval;
}

/*Calculates the Projected Point on the +-Delta Surface given a point*/
void GetProjDelta(float loc[],float GradPhi[], float PhiEval, float delta,float Projected[]) {
    float PhiSign, DistToDelta;

    if (PhiEval==0){
        DistToDelta = delta;
    }
    else{
        PhiSign = PhiEval/fabs(PhiEval);
        DistToDelta = (delta*PhiSign - PhiEval);
    }
    
    Projected[0]=loc[0]+GradPhi[0]*DistToDelta; 
    Projected[1]=loc[1]+GradPhi[1]*DistToDelta; 
}

/*Calculates the Twin Point given a point*/
void GetTwins(float loc[],float PhiEval,float GradPhi[],float Twin[]) {
    Twin[0]=loc[0]-2.0*GradPhi[0]*PhiEval; 
    Twin[1]=loc[1]-2.0*GradPhi[1]*PhiEval;
}

/*Calculate Tangent Vector*/
void TangentVect(float GradPhi[],float TanDir[]) {
    TanDir[0]=GradPhi[1]; 
    TanDir[1]=-GradPhi[0];
}

/*Calculate Tangent Velocity*/
void TangentVel(float Velocity[],float TanVect[],float TanVel[]) {
    float VecDotVel;

    VecDotVel=Velocity[0]*TanVect[0]+Velocity[1]*TanVect[1];

    TanVel[0]=VecDotVel*TanVect[0]; 
    TanVel[1]=VecDotVel*TanVect[1];
}

/*Evaluate Velocity*/
void EvalVelocity(float PhiEval,float Velocity[]) {
    if (PhiEval<0){
        Velocity[0]=3.0;
        Velocity[1]=1.0;
    }
    else if(PhiEval>0){
        Velocity[0]=-3.0;
        Velocity[1]=1.0;
    }
    else{
        Velocity[0]=0.0;
        Velocity[1]=1.0;
    }
}

/*Calculates the Twin Point given a point*/
void LocateInFault(float loc[], bool LocInFault,float GradPhi[], float PhiEval,float delta, float PostLocation[]) {
    if(LocInFault){
        GetProjDelta(loc,GradPhi, PhiEval, delta,PostLocation); 
    }
    else{
        PostLocation[0] = loc[0];
        PostLocation[1] = loc[1];
    }
}

/*Calculate Slip Rate for a location*/
void CalcSlipRate(float VelPlus[], float VelMinus[], float SlipRate[]){ 
    SlipRate[0]=VelPlus[0]-VelMinus[0];
    SlipRate[1]=VelPlus[1]-VelMinus[1];
}

