#include <stdio.h>
#include <stdbool.h>
#include <math.h>


/* function declaration */

float Phi(float *loc);
float *NablaPhi(float *loc);
float *GetProj(float *loc);
float *GetProjDelta(float *loc,float delta);
bool PointInFaut(float *loc,float delta);
float *GetTwins(float *loc);
float *TangentVect(float *loc);
float *TangentVel(float *loc,float *Velocity);
float *CalcSlipRate(float *loc, float *Velocity, float delta);
 
int main () {
    float loc[2];
    float ret;
    float delta;
    bool InFault;
    float *GradPhi;
    float *ProjectedLoc;
    float *TwinPoint;
    float Velocity[2];
    float *tVelDir;
    float *tVel;
    float *SlipRate;

    loc[0] = 4.0;
    loc[1] = 5.0;
    Velocity[0] = 3.0;
    Velocity[1] = 1.0;
    delta = 6.0;
    
    printf("Test Start \n");
    printf( "Location given is : %f, %f\n", loc[0],loc[1] );
    ret = Phi(loc);
    printf( "Signed Distance (Phi) is: %f\n", ret );
    ret = fabs(Phi(loc));
    printf( "Distance (|Phi|) is: %f\n", ret );
    InFault = PointInFaut(loc,delta);
    printf("Is point in the fault? %s\n",InFault ? "true": "false");
    GradPhi = NablaPhi(loc);
    printf( "GradPhi direction is : %f, %f\n", GradPhi[0],GradPhi[1] );
    ProjectedLoc = GetProj(loc);
    printf( "Projected location is : %f, %f\n", ProjectedLoc[0],ProjectedLoc[1] );
    TwinPoint = GetTwins(loc);
    printf( "Twin Point Location is : %f, %f\n", TwinPoint[0],TwinPoint[1] );
    tVelDir = TangentVect(loc);
    printf( "Tangent direction is : %f, %f\n", tVelDir[0],tVelDir[1] );
    tVel = TangentVel(loc,Velocity);
    printf( "Tangent Velocity is : %f, %f\n", tVel[0],tVel[1] );
    SlipRate = CalcSlipRate(loc,Velocity,delta);
    printf( "Slip Rate is : %f, %f\n", SlipRate[0],SlipRate[1] );
    return 0;
}


/*Calculates Phi given a point x,y. Returns scalar*/
float Phi(float *loc) {
    float result;

    result=loc[1]-10.0; /*y - 10*/

    return result; 
}

/*Checks if the location is within the fault*/
bool PointInFaut(float *loc,float delta){
    bool InFault;
    float distLoc;

    distLoc=fabs(Phi(loc)); 

    if(distLoc <= delta){
        InFault=true;
    }
    else{
        InFault=false;
    }
    return InFault;
}
/*Calculates GradPhi and returns an array with direction, normalized by def.*/
float *NablaPhi(float *loc) {
    float GradPhi[2]; 
    float *GradPhi_ptr=GradPhi;

    GradPhi[0]=0.0; /*as in example*/
    GradPhi[1]=1.0;
    return GradPhi_ptr; 
}
/*Calculates the Projected Point given a point*/
float *GetProj(float *loc) {
    float *GradPhi;
    float PhiEval; 
    float Projected[2]; 
    float *Proj=Projected;

    GradPhi=NablaPhi(loc);
    PhiEval=Phi(loc);

    Projected[0]=loc[0]-GradPhi[0]*PhiEval; 
    Projected[1]=loc[1]-GradPhi[1]*PhiEval;
    return Proj; 
}

/*Calculates the Projected Point on the +-Delta Surface given a point*/
float *GetProjDelta(float *loc,float delta) {
    float *GradPhi;
    float PhiEval; 
    float Projected[2]; 
    float *Proj=Projected;
    float DistToDelta; 
    float PhiSign;

    GradPhi = NablaPhi(loc);
    PhiEval = Phi(loc);

    if (PhiEval==0){
        DistToDelta = delta;
    }
    else{
        PhiSign = PhiEval/fabs(PhiEval);
        DistToDelta = (delta*PhiSign - PhiEval);
    }
    

    Projected[0]=loc[0]+GradPhi[0]*DistToDelta; 
    Projected[1]=loc[1]+GradPhi[1]*DistToDelta;
    return Proj; 
}


/*Calculates the Twin Point given a point*/
float *GetTwins(float *loc) {
    float *GradPhi;
    float PhiEval; 
    float Twin[2]; 
    float *Twin_ptr=Twin; 

    GradPhi=NablaPhi(loc);
    PhiEval=Phi(loc);

    Twin[0]=loc[0]-2.0*GradPhi[0]*PhiEval; 
    Twin[1]=loc[1]-2.0*GradPhi[1]*PhiEval;
    return Twin_ptr; 
}

/*Calculate Tangent Vector*/
float *TangentVect(float *loc) {
    float TanDir[2]; 
    float *TanDir_ptr=TanDir; 
    float *GradPhi;

    GradPhi=NablaPhi(loc);

    TanDir[0]=GradPhi[1]; 
    TanDir[1]=-GradPhi[0];
    return TanDir_ptr; 
}

/*Calculate Tangent Velocity*/
float *TangentVel(float *loc,float *Velocity) {
    float TanVel[2];
    float *TanVel_ptr = TanVel;
    float *TanVect;
    float VecDotVel;

    TanVect=TangentVect(loc);
    VecDotVel=Velocity[0]*TanVect[0]+Velocity[1]*TanVect[1];

    TanVel[0]=VecDotVel*TanVect[0]; 
    TanVel[1]=VecDotVel*TanVect[1];
    return TanVel_ptr; 
}

/*Calculate Slip Rate for a location*/
float *CalcSlipRate(float *loc, float *Velocity, float delta){
    bool LocInFault;
    float *PostLocation; 
    float *TwinNew; 
    float TwinVel[2];
    float *VelExtract;
    float *VelPlus;
    float VelP[2];
    float *VelMinus;
    float VelM[2];
    float SlipRate[2];
    float *SlipRate_ptr=SlipRate;


    LocInFault = PointInFaut(loc,delta);

    if(LocInFault){
        PostLocation = GetProjDelta(loc,delta); 
    }
    else{
        PostLocation = loc;
    }
    TwinNew=GetTwins(PostLocation);

    VelPlus = TangentVel(PostLocation,Velocity);
    VelP[0]=(VelPlus[0]);
    VelP[1]=(VelPlus[1]);
    
    /*For the twin vel. it would be good to put it in terms of delta +, flipped 180 deg; basically multiplying by -*/
    TwinVel[0]=-Velocity[0];
    TwinVel[1]=-Velocity[1];    

    VelMinus=TangentVel(TwinNew,TwinVel);
    VelM[0]=(VelMinus[0]);
    VelM[1]=(VelMinus[1]);
    printf( "VelPlus is : %f, %f\n", VelP[0], VelP[1] );
    printf( "VelMinus is : %f, %f\n", VelM[0], VelM[1] );

    SlipRate[0]=(VelP[0])-(VelM[0]);
    SlipRate[1]=(VelP[1])-(VelM[1]);
    return SlipRate_ptr;
}
