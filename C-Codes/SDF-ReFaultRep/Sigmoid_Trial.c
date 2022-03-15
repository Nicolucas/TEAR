#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


void sigmoidPhi(double loc[],double *Phi)
{
    *Phi = -loc[1] + 1.0/(1.0 + exp(-loc[0]));
}


void sigmoidNablaPhi(double loc[],double GradPhi[])
{
    GradPhi[0] = exp(loc[0])/(1.0 + exp(loc[0]));
    GradPhi[1] = -1.0;
}


void SDFEvaluate(double loc[],double GradPhi[])
{
    
}