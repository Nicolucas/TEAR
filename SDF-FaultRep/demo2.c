
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sdf.h"


void test_method(double c[],int tid)
{
  SDF sdf;
  double phi,g[2],n[2],t[2];
  
  printf("=== %s (2d) sdf-id %d ===\n",__FUNCTION__,tid);
  SDFCreate(&sdf);
  SDFSetup(sdf,2,tid);
  SDFEvaluate(sdf,c,&phi); printf("%+1.4e %+1.4e => sdf(x) %+1.4e\n",c[0],c[1],phi);
  SDFEvaluateGradient(sdf,c,g); printf("%+1.4e %+1.4e => grad_sdf(x) %+1.4e , %+1.4e\n",c[0],c[1],g[0],g[1]);
  SDFEvaluateNormal(sdf,c,n); printf("%+1.4e %+1.4e => n(x) %+1.4e , %+1.4e\n",c[0],c[1],n[0],n[1]);
  SDFEvaluateTangent(sdf,c,t); printf("%+1.4e %+1.4e => t(x) %+1.4e , %+1.4e\n",c[0],c[1],t[0],t[1]);
  
  SDFDestroy(&sdf);
}

void test_plot(int tid)
{
  SDF sdf;
  char name[256];
  double xs[] = {-2,-2};
  double xe[] = {2,2};
  
  
  SDFCreate(&sdf);
  SDFSetup(sdf,2,tid);

  sprintf(name,"sdf-%d.gp",tid);
  SDFViewGP(sdf,xs,xe,name);
  
  SDFDestroy(&sdf);
}

int main(int nargs,char *args[])
{
  {
    double c[2];
    
    c[0] = 0.3; c[1] = 9.5;
    test_method(c,0);
    
    c[0] = 0.3; c[1] = 10.0;
    test_method(c,0);
    
    c[0] = 0.3; c[1] = 10.5;
    test_method(c,0);
  }
  
  {
    double c[2];

    c[0] = 0.5; c[1] = 0.5;
    test_method(c,1);
    
    c[0] = 1.4; c[1] = 1.4;
    test_method(c,1);
  }
  
  test_plot(1);
  
  return(0);
}
