

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sdf.h"

static void SDF_EvaluateNormal2_private(SDF s,double coor[],double n[])
{
  double grad[2],mag;
  
  s->evaluate_gradient(coor,grad);
  mag = sqrt(grad[0]*grad[0] + grad[1]*grad[1]);
  n[0] = grad[0]/mag;
  n[1] = grad[1]/mag;
}

static void SDF_EvaluateTangent2_private(SDF s,double coor[],double t[])
{
  double grad[2],mag;
  
  s->evaluate_gradient(coor,grad);
  mag = sqrt(grad[0]*grad[0] + grad[1]*grad[1]);
  t[0] = -grad[1]/mag;
  t[1] =  grad[0]/mag;
}

/* ==== SDF Implementations ==== */

/*
 Defines the SDF given by
  \phi = y - 10
*/
void EvalPhi_planar_y_2d(double coor[], double *phi)
{
  *phi = coor[1] - 10.0;
}

void EvalGradPhi_planar_y_2d(double coor[],double grad_phi[])
{
  grad_phi[0] = 0.0;
  grad_phi[1] = 1.0;
}


/*
 Defines the SDF given by
   \phi = r - 1.333
 where r = sqrt(x^2 + y^2)
 */
void EvalPhi_circle_2d(double coor[], double *phi)
{
  double r;
  r = sqrt(coor[0]*coor[0]+coor[1]*coor[1]);
  *phi = r - 1.333;
}

void EvalGradPhi_circle_2d(double coor[],double grad_phi[])
{
  double or2;
  or2 = 1.0/sqrt(coor[0]*coor[0]+coor[1]*coor[1]);
  grad_phi[0] = or2 * coor[0];
  grad_phi[1] = or2 * coor[1];
}



/* ==== SDF API ==== */
void SDFCreate(SDF *_s)
{
  SDF s;
  
  s = malloc(sizeof(struct _p_SDF));
  memset(s,0,sizeof(struct _p_SDF));
  *_s = s;
}

void SDFDestroy(SDF *_s)
{
  SDF s;
  
  if (!_s) return;
  s = *_s;
  if (!s) return;
  free(s);
  *_s = NULL;
}

void SDFSetup(SDF s,int dim,int type)
{
  
  switch (dim) {
    case 2:

      s->evaluate_normal  = SDF_EvaluateNormal2_private;
      s->evaluate_tangent = SDF_EvaluateTangent2_private;
      
      switch (type) {
        case 0:
          s->evaluate          = EvalPhi_planar_y_2d;
          s->evaluate_gradient = EvalGradPhi_planar_y_2d;
          break;
        case 1:
          s->evaluate          = EvalPhi_circle_2d;
          s->evaluate_gradient = EvalGradPhi_circle_2d;
          break;
          
        default:
          printf("Error[SDFSetup]: No support for dim = 2. SDF type = %d\n",type);
          exit(1);
          break;
      }
      break;
      
    default:
      printf("Error[SDFSetup]: No support for dimension = %d. dim must be 2\n",dim);
      exit(1);
      break;
  }
  s->dim = dim;
  s->type_id = type;
}

void SDFEvaluate(SDF s,double c[],double *phi)
{
  if (!s->evaluate) {
    printf("Error[SDFEvaluate]: SDF evaluator not set - must call SDFSetup() first\n");
    exit(1);
  }
  s->evaluate(c,phi);
}

void SDFEvaluateGradient(SDF s,double c[],double g[])
{
  if (!s->evaluate_gradient) {
    printf("Error[SDFEvaluateGradient]: SDF gradient valuator not set - must call SDFSetup() first\n");
    exit(1);
  }
  s->evaluate_gradient(c,g);
}

void SDFEvaluateNormal(SDF s,double c[],double n[])
{
  if (!s->evaluate_normal) {
    printf("Error[SDFEvaluateNormal]: Normal vector evaluator not set - must call SDFSetup() first\n");
    exit(1);
  }
  s->evaluate_normal(s,c,n);
}

void SDFEvaluateTangent(SDF s,double c[],double t[])
{
  if (!s->evaluate_tangent) {
    printf("Error[SDFEvaluateTangent]: Tangent vector evaluator not set - must call SDFSetup() first\n");
    exit(1);
  }
  s->evaluate_tangent(s,c,t);
}

/*
 
 // plot signed distance function
 splot "sdf-1.gp" u 1:2:3 w l
 or
 set contour base
 splot "sdf-1.gp" u 1:2:3
 
 // plot normals
 plot "sdf-1.gp" u 1:2:($4/10):($5/10) w vec

 // plot tangents
 plot "sdf-1.gp" u 1:2:($6/10):($7/10) w vec

 NOTE: Don't forget to rescale the gnuplot window to have a 1:1 aspect ratio.
       If you forgot to do this, the normals / tangenets will visually appear 
       to point in the wrong direction
 
*/
static void SDFViewGP_2d(SDF s,double xs[],double xe[],const char filename[])
{
  int ns = 100;
  double dx[2];
  int i,j;
  FILE *fp = NULL;
  
  dx[0] = (xe[0] - xs[0])/((double)(ns-1));
  dx[1] = (xe[1] - xs[1])/((double)(ns-1));
  
  
  fp = fopen(filename,"w");
  if (!fp) {
    printf("Error[SDFViewGP_2d]: Could not open file %s\n",filename);
    exit(1);
  }
  
  for (j=0; j<ns; j++) {
    for (i=0; i<ns; i++) {
      double c[2],n[2],t[2],phi;
      
      c[0] = xs[0] + i * dx[0];
      c[1] = xs[1] + j * dx[1];
      
      SDFEvaluate(s,c,&phi);
      SDFEvaluateNormal(s,c,n);
      SDFEvaluateTangent(s,c,t);
      fprintf(fp,"%+1.6e %+1.6e %+1.6e %+1.6e %+1.6e %+1.6e %+1.6e\n",c[0],c[1],phi,n[0],n[1],t[0],t[1]);
    }
    fprintf(fp,"\n");
  }
  
  fclose(fp);
}

void SDFViewGP(SDF s,double xs[],double xe[],const char filename[])
{
  switch (s->dim) {
    case 2:
      SDFViewGP_2d(s,xs,xe,filename);
      break;
      
    default:
      printf("Error[SDFViewGP]: No support for dimension = %d. dim must be 2\n",s->dim);
      exit(1);
      break;
  }
}

