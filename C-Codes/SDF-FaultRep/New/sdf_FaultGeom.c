#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sdf_FaultGeom.h"

/* ==== SDF API ==== */
void SDFCreate(SDF *_s)
{
    SDF s;
    
    s = malloc(sizeof(struct _SDF));
    memset(s,0,sizeof(struct _SDF));
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

void GeoParamsCreate(GeometryParams *_g)
{
    GeometryParams g;
    
    g = malloc(sizeof(struct _GeometryParams));
    memset(g,0,sizeof(struct _GeometryParams));
    *_g = g;
}
void GeoParamsDestroy(GeometryParams *_g)
{
    GeometryParams g;
    
    if (!_g) return;
    g = *_g;
    if (!g) return;
    free(g);
    *_g = NULL;
}


void SDFEvaluate(SDF s,double c[],double *phi)
{
    if (!s->evaluate) {
        printf("Error[SDFEvaluate]: SDF evaluator not set - must call SDFSetup() first\n");
        exit(1);
    }
    
    s->evaluate(s->data, c, phi);
}

void SDFEvaluateGradient(SDF s,double c[],double g[])
{
    if (!s->evaluate_gradient) {
        printf("Error[SDFEvaluateGradient]: SDF gradient valuator not set - must call SDFSetup() first\n");
        exit(1);
    }
    s->evaluate_gradient(s->data, c, g);
}

void SDFEvaluateNormal(double c[],SDF s,double n[])
{
    if (!s->evaluate_normal) {
        printf("Error[SDFEvaluateNormal]: Normal vector evaluator not set - must call SDFSetup() first\n");
        exit(1);
    }
    s->evaluate_normal(c,s,n);
}

void SDFEvaluateTangent(double c[],SDF s,double t[])
{
    if (!s->evaluate_tangent) {
        printf("Error[SDFEvaluateTangent]: Tangent vector evaluator not set - must call SDFSetup() first\n");
        exit(1);
    }
    s->evaluate_tangent(c,s,t);
}

void EvaluateDistOnFault(SDF s,double c[],double *distVal)
{
    if (!s->evaluate_DistOnFault) {
        printf("Error[EvaluateDistOnFault]: Distance on fault function not set  - must call SDFSetup() first\n");
        exit(1);
    }
    s->evaluate_DistOnFault(s->data, c, distVal);
}


/**==================  
 * Normal and Tangent
 * ==================*/

void FaultSDFNormal(double coor[],SDF s,double n[])
{
  double gradphi[2],mag;
  
  SDFEvaluateGradient(s,coor,gradphi);
  mag = sqrt(gradphi[0]*gradphi[0] + gradphi[1]*gradphi[1]);
  n[0] = gradphi[0] / mag;
  n[1] = gradphi[1] / mag;
  
}

void FaultSDFTangent(double coor[],SDF s,double t[])
{
  double gradphi[2],mag;
  
  SDFEvaluateGradient(s,coor,gradphi);
  mag = sqrt(gradphi[0]*gradphi[0] + gradphi[1]*gradphi[1]);
  t[0] = -gradphi[1] / mag;
  t[1] =  gradphi[0] / mag;
  
}


/** ===================================================================
 *                          List of fault geometries
 *  ===================================================================
*/ 

/** 00. Default function, sdf geometry and gradient for a horizontal fault*/
void horizontal_sdf( void * ctx, double coor[], double *phi)
{
  *phi = coor[1];
}

void horizontal_grad_sdf(void * ctx, double coor[], double grad[])
{
  grad[0] = 0.0;
  grad[1] = 1.0;
}

/** Get distance from a coordinate projected onto a tilted (placed here to leave it in a single spot)*/
void Horizontal_DistOnFault(void * ctx, double coor[], double *DistOnFault)
{
  *DistOnFault = coor[0];
}


/** 01. Counterclock-wise Tilted Function: sdf geometry and gradient*/
void tilted_sdf(void * ctx,  double coor[], double *phi)
{
  GeometryParams GeoParamList = (GeometryParams) ctx;
  printf("%f\n",GeoParamList->angle);
  *phi = -sin(GeoParamList->angle* M_PI/180.0) * coor[0] + cos(GeoParamList->angle* M_PI/180.0) * coor[1];
}

void tilted_grad_sdf(void * ctx, double coor[], double grad[])
{
  GeometryParams GeoParamList = (GeometryParams) ctx;
  grad[0] = -sin(GeoParamList->angle* M_PI/180.0);
  grad[1] = cos(GeoParamList->angle* M_PI/180.0);
}

/** Get distance from a coordinate projected onto a tilted (placed here to leave it in a single spot)*/
void Tilted_DistOnFault(void *ctx, double coor[],  double *DistOnFault)
{
  double Fault_angle_deg = *(double*) ctx;
  *DistOnFault = cos(Fault_angle_deg * M_PI/180.0) * coor[0] + sin(Fault_angle_deg* M_PI/180.0) * coor[1];
}


void SDFSetup(SDF s,int dim,int type)
{ 
  
  switch (dim) {
    case 2:
        switch (type) {
            // Horizontal Fault
            case 0:
                printf("Horizontal Fault\n");
                s->evaluate          = horizontal_sdf;
                s->evaluate_gradient = horizontal_grad_sdf;
                s->evaluate_DistOnFault = Horizontal_DistOnFault;
                s->evaluate_normal  = FaultSDFNormal;
                s->evaluate_tangent = FaultSDFTangent;
                s->data = NULL;
                break;
            
            // Tilted Fault
            case 1:
                {
                GeometryParams g;
                GeoParamsCreate(&g);

                printf("Tilted Fault\n");
                g->angle = CONST_FAULT_ANGLE_DEG;
                printf("%f\n",g->angle);
                s->data = (void *) g;

                s->evaluate          = tilted_sdf;
                s->evaluate_gradient = tilted_grad_sdf;
                s->evaluate_DistOnFault = Tilted_DistOnFault;
                s->evaluate_normal  = FaultSDFNormal;
                s->evaluate_tangent = FaultSDFTangent;

                //GeoParamsDestroy(&g);
                }
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
    s->type = type;
}