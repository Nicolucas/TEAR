#ifndef __sdf_faultgeom_h__
#define __sdf_faultgeom_h__


#define CONST_FAULT_ANGLE_DEG 45.0
#define M_PI 3.14159265358979323846264338327950288

typedef struct _SDF *SDF;
typedef struct _GeometryParams *GeometryParams;

struct _GeometryParams {
  double angle;
  double radius;
};

struct _SDF {
  int  type;
  int dim;
  void *data;
  void (*evaluate)(void *, double*, double*);
  void (*evaluate_gradient)(void *, double*,double*);
  void (*evaluate_normal)(double*, SDF, double*);
  void (*evaluate_tangent)(double*, SDF,double*);
  void (*evaluate_DistOnFault)(void *, double*, double*);
};




void SDFCreate(SDF*);
void SDFDestroy(SDF *_s);
void SDFSetup(SDF,int,int);

void GeoParamsCreate(GeometryParams *_g);
void GeoParamsDestroy(GeometryParams *_g);

void SDFEvaluate(SDF s,double c[],double *phi);
void SDFEvaluateGradient(SDF s,double c[],double g[]);
void SDFEvaluateNormal(double c[],SDF s,double n[]);
void SDFEvaluateTangent(double c[],SDF s,double t[]);
void EvaluateDistOnFault(SDF s,double c[],double *distVal);

void FaultSDFNormal(double coor[],SDF s,double n[]);
void FaultSDFTangent(double coor[],SDF s,double t[]);

void horizontal_sdf(void * ctx, double coor[],  double *phi);
void horizontal_grad_sdf(void * ctx, double coor[], double grad[]);
void Horizontal_DistOnFault(void * ctx, double coor[], double *DistOnFault);

void tilted_sdf(void * ctx, double coor[], double *phi);
void tilted_grad_sdf(void * ctx, double coor[], double grad[]);
void Tilted_DistOnFault(void * ctx, double coor[], double *DistOnFault);





#endif