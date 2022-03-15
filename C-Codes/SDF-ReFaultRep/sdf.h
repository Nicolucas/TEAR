
#ifndef __sdf_h__
#define __sdf_h__

typedef struct _p_SDF *SDF;

struct _p_SDF {
  int  dim;
  int  type_id;
  void (*evaluate)(double*,double*);
  void (*evaluate_gradient)(double*,double*);
  void (*evaluate_normal)(SDF,double*,double*);
  void (*evaluate_tangent)(SDF,double*,double*);
};

void SDFCreate(SDF*);
void SDFDestroy(SDF *_s);
void SDFSetup(SDF,int,int);
void SDFEvaluate(SDF s,double c[],double *phi);
void SDFEvaluateGradient(SDF s,double c[],double g[]);
void SDFEvaluateNormal(SDF s,double c[],double n[]);
void SDFEvaluateTangent(SDF s,double c[],double t[]);
void SDFViewGP(SDF s,double xs[],double xe[],const char filename[]);

#endif
