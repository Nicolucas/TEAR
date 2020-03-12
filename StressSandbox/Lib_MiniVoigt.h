#ifndef __Lib_MiniVoigt_h__
#define __Lib_MiniVoigt_h__

void CalcStrainTensor(double e_vect[], double dx, double u_vect[]);
void CalcStress(double lambda, double mu, double e_vect[], double sigma[]);


#endif
