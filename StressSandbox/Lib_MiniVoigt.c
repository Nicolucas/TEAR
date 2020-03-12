#include <stdio.h>
#include "Lib_MiniVoigt.h"


void CalcStrainTensor(double e_vect[], double dx, double u_vect[])
{
    e_vect[0] = dx * u_vect[0];
    e_vect[1] = dx * u_vect[1];
    e_vect[2] = dx * u_vect[1] + dx * u_vect[0];
}

void CalcStress(double lambda, double mu, double e_vect[], double sigma[])
{
    double c11 = 2.0 * mu + lambda;
    double c12 = lambda;
    double c21 = lambda;
    double c22 = 2.0 * mu + lambda;
    double c33 = mu;

    sigma[0] = c11 * e_vect[0] + c12 * e_vect[1];
    sigma[1] = c22 * e_vect[1] + c21 * e_vect[0];
    sigma[2] = c33 * e_vect[2];
}