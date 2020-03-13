
#ifndef __Lib_SetOfFrictionLaws_h__
#define __Lib_SetOfFrictionLaws_h__

void FricSD(double* Fric, double mu);
void EvalStaticDynamic(double Tau[], double sigma[], double mu);
void FricSW(double *Fric, double mu_s, double mu_d, double D_c, double *Slip);
void EvalSlipWeakening(double Tau[], double sigma_n[], double mu_s, double mu_d, double D_c,double Slip[]);
void FricRS(double *Fric, double Sdot, double Theta, double ListOfParameters[]);
void EvalRateStateFriction(double Tau[], double sigma_n[], double Sdot, double Theta, double ListOfParameters[]);
void FricModRS(double Fric[], double Sdot, double Theta, double ListOfParameters[]);
void EvalModRateStateFriction(double Tau[], double sigma_n[], double Sdot, double Theta, double ListOfParameters[]);
void DotState_AgingLaw(double ListOfParameters[], double Sdot, double Theta, double* ThetaDot);
void State_AgingLaw(double theta_o, double Sdot, double ListOfParameters[], double time,double* Theta);


#endif
