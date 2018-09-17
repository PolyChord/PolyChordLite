#ifndef CC_ini_likelihood_HPP
#define CC_ini_likelihood_HPP
#include <string>

double loglikelihood (double theta[], int nDims, double phi[], int nDerived);
void prior (double cube[], double theta[], int nDims);
void dumper(int,int,int,double*,double*,double*,double, double);
void setup_loglikelihood();
void set_ini(std::string);

#endif
