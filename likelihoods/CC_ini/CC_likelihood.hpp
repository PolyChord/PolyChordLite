#ifndef CC_likelihood_HPP
#define CC_likelihood_HPP

double loglikelihood (double theta[], int nDims, double phi[], int nDerived);
void prior (double cube[], double theta[], int nDims);
void dumper(int,int,int,double*,double*,double*,double, double);
void setup_loglikelihood();

#endif
