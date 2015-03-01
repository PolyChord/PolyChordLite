# include <math.h>
# include <iostream>
# include "my_cpp_likelihood.hpp"

double sigma;
double mu;


double cpp_loglikelihood (double theta[], int& nDims, double phi[], int& nDerived)
{

    double value = - nDims * ( log(8*atan(1.0))/2 + log(sigma) );

    int i;
    for (i=0;i<nDims;i++) 
    {
        value = value - (theta[i]-mu) * (theta[i]-mu) / sigma / sigma / 2 ;
    }


    for (i=0;i<nDerived;i++)
    {
        phi[i] = i;
    }


    return value;

}

void cpp_loglikelihood_setup ()
{
    sigma = 0.1;
    mu = 0.5;
}
