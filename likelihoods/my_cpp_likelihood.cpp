# include <math.h>
# include "my_cpp_likelihood.hpp"

double sigma;
double mu;


double cpp_loglikelihood (double theta[], double phi[])
{

    int    nDims = sizeof(theta)/sizeof(*theta) + 1;
    double value = - nDims * ( log(8*atan(1.0))/2 + log(sigma) );

    int i;
    for (i=0;i<nDims;i++) 
    {
        value = value - (theta[i]-mu) * (theta[i]-mu) / sigma / sigma / 2 ;
    }
    return value;

}

void cpp_loglikelihood_setup ()
{
    sigma = 0.1;
    mu = 0.5;
}
