#include "interfaces.h"
#include <math.h>

const double pi = atan(1.0)*4;
const double mu = 0.5;
const double sigma = 0.1;
const double norm = log(2*pi*sigma*sigma)/2;

double loglikelihood(double* theta, int ndims, double* phi, int nderived);
void prior(double* cube, double* theta, int ndims);

int main()
{
    int nDims = 20;
    int nDerived = 1;

    Settings settings(nDims,nDerived);

    settings.file_root = "c_test";

    run_polychord(loglikelihood, prior, settings) ;
}



double loglikelihood(double* theta, int ndims, double* phi, int nderived)
{
    double logl = norm*ndims;
    double radius_squared = 0.0;

    for(int i=0;i<ndims;i++)
        radius_squared -= (theta[i]-mu)*(theta[i]-mu)/2/sigma/sigma;

    logl = radius_squared/2/sigma/sigma;
    phi[0] = radius_squared;

    return logl;
}

void prior(double* cube, double* theta, int ndims)
{
    for( int i=0; i<ndims; i++) theta[i] = cube[i];
}
