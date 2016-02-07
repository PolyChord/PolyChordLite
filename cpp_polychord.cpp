#include "interfaces.h"
#include <math.h>

const double pi = atan(1.0)*4;
const double mu = 0.5;
const double sigma = 0.1;
const double norm = log(2*pi*sigma*sigma)/2;

double loglikelihood(double* theta, int ndims, double* phi, int nderived);

int main()
{
        const int nDims = 10;
        const int nDerived = 1;

        const int nlive = 500;
        const int num_repeats = nDims*5;
        const bool do_clustering = false;
        const int feedback = 1;
        const double precision_criterion = 0.001;
        const int max_ndead = -1;
        const double boost_posterior = 0.0;
        const bool posteriors = false;
        const bool equals = false;
        const bool cluster_posteriors = false;
        const bool write_resume = false;
        const bool write_paramnames = false;
        const bool read_resume = false;
        const bool write_stats = false;
        const bool write_live = false;
        const bool write_dead = false;
        const int update_files = nlive;

        run( loglikelihood, nlive, num_repeats, do_clustering, feedback, precision_criterion, max_ndead, boost_posterior, posteriors, equals, cluster_posteriors, write_resume, write_paramnames, read_resume, write_stats, write_live, write_dead, update_files, nDims, nDerived ) ;


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
