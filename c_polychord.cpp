#include "interfaces.h"

double loglikelihood(double* theta, int nDims, double* phi, int nDerived);

int main()
{

    //double (*c_loglikelihood_ptr)(double[],int&,double[],int&) = loglikelihood;

    int nlive=500;      
    int num_repeats=20;
    bool do_clustering = false;
    int feedback = 1;
    double precision_criterion = 0.0001 ;
    int max_ndead = -1;
    double boost_posterior = 0.0;
    bool posteriors = false ;
    bool equals = false ;
    bool cluster_posteriors = false ;
    bool write_resume = false ;
    bool write_paramnames = false ;
    bool read_resume = false ;
    bool write_stats = false ;
    bool write_live = false ;
    bool write_dead = false;
    int update_files = 500;
    int nDims = 5;
    int nDerived = 0;

    run(loglikelihood, nlive, num_repeats, do_clustering, feedback, precision_criterion, max_ndead, boost_posterior, posteriors, equals, cluster_posteriors, write_resume, write_paramnames, read_resume, write_stats, write_live, write_dead, update_files, nDims, nDerived ) ;

}

double loglikelihood(double* theta, int nDims, double* phi, int nDerived)
{
    double loglike = 0;
    double mu = 0.5;
    double sigma = 10;

    for(int i=0;i<nDims;i++)
        loglike += -(theta[i]-mu)*(theta[i]-mu)/2/sigma/sigma;

    return loglike;

}
