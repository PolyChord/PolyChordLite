#pragma once
#include <string>
#include <vector>
#ifdef USE_MPI
#include "mpi.h"
#endif

struct Settings
{
    int nDims;
    int nDerived;
    int nlive;
    int num_repeats;
    int nprior;
    int nfail;
    bool do_clustering;
    int feedback;
    double precision_criterion;
    double logzero;
    int max_ndead;
    double boost_posterior;
    bool posteriors;
    bool equals;
    bool cluster_posteriors;
    bool write_resume;
    bool write_paramnames;
    bool read_resume;
    bool write_stats;
    bool write_live;
    bool write_dead;
    bool write_prior;
    bool maximise;
    double compression_factor;
    bool synchronous;
    std::string base_dir;
    std::string file_root;
    std::vector<double> grade_frac;
    std::vector<int> grade_dims;
    std::vector<double> loglikes;
    std::vector<int> nlives;
    int seed;

    Settings(int _nDims=0,int _nDerived=0);
};

void run_polychord( 
        double (*c_loglikelihood_ptr)(double*,int,double*,int), 
        void (*c_prior_ptr)(double*,double*,int), 
        void (*c_dumper_ptr)(int,int,int,double*,double*,double*,double,double), 
        void (*c_cluster_ptr)(double*,int*,int,int), 
        Settings s);
void run_polychord( 
        double (*c_loglikelihood_ptr)(double*,int,double*,int), 
        void (*c_prior_ptr)(double*,double*,int), 
        void (*c_dumper_ptr)(int,int,int,double*,double*,double*,double,double), 
        Settings s);
void run_polychord( 
        double (*loglikelihood)(double*,int,double*,int),
        void (*dumper)(int,int,int,double*,double*,double*,double,double), 
        Settings);
void run_polychord( 
        double (*loglikelihood)(double*,int,double*,int),
        void (*prior)(double*,double*,int),
        Settings);
void run_polychord(
        double (*loglikelihood)(double*,int,double*,int),
        Settings);  
void run_polychord( 
        double (*loglikelihood)(double*,int,double*,int), 
        void (*setup_loglikelihood)(), 
        std::string);

#ifdef USE_MPI
void run_polychord( 
        double (*loglikelihood)(double*,int,double*,int), 
        void (*prior)(double*,double*,int), 
        void (*dumper)(int,int,int,double*,double*,double*,double,double), 
        void (*c_cluster_ptr)(double*,int*,int,int), 
        Settings, MPI_Comm &comm);
void run_polychord( 
        double (*loglikelihood)(double*,int,double*,int), 
        void (*prior)(double*,double*,int), 
        void (*dumper)(int,int,int,double*,double*,double*,double,double), 
        Settings, MPI_Comm &comm);
void run_polychord( 
        double (*loglikelihood)(double*,int,double*,int),
        void (*dumper)(int,int,int,double*,double*,double*,double,double), 
        Settings, MPI_Comm &comm);
void run_polychord( 
        double (*loglikelihood)(double*,int,double*,int),
        void (*prior)(double*,double*,int),
        Settings, MPI_Comm &comm);
void run_polychord(
        double (*loglikelihood)(double*,int,double*,int),
        Settings, MPI_Comm &comm);  
void run_polychord( 
        double (*loglikelihood)(double*,int,double*,int), 
        void (*setup_loglikelihood)(), 
        std::string, MPI_Comm &comm);
#endif

double default_loglikelihood(double*,int,double*,int); 
void default_prior(double*,double*,int); 
void default_dumper(int,int,int,double*,double*,double*,double,double);
void default_cluster(double*,int*,int,int);
