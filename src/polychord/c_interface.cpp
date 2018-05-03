#include "interfaces.hpp"
#include "interfaces.h"
#include <iostream>

// Constructor for settings struct
Settings::Settings(int _nDims,int _nDerived): 
    nDims               {_nDims},
    nDerived            {_nDerived},
    nlive               {500}, 
    num_repeats         {nDims*5},
    nprior              {-1},
    do_clustering       {false},
    feedback            {1},
    precision_criterion {0.001},
    logzero             {-1e30},
    max_ndead           {-1},
    boost_posterior     {0.0},
    posteriors          {false},
    equals              {false},
    cluster_posteriors  {false},
    write_resume        {false},
    write_paramnames    {false},
    read_resume         {false},
    write_stats         {false},
    write_live          {false},
    write_dead          {false},
    write_prior         {true},
    compression_factor  {0.36787944117144233},
    base_dir            {"chains"},
    file_root           {"test"},
    grade_frac          {1.0},
    grade_dims          {nDims},
    loglikes            {},
    nlives              {},
    seed                {-1}
{}




void run_polychord( 
        double (*c_loglikelihood_ptr)(double*,int,double*,int), 
        void (*c_prior_ptr)(double*,double*,int), 
        void (*c_dumper_ptr)(int,int,int,double*,double*,double*,double,double), 
        Settings s)
{
    // Ridiculous gubbins for passing strings between C and FORTRAN
    char * base_dir = new char[s.base_dir.size()+1];
    std::copy(s.base_dir.begin(),s.base_dir.end(),base_dir);
    base_dir[s.base_dir.size()] = '\0';

    char * file_root = new char[s.file_root.size()+1];
    std::copy(s.file_root.begin(),s.file_root.end(),file_root);
    file_root[s.file_root.size()] = '\0';

    polychord_c_interface( 
            c_loglikelihood_ptr, 
            c_prior_ptr, 
            c_dumper_ptr, 
            s.nlive, 
            s.num_repeats,
            s.nprior,
            s.do_clustering,
            s.feedback,
            s.precision_criterion,
            s.logzero,
            s.max_ndead,
            s.boost_posterior,
            s.posteriors,
            s.equals,
            s.cluster_posteriors,
            s.write_resume,
            s.write_paramnames,
            s.read_resume,
            s.write_stats,
            s.write_live,
            s.write_dead,
            s.write_prior,
            s.compression_factor,
            s.nDims,
            s.nDerived,
            base_dir,
            file_root,
            s.grade_frac.size(),
            &s.grade_frac[0],
            &s.grade_dims[0],
            s.loglikes.size(),
            &s.loglikes[0],
            &s.nlives[0],
            s.seed
                );

    delete[] base_dir;
    delete[] file_root;
}

void run_polychord( 
        double (*loglikelihood)(double*,int,double*,int),
        void (*dumper)(int,int,int,double*,double*,double*,double,double), 
        Settings S)
{ run_polychord(loglikelihood,default_prior,dumper,S); } 
void run_polychord( 
        double (*loglikelihood)(double*,int,double*,int),
        void (*prior)(double*,double*,int),
        Settings S)
{ run_polychord(loglikelihood,prior,default_dumper,S); } 
void run_polychord(
        double (*loglikelihood)(double*,int,double*,int),
        Settings S)
{ run_polychord(loglikelihood,default_prior,default_dumper,S); } 

void run_polychord( 
        double (*c_loglikelihood_ptr)(double*,int,double*,int), 
        void (*c_setup_loglikelihood_ptr)(), 
        std::string inifile)
{
    // Ridiculous gubbins for passing strings between C and FORTRAN
    char * inifile_c = new char[inifile.size()+1];
    std::copy(inifile.begin(),inifile.end(),inifile_c);
    inifile_c[inifile.size()] = '\0';
    polychord_c_interface_ini( c_loglikelihood_ptr, c_setup_loglikelihood_ptr, inifile_c);
    delete[] inifile_c;
}





void default_prior(double* cube, double* theta, int nDims)
{ for(int i=0;i<nDims;i++) theta[i] = cube[i]; }

void default_dumper(int,int,int,double*,double*,double*,double,double) {}
