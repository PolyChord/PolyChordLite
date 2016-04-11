#include "interfaces.hpp"
#include "interfaces.h"
#include <iostream>

// Constructor for settings struct
Settings::Settings(int _nDims,int _nDerived): 
    nDims               {_nDims},
    nDerived            {_nDerived},
    nlive               {500}, 
    num_repeats         {nDims*5},
    do_clustering       {false},
    feedback            {1},
    precision_criterion {0.001},
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
    update_files        {nlive},
    base_dir            {"chains"},
    file_root           {"test"}
{}




void run_polychord( 
        double (*c_loglikelihood_ptr)(double*,int,double*,int), 
        void (*c_prior_ptr)(double*,double*,int), 
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
            s.nlive, 
            s.num_repeats,
            s.do_clustering,
            s.feedback,
            s.precision_criterion,
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
            s.update_files,
            s.nDims,
            s.nDerived,
            base_dir,
            file_root);



    delete[] base_dir;
    delete[] file_root;
}

double default_loglikelihood(double* theta, int nDims, double* phi, int nDerived)
{
    return 0.0;
}

void default_prior(double* cube, double* theta, int nDims)
{
    for(int i=0;i<nDims;i++) theta[i] = cube[i]; 
}

void run_polychord( double (*loglikelihood)(double*,int,double*,int), Settings S) 
{
    run_polychord(loglikelihood,default_prior,S);
} 
