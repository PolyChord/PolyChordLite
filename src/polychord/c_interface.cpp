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
    nfail               {-1},
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
    maximise            {true},
    compression_factor  {0.36787944117144233},
    synchronous         {true},
    base_dir            {"chains"},
    file_root           {"test"},
    grade_frac          {1.0},
    grade_dims          {nDims},
    loglikes            {},
    nlives              {},
    seed                {-1}
{}




#ifdef USE_MPI
void run_polychord( 
        double (*c_loglikelihood_ptr)(double*,int,double*,int), 
        void (*c_prior_ptr)(double*,double*,int), 
        void (*c_dumper_ptr)(int,int,int,double*,double*,double*,double,double), 
        void (*c_cluster_ptr)(double*,int*,int,int),
        Settings s, MPI_Comm& comm)
#else
void run_polychord( 
        double (*c_loglikelihood_ptr)(double*,int,double*,int), 
        void (*c_prior_ptr)(double*,double*,int), 
        void (*c_dumper_ptr)(int,int,int,double*,double*,double*,double,double), 
        void (*c_cluster_ptr)(double*,int*,int,int),
        Settings s)
#endif
{
    // Ridiculous gubbins for passing strings between C and FORTRAN
    char * base_dir = new char[s.base_dir.size()+1];
    std::copy(s.base_dir.begin(),s.base_dir.end(),base_dir);
    base_dir[s.base_dir.size()] = '\0';

    char * file_root = new char[s.file_root.size()+1];
    std::copy(s.file_root.begin(),s.file_root.end(),file_root);
    file_root[s.file_root.size()] = '\0';

#ifdef USE_MPI
    MPI_Fint fortran_comm = MPI_Comm_c2f(comm);
#else
	 int fortran_comm = 0;
#endif

    polychord_c_interface( 
            c_loglikelihood_ptr, 
            c_prior_ptr, 
            c_dumper_ptr, 
            c_cluster_ptr, 
            s.nlive, 
            s.num_repeats,
            s.nprior,
            s.nfail,
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
            s.maximise,
            s.compression_factor,
            s.synchronous,
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
            s.seed,
				fortran_comm
                );

    delete[] base_dir;
    delete[] file_root;
}

#ifdef USE_MPI
// In this function no MPI communicator is given, so a communicator is automatically
// made from MPI_COMM_WORLD
void run_polychord( 
        double (*c_loglikelihood_ptr)(double*,int,double*,int), 
        void (*c_prior_ptr)(double*,double*,int), 
        void (*c_dumper_ptr)(int,int,int,double*,double*,double*,double,double), 
        void (*c_cluster_ptr)(double*,int*,int,int), 
        Settings s)
{
	int flag;
	MPI_Initialized(&flag);
	if (flag==0) MPI_Init(NULL,NULL);
	MPI_Comm world_comm;
	MPI_Comm_dup(MPI_COMM_WORLD,&world_comm);
	run_polychord(c_loglikelihood_ptr,c_prior_ptr,c_dumper_ptr,c_cluster_ptr,s,world_comm);
	if (flag==0) MPI_Finalize();
}
void run_polychord( 
        double (*loglikelihood)(double*,int,double*,int),
        void (*prior)(double*,double*,int), 
        void (*dumper)(int,int,int,double*,double*,double*,double,double), 
        Settings S, MPI_Comm &comm)
{ run_polychord(loglikelihood,prior,dumper,default_cluster,S,comm); } 
void run_polychord( 
        double (*loglikelihood)(double*,int,double*,int),
        void (*dumper)(int,int,int,double*,double*,double*,double,double), 
        Settings S, MPI_Comm &comm)
{ run_polychord(loglikelihood,default_prior,dumper,default_cluster,S,comm); } 
void run_polychord( 
        double (*loglikelihood)(double*,int,double*,int),
        void (*prior)(double*,double*,int),
        Settings S, MPI_Comm &comm)
{ run_polychord(loglikelihood,prior,default_dumper,default_cluster,S,comm); } 
void run_polychord(
        double (*loglikelihood)(double*,int,double*,int),
        Settings S, MPI_Comm &comm)
{ run_polychord(loglikelihood,default_prior,default_dumper,default_cluster,S,comm); } 
#endif

void run_polychord( 
        double (*loglikelihood)(double*,int,double*,int),
        void (*prior)(double*,double*,int), 
        void (*dumper)(int,int,int,double*,double*,double*,double,double), 
        Settings S)
{ run_polychord(loglikelihood,prior,dumper,default_cluster,S); } 
void run_polychord( 
        double (*loglikelihood)(double*,int,double*,int),
        void (*dumper)(int,int,int,double*,double*,double*,double,double), 
        Settings S)
{ run_polychord(loglikelihood,default_prior,dumper,default_cluster,S); } 
void run_polychord( 
        double (*loglikelihood)(double*,int,double*,int),
        void (*prior)(double*,double*,int),
        Settings S)
{ run_polychord(loglikelihood,prior,default_dumper,default_cluster,S); } 
void run_polychord(
        double (*loglikelihood)(double*,int,double*,int),
        Settings S)
{ run_polychord(loglikelihood,default_prior,default_dumper,default_cluster,S); } 


#ifdef USE_MPI
void run_polychord( 
        double (*c_loglikelihood_ptr)(double*,int,double*,int), 
        void (*c_setup_loglikelihood_ptr)(), 
        std::string inifile, MPI_Comm &comm)
#else
void run_polychord( 
        double (*c_loglikelihood_ptr)(double*,int,double*,int), 
        void (*c_setup_loglikelihood_ptr)(), 
        std::string inifile)
#endif
{
    // Ridiculous gubbins for passing strings between C and FORTRAN
    char * inifile_c = new char[inifile.size()+1];
    std::copy(inifile.begin(),inifile.end(),inifile_c);
    inifile_c[inifile.size()] = '\0';
#ifdef USE_MPI
    MPI_Fint fortran_comm = MPI_Comm_c2f(comm);
#else
	 int fortran_comm = 0;
#endif
    polychord_c_interface_ini( c_loglikelihood_ptr, c_setup_loglikelihood_ptr, inifile_c, fortran_comm);
    delete[] inifile_c;
}

#ifdef USE_MPI
// In this function no MPI communicator is given, so a communicator is automatically
// made from MPI_COMM_WORLD
void run_polychord( 
        double (*c_loglikelihood_ptr)(double*,int,double*,int), 
        void (*c_setup_loglikelihood_ptr)(), 
        std::string inifile)
{
	int flag;
	MPI_Initialized(&flag);
	if (flag==0) MPI_Init(NULL,NULL);
	MPI_Comm world_comm;
	MPI_Comm_dup(MPI_COMM_WORLD,&world_comm);
	run_polychord(c_loglikelihood_ptr,c_setup_loglikelihood_ptr,inifile,world_comm);
	if (flag==0) MPI_Finalize();
}
#endif

void default_prior(double* cube, double* theta, int nDims)
{ for(int i=0;i<nDims;i++) theta[i] = cube[i]; }

void default_dumper(int,int,int,double*,double*,double*,double,double) {}

void default_cluster(double* points, int* cluster_list, int nDims, int nPoints)
{ for(int i=0;i<nPoints;i++) cluster_list[i] = 0; }
