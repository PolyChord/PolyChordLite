#ifndef INTERFACES_H
#define INTERFACES_H
#include <string>

extern "C" void polychord_c_interface( double (*)(double*,int,double*,int), void (*)(double*,double*,int), int, int, bool, int, double, int, double, bool, bool, bool, bool, bool, bool, bool, bool, bool, int, int, int, char*, char*); 

struct Settings
{
    int nDims;
    int nDerived;
    int nlive;
    int num_repeats;
    bool do_clustering;
    int feedback;
    double precision_criterion;
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
    int update_files;
    std::string base_dir;
    std::string file_root;

    Settings(int _nDims,int _nDerived);
};

void run_polychord( double (*loglikelihood)(double*,int,double*,int), void (*prior)(double*,double*,int), Settings);
void run_polychord( double (*loglikelihood)(double*,int,double*,int), Settings);  

double default_loglikelihood(double*,int,double*,int); 
void default_prior(double*,double*,int); 




#endif
