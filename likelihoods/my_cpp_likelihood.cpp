# include "my_cpp_likelihood.hpp"

// This module is where your likelihood code should be placed.
//
// * The loglikelihood is called by the subroutine loglikelihood.
// * The likelihood is set up by setup_loglikelihood.
// * You can store any global/saved variables in the module.


//============================================================
// insert likelihood variables here
//
//
//============================================================





// Main loglikelihood function
//
// Either write your likelihood code directly into this function, or call an
// external library from it. This should return the logarithm of the
// likelihood, i.e:
//
// loglikelihood = log_e ( P (data | parameters, model ) )
//
// theta are the values of the input parameters (in the physical space 
// NB: not the hypercube space).
//
// nDims is the size of the theta array, i.e. the dimensionality of the parameter space
//
// phi are any derived parameters that you would like to save with your
// likelihood.
//
// nDerived is the size of the phi array, i.e. the number of derived parameters
//
// The return value should be the loglikelihood
//
// This function is called from likelihoods/fortran_cpp_wrapper.f90
// If you would like to adjust the signature of this call, then you should adjust it there,
// as well as in likelihoods/my_cpp_likelihood.hpp
// 
double cpp_loglikelihood (double theta[], int& nDims, double phi[], int& nDerived)
{
    //============================================================
    // insert likelihood code here
    //
    //
    //============================================================
    
    return 0.0;

}



// Setup of the loglikelihood
// 
// This is called before nested sampling, but after the priors and settings
// have been set up.
// 
// This is the time at which you should load any files that the likelihoods
// need, and do any initial calculations.
// 
// This module can be used to save variables in between calls
// (at the top of the file).
// 
// All MPI threads will call this function simultaneously, but you may need
// to use mpi utilities to synchronise them. This should be done through the
// integer mpi_communicator (which is normally MPI_COMM_WORLD).
//
// This function is called from likelihoods/fortran_cpp_wrapper.f90
// If you would like to adjust the signature of this call, then you should adjust it there,
// as well as in likelihoods/my_cpp_likelihood.hpp
//
void cpp_loglikelihood_setup ()
{
    //============================================================
    // insert likelihood setup here
    //
    //
    //============================================================
}
