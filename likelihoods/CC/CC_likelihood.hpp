#ifndef CC_likelihood_HPP
#define CC_likelihood_HPP
extern "C" {
    double cpp_loglikelihood (double theta[], int& nDims, double phi[], int& nDerived);
    void cpp_prior (double cube[], double theta[], int& nDims);
    void cpp_loglikelihood_setup ();
}
#endif
