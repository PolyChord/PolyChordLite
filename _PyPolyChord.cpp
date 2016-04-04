#include <Python.h>
#include "interfaces.h"
#include <iostream>

double default_loglikelihood(double* theta, int nDims, double* phi, int nDerived)
{
    double logL = 0.0;
    for(int i=0;i<nDims;i++) logL -= theta[i]*theta[i];
    return logL;
};

void default_prior(double* cube, double* theta, int nDims)
{
    for(int i=0;i<nDims;i++) theta[i] = cube[i]; 
};


/* Docstrings */
static char module_docstring[] =
    "PyPolyChord: This module provides a Python interface to PolyChord.";
static char run_docstring[] =
    "Runs PyPolyChord";


/* Available functions */
static PyObject *run_PyPolyChord(PyObject *self, PyObject *args);

static PyMethodDef module_methods[] = {
    {"run", run_PyPolyChord, METH_VARARGS, run_docstring},
    {NULL, NULL, 0, NULL}
};

/* Initialize the module */
PyMODINIT_FUNC init_PyPolyChord(void)
{
    PyObject *m = Py_InitModule3("_PyPolyChord", module_methods, module_docstring);
    if (m == NULL)
        return;
}

/* Function to run PyPolyChord */
static PyObject *run_PyPolyChord(PyObject *self, PyObject *args)
{
    int nDims;
    int nDerived;
    int nlive;
    int num_repeats;
    int do_clustering;
    int feedback;
    double precision_criterion;
    int max_ndead;
    double boost_posterior;
    int posteriors;
    int equals;
    int cluster_posteriors;
    int write_resume;
    int write_paramnames;
    int read_resume;
    int write_stats;
    int write_live;
    int write_dead;
    int update_files;


    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "iiiiiididiiiiiiiiii", &nDims, &nDerived, &nlive, &num_repeats, &do_clustering, &feedback, &precision_criterion, &max_ndead, &boost_posterior, &posteriors, &equals, &cluster_posteriors, &write_resume, &write_paramnames, &read_resume, &write_stats, &write_live, &write_dead, &update_files))
        return NULL;

    char * base_dir = "chains";
    char * file_root = "PyPolyChord_test";

    polychord_c_interface( 
            default_loglikelihood, 
            default_prior, 
            nlive, 
            num_repeats,
            do_clustering,
            feedback,
            precision_criterion,
            max_ndead,
            boost_posterior,
            posteriors,
            equals,
            cluster_posteriors,
            write_resume,
            write_paramnames,
            read_resume,
            write_stats,
            write_live,
            write_dead,
            update_files,
            nDims,
            nDerived,
            base_dir,
            file_root);


    /* Build the output tuple */
    PyObject *ret = Py_None;
    Py_INCREF(Py_None);
    return ret;
}

