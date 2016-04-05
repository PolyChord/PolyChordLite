#include <Python.h>
#include "interfaces.h"
#include <iostream>

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


/* Functions for converting from lists to arrays */
void convert_list(double* array, PyObject* list) {
    int n = PyList_Size(list);
    for (int i=0; i<n; i++) 
        PyList_SET_ITEM(list, i, PyFloat_FromDouble(array[i]));
}

void convert_list(PyObject* list, double* array) {
    PyObject* temp;
    int n = PyList_Size(list);
    for (int i=0; i<n; i++) 
    {
        temp = PyList_GET_ITEM(list, i);
        array[i] = PyFloat_AsDouble(temp);
    }
}

/* Callback to the likelihood and prior */
static PyObject *python_loglikelihood = NULL;

double loglikelihood(double* theta, int nDims, double* phi, int nDerived)
{
    PyObject* list_theta = PyList_New(nDims);
    convert_list(theta,list_theta);

    PyObject* list_phi = PyList_New(nDerived);
    convert_list(phi,list_phi);

    PyObject* pylogL = PyObject_CallFunctionObjArgs(python_loglikelihood,list_theta,list_phi,NULL);

    double logL = PyFloat_AsDouble(pylogL);
    return logL;
}

static PyObject *python_prior = NULL;

void prior(double* cube, double* theta, int nDims)
{
    PyObject* list_cube = PyList_New(nDims);
    convert_list(cube,list_cube);

    PyObject* list_theta = PyObject_CallFunctionObjArgs(python_prior,list_cube,NULL);
    convert_list(list_theta,theta);
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
    char* base_dir;
    char* file_root;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "OOiiiiiididiiiiiiiiiiss", &python_loglikelihood, &python_prior, &nDims, &nDerived, &nlive, &num_repeats, &do_clustering, &feedback, &precision_criterion, &max_ndead, &boost_posterior, &posteriors, &equals, &cluster_posteriors, &write_resume, &write_paramnames, &read_resume, &write_stats, &write_live, &write_dead, &update_files, &base_dir, &file_root))
        return NULL;

    polychord_c_interface( 
            loglikelihood, 
            prior, 
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
    return ret;
}

