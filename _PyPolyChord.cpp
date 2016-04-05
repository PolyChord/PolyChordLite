#include <Python.h>
#include "interfaces.h"
#include <iostream>

//double default_loglikelihood(double* theta, int nDims, double* phi, int nDerived)
//{
//    double logL = 0.0;
//    for(int i=0;i<nDims;i++) logL -= theta[i]*theta[i];
//    return logL;
//}

void default_prior(double* cube, double* theta, int nDims)
{
    for(int i=0;i<nDims;i++) theta[i] = cube[i]; 
}


/* Docstrings */
static char module_docstring[] =
    "PyPolyChord: This module provides a Python interface to PolyChord.";
static char run_docstring[] =
    "Runs PyPolyChord";
static char set_loglikelihood_docstring[] =
    "Sets the loglikelihood function";


/* Available functions */
static PyObject *run_PyPolyChord(PyObject *self, PyObject *args);
static PyObject *set_loglikelihood_PyPolyChord(PyObject *self, PyObject *args);

static PyMethodDef module_methods[] = {
    {"run", run_PyPolyChord, METH_VARARGS, run_docstring},
    {"set_loglikelihood", set_loglikelihood_PyPolyChord, METH_VARARGS, set_loglikelihood_docstring},
    {NULL, NULL, 0, NULL}
};


/* Initialize the module */
PyMODINIT_FUNC init_PyPolyChord(void)
{
    PyObject *m = Py_InitModule3("_PyPolyChord", module_methods, module_docstring);
    if (m == NULL)
        return;
}


PyObject* makelist(double* array, int n) {
    PyObject *list = PyList_New(n);
    for (int i=0; i<n; i++) 
        PyList_SET_ITEM(list, i, PyFloat_FromDouble(array[i]));
    return list;
}

/* Callback to the likelihood */
static PyObject *python_loglikelihood = NULL;

static PyObject *set_loglikelihood_PyPolyChord(PyObject *self, PyObject *args)
{
    PyObject *result = NULL;
    PyObject *temp;

    if (PyArg_ParseTuple(args, "O", &temp)) {
        if (!PyCallable_Check(temp)) {
            PyErr_SetString(PyExc_TypeError, "parameter must be callable");
            return NULL;
        }
        Py_XINCREF(temp);                  /* Add a reference to new callback */
        Py_XDECREF(python_loglikelihood);  /* Dispose of previous callback */
        python_loglikelihood = temp;       /* Remember new callback */

        /* Boilerplate to return "None" */
        Py_INCREF(Py_None);
        result = Py_None;
    }
    return result;
}


double loglikelihood(double* theta, int nDims, double* phi, int nDerived)
{
    PyObject* list_theta = makelist(theta,nDims);
    PyObject* list_phi = makelist(phi,nDerived);
    PyObject* pylogL = PyObject_CallFunctionObjArgs(python_loglikelihood,list_theta,list_phi,NULL);
    Py_DECREF(list_theta);
    Py_DECREF(list_phi);
    double logL = PyLong_AsDouble(pylogL);
    Py_DECREF(pylogL);
    return logL;
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
    if (!PyArg_ParseTuple(args, "iiiiiididiiiiiiiiiiss", &nDims, &nDerived, &nlive, &num_repeats, &do_clustering, &feedback, &precision_criterion, &max_ndead, &boost_posterior, &posteriors, &equals, &cluster_posteriors, &write_resume, &write_paramnames, &read_resume, &write_stats, &write_live, &write_dead, &update_files, &base_dir, &file_root))
        return NULL;


    polychord_c_interface( 
            loglikelihood, 
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

