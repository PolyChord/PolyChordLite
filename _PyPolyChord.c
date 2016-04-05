#include <Python.h>
#include <stdbool.h>
#include <stdio.h>

extern void polychord_c_interface(double (*)(double*,int,double*,int), void (*)(double*,double*,int), int, int, bool, int, double, int, double, bool, bool, bool, bool, bool, bool, bool, bool, bool, int, int, int, char*, char*); 

/* Docstrings */
static char module_docstring[] =
    "PyPolyChord: This module provides a Python interface to PolyChord.";
static char run_docstring[] =
    "Runs PyPolyChord";


/* Available functions */
static PyObject *run_PyPolyChord(PyObject *self, PyObject *args);


/* Module interface */
static PyMethodDef module_methods[] = {
    {"run", run_PyPolyChord, METH_VARARGS, run_docstring},
    {NULL, NULL, 0, NULL}
};


/* Initialize the module */
PyMODINIT_FUNC init_PyPolyChord(void)
{
    Py_InitModule3("_PyPolyChord", module_methods, module_docstring);
}


/* Convert from C array to Python list */
void list_C2Py(double* array, PyObject* list) {
    int i;
    for (i=0; i<PyList_Size(list); i++) 
        PyList_SET_ITEM(list, i, PyFloat_FromDouble(array[i]));
}

/* Convert from Python list to C array */
void list_Py2C(PyObject* list, double* array) {
    int i;
    for (i=0; i<PyList_Size(list); i++) 
        array[i] = PyFloat_AsDouble(PyList_GET_ITEM(list, i));
}

/* Callback to the likelihood and prior */
static PyObject *python_loglikelihood = NULL;

double loglikelihood(double* theta, int nDims, double* phi, int nDerived)
{
    /* Create a python version of theta */
    PyObject* list_theta = PyList_New(nDims);
    list_C2Py(theta,list_theta);

    /* Create a python version of phi */
    PyObject* list_phi = PyList_New(nDerived);
    list_C2Py(phi,list_phi);

    /* Compute the likelihood and phi from theta  */
    PyObject* pylogL = PyObject_CallFunctionObjArgs(python_loglikelihood,list_theta,list_phi,NULL);

    /* Convert the python answer back to a C double */
    double logL = PyFloat_AsDouble(pylogL);

    /* Garbage collect */
    Py_DECREF(list_theta);
    Py_DECREF(list_phi);
    Py_DECREF(pylogL);

    return logL;
}

static PyObject *python_prior = NULL;

void prior(double* cube, double* theta, int nDims)
{
    /* create a python version of cube */
    PyObject* list_cube = PyList_New(nDims);
    list_C2Py(cube,list_cube);

    /* Compute theta from the prior */
    PyObject* list_theta = PyObject_CallFunctionObjArgs(python_prior,list_cube,NULL);

    /* Convert the python answer back to a C array */
    list_Py2C(list_theta,theta);

    /* Garbage collect */
    Py_DECREF(list_cube);
    Py_DECREF(list_theta);
}

/* Function to run PyPolyChord */
static PyObject *run_PyPolyChord(PyObject *self, PyObject *args)
{
    /* Inputs to PolyChord in the order that they are passed to python */
    int nDims;
    int nDerived;
    int nlive;
    int num_repeats;
    bool do_clustering;
    bool feedback;
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
    char* base_dir;
    char* file_root;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "OOiiiiiididiiiiiiiiiiss", &python_loglikelihood, &python_prior, &nDims, &nDerived, &nlive, &num_repeats, &do_clustering, &feedback, &precision_criterion, &max_ndead, &boost_posterior, &posteriors, &equals, &cluster_posteriors, &write_resume, &write_paramnames, &read_resume, &write_stats, &write_live, &write_dead, &update_files, &base_dir, &file_root))
        return NULL;

    /* Run PolyChord */
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
            file_root
            );

    /* Return None */
    PyObject *ret = Py_None;
    Py_INCREF(Py_None);
    return ret;
}

