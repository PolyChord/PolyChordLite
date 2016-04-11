#include <Python.h>
#include <stdbool.h>

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

    /* Compute the likelihood and phi from theta  */
    PyObject* answer = PyObject_CallFunctionObjArgs(python_loglikelihood,list_theta,NULL);

    /* Convert the python answer back to a C double */
    double logL = PyFloat_AsDouble( PyTuple_GetItem(answer,0));
    list_Py2C(PyTuple_GetItem(answer,1),phi);

    /* Garbage collect */
    Py_DECREF(list_theta);
    Py_DECREF(answer);

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
    /* Parse the input tuple */
    int i=0;
    python_loglikelihood = PyTuple_GetItem(args,i++);
    python_prior         = PyTuple_GetItem(args,i++);

    int nDims                  = (int)    PyInt_AsLong(     PyTuple_GetItem(args,i++));
    int nDerived               = (int)    PyInt_AsLong(     PyTuple_GetItem(args,i++)); 
    int nlive                  = (int)    PyInt_AsLong(     PyTuple_GetItem(args,i++));  
    int num_repeats            = (int)    PyInt_AsLong(     PyTuple_GetItem(args,i++));
    bool do_clustering         = (bool)   PyObject_IsTrue(  PyTuple_GetItem(args,i++));
    int feedback               = (int)    PyInt_AsLong(     PyTuple_GetItem(args,i++));
    double precision_criterion = (double) PyFloat_AsDouble( PyTuple_GetItem(args,i++));
    int max_ndead              = (int)    PyInt_AsLong(     PyTuple_GetItem(args,i++));
    double boost_posterior     = (double) PyFloat_AsDouble( PyTuple_GetItem(args,i++));
    bool posteriors            = (bool)   PyObject_IsTrue(  PyTuple_GetItem(args,i++));
    bool equals                = (bool)   PyObject_IsTrue(  PyTuple_GetItem(args,i++));
    bool cluster_posteriors    = (bool)   PyObject_IsTrue(  PyTuple_GetItem(args,i++));
    bool write_resume          = (bool)   PyObject_IsTrue(  PyTuple_GetItem(args,i++));
    bool write_paramnames      = (bool)   PyObject_IsTrue(  PyTuple_GetItem(args,i++));
    bool read_resume           = (bool)   PyObject_IsTrue(  PyTuple_GetItem(args,i++));
    bool write_stats           = (bool)   PyObject_IsTrue(  PyTuple_GetItem(args,i++));
    bool write_live            = (bool)   PyObject_IsTrue(  PyTuple_GetItem(args,i++));
    bool write_dead            = (bool)   PyObject_IsTrue(  PyTuple_GetItem(args,i++));
    int update_files           = (int)    PyInt_AsLong(     PyTuple_GetItem(args,i++)); 
    char* base_dir             =          PyString_AsString(PyTuple_GetItem(args,i++));
    char* file_root            =          PyString_AsString(PyTuple_GetItem(args,i++));


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

