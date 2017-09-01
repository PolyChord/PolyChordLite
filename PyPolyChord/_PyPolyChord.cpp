#include <Python.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdexcept>

#if PY_MAJOR_VERSION >=3
#define PYTHON3
#endif

extern void polychord_c_interface(double (*)(double*,int,double*,int), void (*)(double*,double*,int), int, int, bool, int, double, int, double, bool, bool, bool, bool, bool, bool, bool, bool, bool, int, int, int, char*, char*, int, double*, int* ); 

/* Exception */
class PythonException : std::exception {};


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

#ifdef PYTHON3
static struct PyModuleDef cModPyDem =
{
    PyModuleDef_HEAD_INIT,
    "_PyPolyChord", /* name of module */
    module_docstring,          /* module documentation, may be NULL */
    -1,          /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    module_methods
};
#endif

/* Initialize the module */
#ifdef PYTHON3
PyMODINIT_FUNC PyInit__PyPolyChord(void)
{
    return PyModule_Create(&cModPyDem);
}
#else
PyMODINIT_FUNC init_PyPolyChord(void)
{
    Py_InitModule3("_PyPolyChord", module_methods, module_docstring);
}
#endif


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
    if (answer==NULL){
        Py_DECREF(list_theta);
        throw PythonException();
    }

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

    int nDims                  = (int)    PyLong_AsLong(     PyTuple_GetItem(args,i++));
    int nDerived               = (int)    PyLong_AsLong(     PyTuple_GetItem(args,i++)); 
    int nlive                  = (int)    PyLong_AsLong(     PyTuple_GetItem(args,i++));  
    int num_repeats            = (int)    PyLong_AsLong(     PyTuple_GetItem(args,i++));
    bool do_clustering         = (bool)   PyObject_IsTrue(  PyTuple_GetItem(args,i++));
    int feedback               = (int)    PyLong_AsLong(     PyTuple_GetItem(args,i++));
    double precision_criterion = (double) PyFloat_AsDouble( PyTuple_GetItem(args,i++));
    int max_ndead              = (int)    PyLong_AsLong(     PyTuple_GetItem(args,i++));
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
    int update_files           = (int)    PyLong_AsLong(     PyTuple_GetItem(args,i++)); 
#ifdef PYTHON3
    char* base_dir             =          _PyUnicode_AsString(PyTuple_GetItem(args,i++));
    char* file_root            =          _PyUnicode_AsString(PyTuple_GetItem(args,i++));
#else
    char* base_dir             =          PyString_AsString(PyTuple_GetItem(args,i++));
    char* file_root            =          PyString_AsString(PyTuple_GetItem(args,i++));
#endif

    PyObject * py_grade_frac = PyTuple_GetItem(args,i++);
    PyObject * py_grade_dims = PyTuple_GetItem(args,i++);
    int nGrade = (int) PyList_Size(py_grade_frac);
    double grade_frac[nGrade];
    int grade_dims[nGrade];

    int j;
    for(j=0; j<nGrade; j++)
    {
        grade_frac[j] = (double) PyFloat_AsDouble( PyList_GET_ITEM(py_grade_frac,j) );
        grade_dims[j] = (int) PyLong_AsLong( PyList_GET_ITEM(py_grade_dims,j) );
    }




    /* Run PolyChord */
    try{
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
            file_root,
            nGrade,
            grade_frac,
            grade_dims
            );
    }
    catch (PythonException& e)
    {
        Py_DECREF(python_loglikelihood);
        return NULL;
    }

    /* Return None */
    Py_DECREF(python_loglikelihood);
    Py_RETURN_NONE;
}

