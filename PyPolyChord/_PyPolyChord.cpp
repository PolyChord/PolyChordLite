#include <Python.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdexcept>

#if PY_MAJOR_VERSION >=3
#define PYTHON3
#endif

extern "C" void polychord_c_interface(double (*)(double*,int,double*,int), void (*)(double*,double*,int), int, int, bool, int, double, int, double, bool, bool, bool, bool, bool, bool, bool, bool, bool, int, int, int, char*, char*, int, double*, int* ); 

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
PyObject* list_C2Py(double* array, int nDims) {
    PyObject* list = PyList_New(nDims);
    for (int i=0; i<PyList_Size(list); i++) 
        PyList_SET_ITEM(list, i, PyFloat_FromDouble(array[i]));
    return list;
}

/* Convert from Python list to C array */
void list_Py2C(PyObject* list, double* array) {
    int i;
    for (i=0; i<PyList_Size(list); i++) 
    {
        PyObject *obj = PyList_GET_ITEM(list, i);
        if (obj==NULL) throw PythonException();
        if(!PyFloat_Check(obj)) throw PythonException();
        array[i] = PyFloat_AsDouble(obj);
    }
}

/* Callback to the likelihood and prior */
static PyObject *python_loglikelihood = NULL;

double loglikelihood(double* theta, int nDims, double* phi, int nDerived)
{
    /* Create a python version of theta */
    PyObject *list_theta = list_C2Py(theta,nDims);
    if (list_theta==NULL) throw PythonException();

    /* Compute the likelihood and phi from theta  */
    PyObject* answer = PyObject_CallFunctionObjArgs(python_loglikelihood,list_theta,NULL);
    if (answer==NULL){
        Py_DECREF(list_theta);
        throw PythonException();
    }

    if (!PyTuple_Check(answer) or PyTuple_Size(answer) != 2){
        Py_DECREF(list_theta); Py_DECREF(answer);
        PyErr_SetString(PyExc_TypeError,"Return from loglikelihood must be a tuple of (loglikelihood, [derived parameters])");
        throw PythonException();
    }
    
    PyObject *py_logL = PyTuple_GetItem(answer,0);
    if (!PyFloat_Check(py_logL)){
        Py_DECREF(list_theta); Py_DECREF(answer);
        PyErr_SetString(PyExc_TypeError,"loglikelihood must be a float (element 0 of loglikelihood return)");
        throw PythonException();
    }
    double logL = PyFloat_AsDouble( PyTuple_GetItem(answer,0));

    PyObject *py_derived = PyTuple_GetItem(answer,1);
    if (!PyList_Check(py_derived)){
        Py_DECREF(list_theta); Py_DECREF(answer);
        PyErr_SetString(PyExc_TypeError,"Derived parameters must be a list (element 1 of loglikelihood return)");
        throw PythonException();
    } 
    else if(PyList_Size(py_derived) != nDerived)
    {
        Py_DECREF(list_theta); Py_DECREF(answer);
        PyErr_SetString(PyExc_TypeError,"Derived parameters must have length nDerived (element 1 of loglikelihood return)");
        throw PythonException();
    }

    try{ list_Py2C(PyTuple_GetItem(answer,1),phi); }
    catch (PythonException& e) {
        PyErr_SetString(PyExc_TypeError,"Derived parameters must be a list of floats (element 1 of loglikelihood return)");
        throw PythonException(); }

    /* Garbage collect */
    Py_DECREF(list_theta);
    Py_DECREF(answer);

    return logL;
}

static PyObject *python_prior = NULL;

void prior(double* cube, double* theta, int nDims)
{
    /* create a python version of cube */
    PyObject* list_cube = list_C2Py(cube,nDims);
    if (list_cube==NULL) throw PythonException();

    /* Compute theta from the prior */
    PyObject* list_theta = PyObject_CallFunctionObjArgs(python_prior,list_cube,NULL);
    if (list_theta==NULL) {Py_DECREF(list_cube); throw PythonException();}

    if (!PyList_Check(list_theta)){
        Py_DECREF(list_theta); Py_DECREF(list_cube);
        PyErr_SetString(PyExc_TypeError,"Physical parameters must be a list (return from prior)");
        throw PythonException();
    } 
    else if(PyList_Size(list_theta) != nDims)
    {
        Py_DECREF(list_theta); Py_DECREF(list_cube);
        PyErr_SetString(PyExc_TypeError,"Physical parameters must have length nDims (return from prior)");
        throw PythonException();
    }

    /* Convert the python answer back to a C array */
    try{ list_Py2C(list_theta,theta); }
    catch (PythonException& e){
        Py_DECREF(list_theta); Py_DECREF(list_cube);
        PyErr_SetString(PyExc_TypeError,"Physical parameters must be a list of floats (return from prior)");
        throw PythonException(); 
    }

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
    Py_RETURN_NONE;
}

