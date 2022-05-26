#include "_python.hpp"
#include "_pypolychord.hpp"
#include "_array.hpp"
#include "interfaces.hpp"
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>


/* Initialize the module */
#ifdef PYTHON3
PyMODINIT_FUNC PyInit__pypolychord(void)
{
    import_array();
    return PyModule_Create(&_pypolychordmodule);
}
#else
PyMODINIT_FUNC init_pypolychord(void)
{
    Py_InitModule3("_pypolychord", module_methods, module_docstring);
    import_array();
}
#endif



/* Callback to the likelihood */
static PyObject *python_loglikelihood = NULL;

double loglikelihood(double* theta, int nDims, double* phi, int nDerived)
{
    /* Create a python version of theta */
    npy_intp theta_shape[] = {nDims};            
    PyObject *array_theta = PyArray_SimpleNewFromData(1, theta_shape, NPY_DOUBLE, theta);
    if (array_theta==NULL) throw PythonException();
    PyArray_CLEARFLAGS(reinterpret_cast<PyArrayObject*>(array_theta), NPY_ARRAY_WRITEABLE);

    /* Create a python version of phi */
    npy_intp phi_shape[] = {nDerived};            
    PyObject *array_phi = PyArray_SimpleNewFromData(1, phi_shape, NPY_DOUBLE, phi);
    if (array_phi==NULL) {Py_DECREF(array_theta); throw PythonException();}

    /* Compute the likelihood and phi from theta  */
    PyObject* py_logL = PyObject_CallFunctionObjArgs(python_loglikelihood,array_theta,array_phi,NULL);
    if (py_logL==NULL){ Py_DECREF(array_theta); Py_DECREF(array_phi); throw PythonException(); }

    /* Check return is a float  */
    if (!PyFloat_Check(py_logL)){
        Py_DECREF(array_theta); Py_DECREF(array_phi); Py_DECREF(py_logL);
        PyErr_SetString(PyExc_TypeError,"loglikelihood must be a float (element 0 of loglikelihood return)");
        throw PythonException();
    }
    double logL = PyFloat_AsDouble(py_logL);

    /* Garbage collect */
    Py_DECREF(array_theta); Py_DECREF(array_phi); Py_DECREF(py_logL);

    return logL;
}

/* Callback to the prior */
static PyObject *python_prior = NULL;

void prior(double* cube, double* theta, int nDims)
{
    /* create a python version of cube */
    npy_intp shape[] = {nDims};            
    PyObject *array_cube = PyArray_SimpleNewFromData(1, shape, NPY_DOUBLE, cube);
    if (array_cube==NULL) throw PythonException();
    PyArray_CLEARFLAGS(reinterpret_cast<PyArrayObject*>(array_cube), NPY_ARRAY_WRITEABLE);

    /* create a python version of theta */
    PyObject *array_theta = PyArray_SimpleNewFromData(1, shape, NPY_DOUBLE, theta);
    if (array_theta==NULL) {Py_DECREF(array_cube); throw PythonException();}

    /* Compute theta from the prior */
    PyObject_CallFunctionObjArgs(python_prior,array_cube,array_theta,NULL);

    /* Garbage collect */
    Py_DECREF(array_theta); Py_DECREF(array_cube);
}

/* Callback to the dumper */
static PyObject *python_dumper = NULL;

void dumper(int ndead, int nlive, int npars, double* live, double* dead, double* logweights, double logZ, double logZerr)
{
    /* create a python version of live points */
    npy_intp live_shape[] = {nlive, npars};            
    PyObject *array_live = PyArray_SimpleNewFromData(2, live_shape, NPY_DOUBLE, live);
    if (array_live==NULL) throw PythonException();
    PyArray_CLEARFLAGS(reinterpret_cast<PyArrayObject*>(array_live), NPY_ARRAY_WRITEABLE);

    /* create a python version of dead points */
    npy_intp dead_shape[] = {ndead, npars};            
    PyObject *array_dead = PyArray_SimpleNewFromData(2, dead_shape, NPY_DOUBLE, dead);
    if (array_dead==NULL) {Py_DECREF(array_live); throw PythonException();}
    PyArray_CLEARFLAGS(reinterpret_cast<PyArrayObject*>(array_dead), NPY_ARRAY_WRITEABLE);

    /* create a python version of posterior weights */
    npy_intp logweights_shape[] = {ndead};            
    PyObject *array_logweights = PyArray_SimpleNewFromData(1, logweights_shape, NPY_DOUBLE, logweights);
    if (array_logweights==NULL) {Py_DECREF(array_live); Py_DECREF(array_dead); throw PythonException();}
    PyArray_CLEARFLAGS(reinterpret_cast<PyArrayObject*>(array_logweights), NPY_ARRAY_WRITEABLE);

    PyObject *py_logZ = Py_BuildValue("d",logZ);
    if (py_logZ==NULL) {Py_DECREF(array_live); Py_DECREF(array_dead); Py_DECREF(array_logweights); throw PythonException();}
    PyObject *py_logZerr = Py_BuildValue("d",logZerr);
    if (py_logZerr==NULL) {Py_DECREF(array_live); Py_DECREF(array_dead); Py_DECREF(array_logweights); Py_DECREF(py_logZ); throw PythonException();}

    /* Call dumper */
    PyObject_CallFunctionObjArgs(python_dumper,array_live,array_dead,array_logweights,py_logZ,py_logZerr,NULL);

    /* Garbage collect */
    Py_DECREF(array_live); Py_DECREF(array_dead); Py_DECREF(array_logweights); Py_DECREF(py_logZ); Py_DECREF(py_logZerr);
}

/* Callback to the cluster */
static PyObject *python_cluster = NULL;

void cluster(double* points, int* cluster_list, int nDims, int nPoints)
{
    /* create a python version of points */
    npy_intp shape[] = {nPoints,nDims};            
    PyObject *array_points = PyArray_SimpleNewFromData(2, shape, NPY_DOUBLE, points);
    if (array_points ==NULL) throw PythonException();
    PyArray_CLEARFLAGS(reinterpret_cast<PyArrayObject*>(array_points), NPY_ARRAY_WRITEABLE);

    /* create a python version of cluster_list */
    npy_intp shape1[] = {nPoints};            
    PyObject *array_cluster_list = PyArray_SimpleNewFromData(1, shape1, NPY_INT, cluster_list);
    if (array_cluster_list==NULL) {Py_DECREF(array_points); throw PythonException();}

    /* Compute cluster_list from the cluster */
    PyObject_CallFunctionObjArgs(python_cluster,array_points,array_cluster_list,NULL);

    /* Garbage collect */
    Py_DECREF(array_cluster_list); Py_DECREF(array_points);
}


/* Function to run pypolychord */
static PyObject *run_pypolychord(PyObject *, PyObject *args)
{
    Settings S;

    PyObject *temp_logl, *temp_prior, *temp_dumper, *temp_cluster;
    PyObject* py_grade_dims, *py_grade_frac, *py_nlives;
    char* base_dir, *file_root;
        

    if (!PyArg_ParseTuple(args,
                "OOOOiiiiiiiiddidiiiiiiiiiiidissO!O!O!i:run",
                &temp_logl,
                &temp_prior,
                &temp_dumper,
                &temp_cluster,
                &S.nDims,
                &S.nDerived,
                &S.nlive,
                &S.num_repeats,
                &S.nprior,
                &S.nfail,
                &S.do_clustering,
                &S.feedback,
                &S.precision_criterion,
                &S.logzero,
                &S.max_ndead,
                &S.boost_posterior,
                &S.posteriors,
                &S.equals,
                &S.cluster_posteriors,
                &S.write_resume,
                &S.write_paramnames,
                &S.read_resume,
                &S.write_stats,
                &S.write_live,
                &S.write_dead,
                &S.write_prior,
                &S.maximise,
                &S.compression_factor,
                &S.synchronous,
                &base_dir,
                &file_root,
                &PyList_Type,
                &py_grade_frac,
                &PyList_Type,
                &py_grade_dims,
                &PyDict_Type,
                &py_nlives,
                &S.seed
                )
            )
        return NULL;
    S.base_dir = base_dir;
    S.file_root = file_root;
    
    int nGrade = PyList_Size(py_grade_frac);
    S.grade_frac.resize(nGrade);
    S.grade_dims.resize(nGrade);
    Py_INCREF(py_grade_frac); Py_INCREF(py_grade_dims);

    try{ S.grade_frac = list_Py2C_double(py_grade_frac); }
    catch (PythonException& e){
        Py_DECREF(py_grade_frac); Py_DECREF(py_grade_dims);
        PyErr_SetString(PyExc_TypeError,"grade_frac must be a list of doubles");
        return NULL;
    }
    try{ S.grade_dims = list_Py2C_int(py_grade_dims); }
    catch (PythonException& e){
        Py_DECREF(py_grade_frac); Py_DECREF(py_grade_dims);
        PyErr_SetString(PyExc_TypeError,"grade_dims must be a list of integers");
        return NULL;
    }
    if (PyList_Size(py_grade_frac) != PyList_Size(py_grade_dims)){
        Py_DECREF(py_grade_frac); Py_DECREF(py_grade_dims);
        PyErr_SetString(PyExc_ValueError,"grade_dims and grade_frac must have the same size");
        return NULL;
    }
    int tot = 0; for (int i=0;i<nGrade;i++) tot += S.grade_dims[i];
    if (tot != S.nDims) {
        PyErr_SetString(PyExc_ValueError,"grade_dims must sum to nDims");
        return NULL;
    }
    try{ dict_Py2C(py_nlives,S.loglikes,S.nlives); }
    catch (PythonException& e){
        PyErr_SetString(PyExc_TypeError,"nlives must be a dict mapping floats to integers");
        return NULL;
    }

    Py_XINCREF(temp_logl);
    Py_XDECREF(python_loglikelihood);
    python_loglikelihood = temp_logl;

    Py_XINCREF(temp_prior);
    Py_XDECREF(python_prior);
    python_prior = temp_prior;

    Py_XINCREF(temp_dumper);
    Py_XDECREF(python_dumper);
    python_dumper = temp_dumper;

    Py_XINCREF(temp_cluster);
    Py_XDECREF(python_cluster);
    python_cluster = temp_cluster;

    /* Run PolyChord */
    try{ run_polychord(loglikelihood, prior, dumper, cluster, S); }
    catch (PythonException& e)
    { 
        Py_DECREF(py_grade_frac);Py_DECREF(py_grade_dims);Py_DECREF(python_loglikelihood);Py_DECREF(python_prior);Py_DECREF(python_cluster); 
        return NULL; 
    }

    /* Return None */
    Py_RETURN_NONE;
}

