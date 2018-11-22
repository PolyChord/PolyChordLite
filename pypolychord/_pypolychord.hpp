#pragma once
#include "_python.hpp"
#include <stdexcept>

/* Exception */
class PythonException : std::exception {};


/* Docstrings */
static char module_docstring[] =
    "pypolychord: This module provides a Python interface to PolyChord.";
static char run_docstring[] =
    "Runs pypolychord";


/* Available functions */
static PyObject *run_pypolychord(PyObject *self, PyObject *args);


/* Module interface */
static PyMethodDef module_methods[] = {
    {"run", run_pypolychord, METH_VARARGS, run_docstring},
    {NULL, NULL, 0, NULL}
};

#ifdef PYTHON3
static struct PyModuleDef _pypolychordmodule =
{
    PyModuleDef_HEAD_INIT,
    "_pypolychord", /* name of module */
    module_docstring,          /* module documentation, may be NULL */
    -1,          /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    module_methods
};
#endif

// static PyObject *
// moduleinit(void)
// {
//     PyObject *m;

// #ifdef PYTHON3
//     m = PyModule_Create(&_pypolychordmodule);
// #else
//     m = Py_InitModule3("extension",
//                         module_methods, "This is a module");
// #endif
//   return m;
// }

// #ifdef PYTHON3
//     PyMODINIT_FUNC
//     initextension(void)
//     {
//         moduleinit();
//     }
// #else
//     PyMODINIT_FUNC
//     PyInit_extension(void)
//     {
//         moduleinit();
//     }
// #endif