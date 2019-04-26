#pragma once
#include "_python.hpp"
#include "_pypolychord.hpp"
#include <vector>

/* Convert from Python list to C array */
inline std::vector<double> list_Py2C_double(PyObject* list) {
    std::vector<double> array;
    int i;
    for (i=0; i<PyList_Size(list); i++) 
    {
        PyObject *obj = PyList_GET_ITEM(list, i);
        if (obj==NULL) throw PythonException();
        if(!PyNumber_Check(obj)) throw PythonException();
        array.push_back(PyFloat_AsDouble(obj));
    }
    return array;
}
inline std::vector<int> list_Py2C_int(PyObject* list) {
    std::vector<int> array;
    int i;
    for (i=0; i<PyList_Size(list); i++) 
    {
        PyObject *obj = PyList_GET_ITEM(list, i);
        if (obj==NULL) throw PythonException();
#ifdef PYTHON3
        if(!PyLong_Check(obj)) throw PythonException();
        array.push_back(PyLong_AsLong(obj));
#else
        if(PyInt_Check(obj)) array.push_back(PyInt_AsLong(obj)); 
        else if(PyLong_Check(obj)) array.push_back(PyLong_AsLong(obj)); 
        else throw PythonException();
#endif
    }
    return array;
}

inline void dict_Py2C(PyObject* dict, std::vector<double>& loglikes, std::vector<int>& nlives) {
    PyObject *key, *value;
    Py_ssize_t pos = 0;
    nlives = {};
    loglikes = {};
    if (dict==NULL) throw PythonException();

    while (PyDict_Next(dict, &pos, &key, &value)) {
        if (key==NULL) throw PythonException();
        if (value==NULL) throw PythonException();
#ifdef PYTHON3
        if(PyLong_Check(value)) nlives.push_back(PyLong_AsLong(value));
        else throw PythonException();
#else
        if(PyInt_Check(value)) nlives.push_back(PyInt_AsLong(value)); 
        else if(PyLong_Check(value)) nlives.push_back(PyInt_AsLong(value)); 
        else throw PythonException();
#endif
        if(PyFloat_Check(key)) loglikes.push_back(PyFloat_AsDouble(key));
        else throw PythonException();
    }
}
