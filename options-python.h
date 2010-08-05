#ifndef __OPTIONS_PYTHON_H__
#define __OPTIONS_PYTHON_H__

#include <Python.h>

// Macros for Python/C interface

#define getLongAttr(obj,name) PyInt_AS_LONG(PyObject_GetAttrString(obj, "name"))
#define getDoubleAttr(obj,name) PyFloat_AS_DOUBLE(PyObject_GetAttrString(obj, "name"))
#define getStringAttr(obj,name) PyString_AS_STRING(PyObject_GetAttrString(obj, "name"))
#define callFunc_IntToNone(obj,name,arg) PyObject_CallObject(PyObject_GetAttrString(obj, "name"), Py_BuildValue("(i)", arg)
#define callFunc_IntToNone(obj,name,arg) PyObject_CallObject(PyObject_GetAttrString(obj, "name"), Py_BuildValue("(i)", arg)


#endif
