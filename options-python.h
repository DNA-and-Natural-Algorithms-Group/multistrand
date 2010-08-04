#ifndef __OPTIONS_PYTHON_H__
#define __OPTIONS_PYTHON_H__

#include <Python.h>

// Macros for Python/C interface

#define getLongAttr(x,y) PyInt_AS_LONG(PyObject_GetAttrString(x, "y"))
#define getDoubleAttr(x,y) PyFloat_AS_DOUBLE(PyObject_GetAttrString(x, "y"))
#define getStringAttr(x,y) PyString_AS_STRING(PyObject_GetAttrString(x, "y"))


#endif
