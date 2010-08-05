#ifndef __OPTIONS_PYTHON_H__
#define __OPTIONS_PYTHON_H__

#include <python2.6/Python.h>


// Macros for Python/C interface

// Getters
#define getLongAttr(obj, name) PyInt_AS_LONG(PyObject_GetAttrString(obj, "name"))
#define getDoubleAttr(obj, name) PyFloat_AS_DOUBLE(PyObject_GetAttrString(obj, "name"))
#define getStringAttr(obj, name) PyString_AS_STRING(PyObject_GetAttrString(obj, "name"))
#define getListAttr(obj, name) PyObject_GetAttrString(obj, "name")

// Setters
#define setDoubleAttr(obj, name, arg) PyObject_SetAttrString(obj, "name", PyFloat_FromDouble(arg))

// Function calls
#define callFunc_NoArgsToNone(obj, name) PyObject_CallObject(PyObject_GetAttrString(obj, "name"), Py_BuildValue("()"))
#define callFunc_DoubleToNone(obj, name, arg) PyObject_CallObject(PyObject_GetAttrString(obj, "name"), Py_BuildValue("(f)", PyFloat_FromDouble(arg)))

// List indexing
#define getStringItem(list, index) PyString_AS_STRING(PyList_GET_ITEM(list, index))


// Not currently used, but might be a good reference for later
#define callFunc_IntToNone(obj, name, arg) PyObject_CallObject(PyObject_GetAttrString(obj, "name"), Py_BuildValue("(i)", arg))
#define callFunc_IntToString(obj, name, arg) PyString_AS_STRING(PyObject_CallObject(PyObject_GetAttrString(obj, "name"), Py_BuildValue("(i)", arg)))

#endif
