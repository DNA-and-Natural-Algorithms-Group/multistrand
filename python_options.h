#ifndef __PYTHON_OPTIONS_H__
#define __PYTHON_OPTIONS_H__

#include <python2.6/Python.h>


// Macros for Python/C interface

// Import/instantiate
#define newObject(mod, name) PyObject_CallObject(PyObject_GetAttrString(PyImport_ImportModule(#mod), #name), NULL)

// Getters
#define getLongAttr(obj, name) PyInt_AS_LONG(PyObject_GetAttrString(obj, #name))
#define getDoubleAttr(obj, name) PyFloat_AS_DOUBLE(PyObject_GetAttrString(obj, #name))
#define getStringAttr(obj, name) PyString_AS_STRING(PyObject_GetAttrString(obj, #name))
#define getListAttr(obj, name) PyObject_GetAttrString(obj, #name)
#define getStopcomplexList(obj) convert_stopcomplex_list(PyObject_GetAttrString(obj, "stopcomplexes"))

// Setters
#define setDoubleAttr(obj, name, arg) PyObject_SetAttrString(obj, #name, PyFloat_FromDouble(arg))

// Function calls
#define callFunc_NoArgsToNone(obj, name) PyObject_CallMethod(obj, #name, "()")
#define callFunc_DoubleToNone(obj, name, arg) PyObject_CallMethod(obj, #name, "(f)", arg)

// List indexing
#define getStringItem(list, index) PyString_AS_STRING(PyList_GET_ITEM(list, index))
#define getLongItem(list, index) PyInt_AS_LONG(PyList_GET_ITEM(list, index))
#define getLongItemFromTuple(tuple, index) PyInt_AS_LONG(PyTuple_GET_ITEM(tuple, index))


// Not currently used, but might be a good reference for later
#define callFunc_IntToNone(obj, name, arg) PyObject_CallObject(PyObject_GetAttrString(obj, #name), Py_BuildValue("(i)", arg))
#define callFunc_IntToString(obj, name, arg) PyString_AS_STRING(PyObject_CallObject(PyObject_GetAttrString(obj, #name), Py_BuildValue("(i)", arg)))


// Conversion functions
class stopcomplexes *convert_stopcomplex_list(*PyObject stopcomplexes_python);
















#endif
