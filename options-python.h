#ifndef __OPTIONS_PYTHON_H__
#define __OPTIONS_PYTHON_H__

#include <python2.6/Python.h>


// Macros for Python/C interface

// Import/instantiate
#define newObject(mod, name) PyObject_CallObject(PyObject_GetAttrString(PyImport_ImportModule(#mod), #name), NULL)

// Getters
#define getBoolAttr(obj, name) PyInt_AS_LONG(PyObject_GetAttrString(obj, #name))
#define getLongAttr(obj, name) PyInt_AS_LONG(PyObject_GetAttrString(obj, #name))
#define getDoubleAttr(obj, name) PyFloat_AS_DOUBLE(PyObject_GetAttrString(obj, #name))
#define getStringAttr(obj, name) PyString_AS_STRING(PyObject_GetAttrString(obj, #name))
#define getListAttr(obj, name) PyObject_GetAttrString(obj, #name)

// Setters
#define setDoubleAttr(obj, name, arg) PyObject_SetAttrString(obj, #name, PyFloat_FromDouble(arg))

// Function calls
#define callFunc_NoArgsToNone(obj, name) PyObject_CallMethod(obj, #name, "()")
#define callFunc_DoubleToNone(obj, name, arg) PyObject_CallMethod(obj, #name, "(f)", arg)

// List indexing
#define getStringItem(list, index) PyString_AS_STRING(PyList_GET_ITEM(list, index))


// Not currently used, but might be a good reference for later
#define callFunc_IntToNone(obj, name, arg) PyObject_CallObject(PyObject_GetAttrString(obj, #name), Py_BuildValue("(i)", arg))
#define callFunc_IntToString(obj, name, arg) PyString_AS_STRING(PyObject_CallObject(PyObject_GetAttrString(obj, #name), Py_BuildValue("(i)", arg)))

/*****************************************************

 #defines for const values used in python_options.py

******************************************************/

/* WARNING: If you change the following defines, you must also
            change the values in python_options._OptionsConstants.RATEMETHOD
			in the file python_options.py.
*/
#define RATE_METHOD_INVALID         0x00
#define RATE_METHOD_METROPOLIS      0x01
#define RATE_METHOD_KAWASAKI        0x02
#define RATE_METHOD_ENTROPYENTHALPY 0x03


/* WARNING: If you change the following defines, you must also
            change the values in python_options._OptionsConstants.DANGLES
			in the file python_options.py.
*/
#define DANGLES_NONE    0x00
#define DANGLES_SOME    0x01
#define DANGLES_ALL     0x02

/* WARNING: If you change the following defines, you must also
            change the values in python_options._OptionsConstants.ENERGYMODEL_TYPE
			in the file python_options.py.
*/
#define ENERGYMODEL_VIENNA 0x00
#define ENERGYMODEL_NUPACK 0x01

/* WARNING: If you change the following defines, you must also
            change the values in python_options._OptionsConstants.SUBSTRATE_TYPE
			in the file python_options.py.
*/
#define SUBSTRATE_INVALID 0x00
#define SUBSTRATE_RNA     0x01
#define SUBSTRATE_DNA     0x02


#endif
