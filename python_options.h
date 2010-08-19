#ifndef __PYTHON_OPTIONS_H__
#define __PYTHON_OPTIONS_H__

#include <python2.6/Python.h>

// Macros for Python/C interface

// Import/instantiate
#define newObject(mod, name) PyObject_CallObject(PyObject_GetAttrString(PyImport_ImportModule(#mod), #name), NULL)

// Getters
//#define getLongAttr(obj, name) PyInt_AS_LONG(PyObject_GetAttrString(obj, #name))
#define getBoolAttr(obj, name, pvar) _m_getAttr_DECREF( obj, #name, PyInt_AS_LONG, pvar)
#define getLongAttr(obj, name, pvar) _m_getAttr_DECREF( obj, #name, PyInt_AS_LONG, pvar)
#define getDoubleAttr(obj, name) _m_getAttr_DECREF( obj, #name, PyFloat_AS_DOUBLE, pvar)
#define getStringAttr(obj, name) PyString_AS_STRING(PyObject_GetAttrString(obj, #name))
#define getListAttr(obj, name) PyObject_GetAttrString(obj, #name)

// Setters
#define setDoubleAttr(obj, name, arg) PyObject_SetAttrString(obj, #name, PyFloat_FromDouble(arg))

// Testers
#define testLongAttr(obj, name, test, value) _testLongAttr( obj, #name, #test, value )
#define testBoolAttr(obj, name) _testLongAttr( obj, #name, "=", 1 )

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


// New macros
#define _m_getAttr_DECREF( obj, name, function, pvar )               \
 {                                                                   \
       PyObject *_m_attr = PyObject_GetAttrString( obj, name);       \
       *(pvar) = function(_m_attr);                                  \
       Py_DECREF(_m_attr);                                           \
 }

// Functions
class identlist *makeID_list(PyObject *strand_list);
class stopcomplexes *getStopComplexList(PyObject *options, int index);
class identlist *getID_list(PyObject *options, int index);

static inline bool _testLongAttr( PyObject *obj, const char *attrname, const char *test, long value )
{
 PyObject *_m_attr = PyObject_GetAttrString( obj, attrname);                
 long local_val = PyInt_AS_LONG(_m_attr);
 Py_DECREF(_m_attr);                                                                                        
 if( test[0] == '=' )
       return (local_val == value);
 if( test[0] == '<' )
       return (local_val < value );
 if( test[0] == '>' )
       return (local_val > value );
}












#endif
