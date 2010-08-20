#ifndef __PYTHON_OPTIONS_H__
#define __PYTHON_OPTIONS_H__

#include <python2.6/Python.h>

// Macros for Python/C interface

// Import/instantiate
#define newObject(mod, name) PyObject_CallObject(PyObject_GetAttrString(PyImport_ImportModule(#mod), #name), NULL)

#define _m_getAttr_DECREF( obj, name, function, pvar )				\
  {																	\
	PyObject *_m_attr = PyObject_GetAttrString( obj, #name);		\
	*(pvar) = function(_m_attr);									\
	Py_DECREF(_m_attr);												\
  }

bool _testLongAttr( PyObject *obj, const char *attrname, const char *test, long value )
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

#define testLongAttr(obj, name, test, value) \
  _testLongAttr( obj, #name, #test, value )

#define getBoolAttr(obj, name)  _getBoolAttr( #obj, #name)

#define _getBoolAttr(obj, name) \
  (PyRun_SimpleString('"' obj "." name "== True" '"')==0?0:1)

// Getters
//#define getLongAttr(obj, name) PyInt_AS_LONG(PyObject_GetAttrString(obj, #name))
#define getBoolAttr(obj, name, pvar) _m_getAttr_DECREF( obj, #name, PyInt_AS_LONG, pvar)
// TODO: fix 
#define getLongAttr(obj, name, pvar) _m_getAttr_DECREF( obj, #name, PyInt_AS_LONG, pvar)
#define getDoubleAttr(obj, name, pvar) _m_getAttr_DECREF( obj, #name, PyFloat_AS_DOUBLE, pvar)
#define getStringAttr(obj, name, pyo) PyString_AS_STRING(pyo=PyObject_GetAttrString(obj, #name))
#define getListAttr(obj, name) PyObject_GetAttrString(obj, #name)

// Setters
#define setDoubleAttr(obj, name, arg) PyObject_SetAttrString(obj, #name, PyFloat_FromDouble(arg))

// Testers
#define testLongAttr(obj, name, test, value) _testLongAttr( obj, #name, #test, value )
#define testBoolAttr(obj, name) _testLongAttr( obj, #name, "=", 1 )

// Function calls
#define callFunc_NoArgsToNone(obj, name) PyObject_CallMethod(obj, #name, "()")
#define callFunc_NoArgsToLong(obj, name) PyInt_AS_LONG(PyObject_CallMethod(obj, #name, "()"))
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

/* WARNING: If you change the following defines, you must also
            change the values in python_options._OptionsConstants.SIMULATION_MODE
			in the file python_options.py.
*/

#define SIMULATION_MODE_NORMAL              0x00
#define SIMULATION_MODE_FIRST_BIMOLECULAR   0x01
#define SIMULATION_MODE_PYTHON_NORMAL       0x02
#define SIMULATION_MODE_PYTHON_FIRST_BI     0x03

// simulation modes are bitwise -> bit 0 is normal/first bi
//                                 bit 1 is normal interface/python interface
// the following are the bit definitions for tests on those:

#define SIMULATION_MODE_FLAG_FIRST_BIMOLECULAR         0x01
#define SIMULATION_MODE_FLAG_PYTHON                    0x02


#endif
