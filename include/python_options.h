#ifndef __PYTHON_OPTIONS_H__
#define __PYTHON_OPTIONS_H__

#include <python2.6/Python.h>

// Macros for Python/C interface

/***************************************/
/* Helper functions / internal macros. */
/***************************************/

#ifndef DEBUG_MACROS
#define _m_getAttr_DECREF( obj, name, function, pvar, vartype )     \
  {																	\
	PyObject *_m_attr = PyObject_GetAttrString( obj, name);		\
	*(vartype *)(pvar) = function(_m_attr);                         \
	Py_DECREF(_m_attr);												\
  }

#define _m_setAttr_DECREF( obj, name, function, arg )               \
  {																	\
	PyObject *val = function(arg);                                  \
    PyObject_SetAttrString( obj, #name, val);                       \
	Py_DECREF(val);                                                 \
  }

// Import/instantiate
#define newObject(mod, name) _m_newObject( #mod, #name )

// Getters
#define getBoolAttr(obj, name, pvar) _m_getAttr_DECREF( obj, #name, PyInt_AS_LONG, pvar, bool)
#define getLongAttr(obj, name, pvar) _m_getAttr_DECREF( obj, #name, PyInt_AS_LONG, pvar, long)
#define getDoubleAttr(obj, name, pvar) _m_getAttr_DECREF( obj, #name, PyFloat_AS_DOUBLE, pvar, double)
#define getStringAttr(obj, name, pyo) ((char *)PyString_AS_STRING(pyo=PyObject_GetAttrString(obj, #name)))
#define getListAttr(obj, name) PyObject_GetAttrString(obj, #name)

#define pingAttr(obj, name) Py_XDECREF(PyObject_GetAttrString( obj, #name ))
// note: does not do anything crazy on a NULL return from GetAttrString, but if that 
// returned null it might be an error...

// List indexing
#define getStringItem(list, index) PyString_AS_STRING(PyList_GET_ITEM(list, index))
#define getLongItem(list, index) PyInt_AS_LONG(PyList_GET_ITEM(list, index))
#define getLongItemFromTuple(tuple, index) PyInt_AS_LONG(PyTuple_GET_ITEM(tuple, index))

// Setters
#define setDoubleAttr(obj, name, arg) _m_setAttr_DECREF( obj, #name, PyFloat_FromDouble, arg)
#define setLongAttr(obj, name, arg) _m_setAttr_DECREF( obj, #name, PyLong_FromLong, arg)

// Testers
#define testLongAttr(obj, name, test, value) _m_testLongAttr( obj, #name, #test, value )
#define testBoolAttr(obj, name) _m_testLongAttr( obj, #name, "=", 1 )

// Function calls
#define callFunc_NoArgsToNone(obj, name) PyObject_CallMethod(obj, #name, "()")
#define callFunc_NoArgsToLong(obj, name) PyInt_AS_LONG(PyObject_CallMethod(obj, #name, "()"))
#define callFunc_DoubleToNone(obj, name, arg) PyObject_CallMethod(obj, #name, "(f)", arg)


// Not currently used, but might be a good reference for later
#define callFunc_IntToNone(obj, name, arg) PyObject_CallObject(PyObject_GetAttrString(obj, #name), Py_BuildValue("(i)", arg))
// these are refcounting insensitive at the moment.
#define callFunc_IntToString(obj, name, arg) PyString_AS_STRING(PyObject_CallObject(PyObject_GetAttrString(obj, #name), Py_BuildValue("(i)", arg)))





// Print calls, currently null statements.
#define m_printStatusLine( obj,a,b,c)
//#define m_printStatusLine( obj,a)
#define m_printTrajLine(obj,a,b)
#define m_printStatusLine_First_Bimolecular( obj,a,b,c,d,e)
#define m_printStatusLine_Final_First_Bimolecular( obj, a,b,c,d,e,f )
#define m_printStatusLine_Warning( obj, a, b )

#endif  // DEBUG_MACROS is not set.

/***************************************************

  Debug version of macros:

 **************************************************/
#ifdef DEBUG_MACROS

#define _m_printPyError_withLineNumber() \
  fprintf(stderr,"ERROR: Python Interpreter error: file %s, line %d.\nERROR (Python): ", __FILE__, (int)__LINE__), PyErr_PrintEx(1)


#define _m_d_getAttr_DECREF( obj, name, pvar, c_type_name, py_type, py_c_type )     \
  {																	\
	PyObject *_m_attr = PyObject_GetAttrString( obj, name);		    \
    if( _m_attr == NULL && PyErr_Occurred() != NULL)                \
      _m_printPyError_withLineNumber();                              \
    else if (_m_attr == NULL )                                      \
      fprintf(stderr,"WARNING: _m_getAttr_DECREF: No error occurred,\
 but the returned pointer was still NULL!\n"); \
    else                                                            \
      {                                                             \
        if( !Py##py_type##_Check( _m_attr ) )\
          fprintf(stderr,"WARNING: _m_getAttr_DECREF: The value returned by attribute '%s' was not the expected type!\n", name); \
        else                                                            \
          *(c_type_name *)(pvar) = Py##py_type##_AS_##py_c_type(_m_attr);                       \
        Py_DECREF(_m_attr);                                             \
      }                                                                 \
  }

#define _m_d_setAttr_DECREF( obj, name, function, arg )               \
  {																	\
	PyObject *val = function(arg);                                  \
    if( val == NULL && PyErr_Occurred() != NULL)                    \
      _m_printPyError_withLineNumber();                              \
    else if (val == NULL )                                          \
      fprintf(stderr,"WARNING: _m_setAttr_DECREF: No error occurred,\
 but the returned pointer was still NULL!\n"); \
    else                                                            \
      {                                                             \
        PyObject_SetAttrString( obj, name, val);                    \
        Py_DECREF(val);                                             \
      }                                                             \
  }

// Import/instantiate
#define newObject(mod, name) \
  _m_d_newObject( #mod, #name )


// Getters
// NOTE: these three use a different footprint for the _m_getAttr_DECREF, as they need to check
// the python object type.
#define getBoolAttr(obj, name, pvar) _m_d_getAttr_DECREF( obj, #name, pvar, bool, Int,LONG)
#define getLongAttr(obj, name, pvar) _m_d_getAttr_DECREF( obj, #name, pvar, long, Int, LONG)
#define getDoubleAttr(obj, name, pvar) _m_d_getAttr_DECREF( obj, #name, pvar, double, Float, DOUBLE )

// NOTE: caller is responsible for checking return values for strings!
#define getStringAttr(obj, name, pyo) ((char *)PyString_AS_STRING(pyo=PyObject_GetAttrString(obj, #name)))

// caller responsible for checking return values of lists.
#define getListAttr(obj, name) PyObject_GetAttrString(obj, #name)


#define pingAttr(obj, name) { \
  PyObject *_m_attr = PyObject_GetAttrString( obj, #name );\
  if (_m_attr == NULL && PyErr_Occurred() != NULL )        \
    _m_printPyError_withLineNumber();                       \
  else if (_m_attr == NULL )                               \
    fprintf(stderr,"WARNING: pingAttr: No error occurred,\
 but the returned pointer was still NULL!\n"); \
  else { Py_DECREF(_m_attr); }               \
  }
    

// List indexing
// NOTE: caller is responsible for checking return values for strings!

#define getStringItem(list, index) PyString_AS_STRING(PyList_GET_ITEM(list, index))

// The following works without ref counting issues as PyList_GET_ITEM borrows the references.
// Since these macros must be expressions and not statements, we
//cannot raise error conditions here easily, thus it is caller's
//responsibility.
#define getLongItem(list, index) \
  (PyInt_Check(PyList_GET_ITEM(list, index))?PyInt_AS_LONG(PyList_GET_ITEM(list,index)):-1)

#define getLongItemFromTuple(tuple, index) \
  (PyInt_Check(PyTuple_GET_ITEM(tuple, index))?PyInt_AS_LONG(PyTuple_GET_ITEM(tuple,index)):-1)

// Setters
#define setDoubleAttr(obj, name, arg) _m_d_setAttr_DECREF( obj, #name, PyFloat_FromDouble, arg)
#define setLongAttr(obj, name, arg) _m_d_setAttr_DECREF( obj, #name, PyLong_FromLong, arg)

// Testers
#define testLongAttr(obj, name, test, value) _m_d_testLongAttr( obj, #name, #test, value )
#define testBoolAttr(obj, name) _m_d_testLongAttr( obj, #name, "=", 1 )

// Function calls
// TODO: no debug versions of these yet. 
#define callFunc_NoArgsToNone(obj, name) PyObject_CallMethod(obj, #name, "()")
#define callFunc_NoArgsToLong(obj, name) PyInt_AS_LONG(PyObject_CallMethod(obj, #name, "()"))
#define callFunc_DoubleToNone(obj, name, arg) PyObject_CallMethod(obj, #name, "(f)", arg)

// Print calls, currently null statements.
#define m_printStatusLine( obj,a,b,c)
//#define m_printStatusLine( obj,a)
#define m_printTrajLine(obj,a,b)
#define m_printStatusLine_First_Bimolecular( obj,a,b,c,d,e)
#define m_printStatusLine_Final_First_Bimolecular( obj, a,b,c,d,e,f )
#define m_printStatusLine_Warning( obj, a, b )

// Not currently used, but might be a good reference for later
#define callFunc_IntToNone(obj, name, arg) PyObject_CallObject(PyObject_GetAttrString(obj, #name), Py_BuildValue("(i)", arg))
// these are refcounting insensitive at the moment.
#define callFunc_IntToString(obj, name, arg) PyString_AS_STRING(PyObject_CallObject(PyObject_GetAttrString(obj, #name), Py_BuildValue("(i)", arg)))


#endif




/*****************************************************

Inline functions for various python macro operations.

Prefixed by _m_ if used internally by macros, and by 
_m_d_ if it's a debug version for when DEBUG_MACROS is set.

*****************************************************/

static inline bool _m_testLongAttr( PyObject *obj, const char *attrname, const char *test, long value )
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

#ifdef DEBUG_MACROS
static inline bool _m_d_testLongAttr( PyObject *obj, const char *attrname, const char *test, long value )
{
 PyObject *_m_attr = PyObject_GetAttrString( obj, attrname);                
 long local_val;
 if( _m_attr == NULL && PyErr_Occurred() != NULL )
   {
     _m_printPyError_withLineNumber();                           
     return false;
   }
 else if (_m_attr == NULL )
   {
     fprintf(stderr,"WARNING: _m_d_testLongAttr: No error occurred,\
 but the returned object from GetAttrString was still NULL!\n");
     return false;
   }
 else
   {
     if( !PyInt_Check( _m_attr ) )
       {
         fprintf(stderr,"ERROR: _m_d_testLongAttr: Attribute name %s was not an integer type or subclass.\n", attrname );
         Py_DECREF(_m_attr);
         return false;
       }
     local_val = PyInt_AS_LONG(_m_attr);
     Py_DECREF(_m_attr);                                                                                        
     if( test[0] == '=' )
       return (local_val == value);
     if( test[0] == '<' )
       return (local_val < value );
     if( test[0] == '>' )
       return (local_val > value );
   }
}
#endif

static inline PyObject *_m_newObject( const char *mod, const char *name )
{
  PyObject *module = NULL;
  PyObject *class_obj = NULL;
  PyObject *new_obj = NULL;

  module = PyImport_ImportModule( mod ); // new reference
  if(module == NULL )
    return NULL;

  class_obj = PyObject_GetAttrString(module, name ); // new reference
  if (class_obj == NULL )
    {
      Py_DECREF(module);
      return NULL;
    }

  new_obj = PyObject_CallObject( class_obj, NULL );

  Py_DECREF( module );
  Py_DECREF( class_obj );
  // new_obj is a new reference, which we return. Caller is responsible.
  return new_obj;
}


#ifdef DEBUG_MACROS
static inline PyObject *_m_d_newObject( const char *mod, const char *name )
{
  PyObject *module = NULL;
  PyObject *class_obj = NULL;
  PyObject *new_obj = NULL;

  module = PyImport_ImportModule( mod ); // new reference
  if (module == NULL && PyErr_Occurred() != NULL)
    _m_printPyError_withLineNumber();                            
  else if (module == NULL)
    {
      fprintf(stderr,"WARNING: _m_d_newObject: No error occurred,\
 but the returned module was still NULL!\n");
      return NULL;
    }
  class_obj = PyObject_GetAttrString(module, name ); // new reference
  if (class_obj == NULL && PyErr_Occurred() != NULL)
    _m_printPyError_withLineNumber();                                
  else if (class_obj == NULL)
    {
      fprintf(stderr,"WARNING: _m_d_newObject: No error occurred,\
 but the returned class object was still NULL!\n");
      return NULL;
    }

  new_obj = PyObject_CallObject( class_obj, NULL );
  if (new_obj == NULL && PyErr_Occurred() != NULL)
    _m_printPyError_withLineNumber();                                
  else if (new_obj == NULL)
    {
      fprintf(stderr,"WARNING: _m_d_newObject: No error occurred,\
 but the returned class object was still NULL!\n");
      return NULL;
    }
  Py_DECREF( module );
  Py_DECREF( class_obj );
  // new_obj is a new reference, which we return. Caller is responsible.
  return new_obj;
}

#endif // DEBUG_MACROS was set

/******************************************************

 Functions defined in python_options.cc

******************************************************/

// Functions
class identlist *makeID_list(PyObject *strand_list);
class stopcomplexes *getStopComplexList(PyObject *options, int index);
class identlist *getID_list(PyObject *options, int index);



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

#define SIMULATION_MODE_ENERGY_ONLY         0x10

// simulation modes are bitwise -> bit 0 is normal/first bi
//                                 bit 1 is normal interface/python interface
//                                 bit 4 is compute energy mode only, should not 
//                                          be combined with any other flags.
// the following are the bit definitions for tests on those:

#define SIMULATION_MODE_FLAG_FIRST_BIMOLECULAR         0x01
#define SIMULATION_MODE_FLAG_PYTHON                    0x02

// stopconditions used in ssystem. 
// TODO: clean up/add docs.

#define STOPCONDITION_NORMAL           1
#define STOPCONDITION_REVERSE          2
#define STOPCONDITION_TIME            -1
#define STOPCONDITION_FORWARD          3
#define STOPCONDITION_ERROR           -2


#endif
