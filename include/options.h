#ifndef __PYTHON_OPTIONS_H__
#define __PYTHON_OPTIONS_H__

#include <python2.7/Python.h>

// Macros for Python/C interface

/***************************************/
/* Helper functions / internal macros. */
/*                                     */
/* These are used in all build types.  */
/*                                     */
/***************************************/

/* Utility */
#define _m_prepStatusTuple( seed, com_type, time, tag )    \
  Py_BuildValue("(lids)", seed,(int) (com_type), time, tag )

#define _m_prepTrajTuple( tag, time )\
  Py_BuildValue("(sd)", tag, time)

#define _m_prepStatusFirstTuple( seed, com_type, com_time, frate, tag) \
  Py_BuildValue("(lidds)", seed, com_type, com_time, frate, tag )

#define _m_prepComplexStateTuple( seed, id, names, sequence, structure, energy ) \
  Py_BuildValue("(lisssd)", seed, id, names, sequence, structure, energy )
/* These four prep functions return a new reference via Py_BuildValue, error checking and reference counting is the caller's responsibility. */

/* Accessors (ref counting caller responsibility */
#define getStringAttr(obj, name, pyo) ((char *)PyString_AS_STRING(pyo=PyObject_GetAttrString(obj, #name)))
#define getListAttr(obj, name) PyObject_GetAttrString(obj, #name)

/* List indexing (ref counting caller responsibility */
#define getStringItem(list, index) PyString_AS_STRING(PyList_GET_ITEM(list, index))


/***************************************/
/* Helper functions / internal macros. */
/*                                     */
/* NON DEBUG ONLY.                     */
/*                                     */
/***************************************/
#ifndef DEBUG_MACROS

#define _m_getAttr_DECREF( obj, name, function, pvar, vartype )     \
  do {																	\
	PyObject *_m_attr = PyObject_GetAttrString( obj, name);		\
	*(vartype *)(pvar) = function(_m_attr);                         \
	Py_DECREF(_m_attr);												\
  } while(0)

#define _m_setAttr_DECREF( obj, name, function, arg )               \
  do {																	\
	PyObject *val = function(arg);                                  \
    PyObject_SetAttrString( obj, name, val);                       \
	Py_DECREF(val);                                                 \
  } while(0)

#define _m_setStringAttr( obj, name, arg )   \
  do {                                                 \
    PyObject *pyo_str = PyString_FromString( arg ); \
    PyObject_SetAttrString( obj, name, pyo_str );   \
    Py_XDECREF( pyo_str );                          \
  } while(0)

// Import/instantiate
#define newObject(mod, name) _m_newObject( #mod, #name )

// Accessors (no ref counting needed )
#define getBoolAttr(obj, name, pvar) _m_getAttr_DECREF( obj, #name, PyInt_AS_LONG, pvar, bool)
#define getLongAttr(obj, name, pvar) _m_getAttr_DECREF( obj, #name, PyInt_AS_LONG, pvar, long)
#define getDoubleAttr(obj, name, pvar) _m_getAttr_DECREF( obj, #name, PyFloat_AS_DOUBLE, pvar, double)

// Accessors (borrowed refs only)
#define getLongItem(list, index) PyInt_AS_LONG(PyList_GET_ITEM(list, index))
#define getLongItemFromTuple(tuple, index) PyInt_AS_LONG(PyTuple_GET_ITEM(tuple, index))

/* Procedure calling (no ref counts) */
#define pingAttr(obj, name) Py_DECREF(PyObject_GetAttrString( obj, #name ))
// CB: changed this from Py_XDECREF to Py_DECREF because it was accessing the
// attribute twice for each call to pingAttr (a problem for incrementors)
// note: does not do anything crazy on a NULL return from GetAttrString, but if that
// returned null it might be an error...

// Setters
#define setDoubleAttr(obj, name, arg) _m_setAttr_DECREF( obj, #name, PyFloat_FromDouble, (arg))
#define setLongAttr(obj, name, arg) _m_setAttr_DECREF( obj, #name, PyLong_FromLong, (arg))

#define setStringAttr(obj, name, arg) _m_setStringAttr( obj, #name, (arg) )

// Testers
#define testLongAttr(obj, name, test, value) _m_testLongAttr( obj, #name, #test, value )
#define testBoolAttr(obj, name) _m_testLongAttr( obj, #name, "=", 1 )

/* // Function calls */

#define _m_pushList( obj, a, b ) \
  do {                                          \
    PyObject *pyo = a;                          \
    PyObject_SetAttrString( obj, #b, pyo );     \
    Py_DECREF(pyo);                             \
  }while(0)

#define addResultLine_Energy( obj, energy )                    \
  _m_pushList( obj, PyFloat_FromDouble( energy ), add_result_energy)
/* Note: if energy fails to create via PyFloat_FromDouble, it'll be
   NULL we'll probably segfault. The only failure mode I can forsee is
   if you don't correctly pass a double somehow, so that's a compile
   issue rather than runtime.

    The ref is a new one, but we decref always in pushList, as we no
    longer need the ref and the SetAttrString should cause the owning
    object to have a good ref to it.*/


#define printStatusLine( obj, seed, com_type, time, tag )                   \
  _m_pushList( obj, _m_prepStatusTuple( seed, com_type, time,(char *)(tag) ), add_result_status_line)

#define printTrajLine( obj, name, time ) \
  _m_pushList( obj, _m_prepTrajTuple( (char *)(name), time ), print_traj_line )

#define printStatusLine_First_Bimolecular( obj,seed,com_type,com_time,frate,tag)  \
  _m_pushList( obj, _m_prepStatusFirstTuple( seed, com_type, com_time, frate, (char *)(tag)), add_result_status_line_firststep )

#define printComplexStateLine( obj, seed, id, names, sequence, structure, energy ) \
  _m_pushList( obj, _m_prepComplexStateTuple( seed, id, names, sequence, structure, energy ), add_complex_state_line )

#define pushTrajectoryComplex( obj, seed, id, names, sequence, structure, energy ) \
  _m_pushList( obj, _m_prepComplexStateTuple( seed, id, names, sequence, structure, energy ), add_trajectory_complex )

#define pushTrajectoryInfo( obj, time ) \
  setDoubleAttr( obj, add_trajectory_current_time, time )

// This macro DECREFs the passed obj once it's done with it.
#define pushTransitionInfo( options_obj, obj ) \
  _m_pushList( options_obj, obj, add_transition_info )

#endif  // DEBUG_MACROS is FALSE (not set).

/***************************************************

  Debug version of macros:

 **************************************************/
#ifdef DEBUG_MACROS

#define _m_printPyError_withLineNumber() \
  fprintf(stderr,"ERROR: Python Interpreter error: file %s, line %d.\nERROR (Python): ", __FILE__, (int)__LINE__), PyErr_PrintEx(1)


#define printPyError_withLineNumber()                                   \
  do {                                                                  \
    if (PyErr_Occurred() != NULL )                                      \
      {                                                                 \
        fprintf(stderr,"ERROR: Python Interpreter error: file %s, line %d.\nERROR (Python): ", __FILE__, (int)__LINE__); \
        PyErr_PrintEx(1);                                               \
      }                                                                 \
  } while(0)

#define _m_d_getAttr_DECREF( obj, name, pvar, c_type_name, py_type, py_c_type )     \
  do {                                                                     \
	PyObject *_m_attr = PyObject_GetAttrString( obj, name);		    \
    if( _m_attr == NULL && PyErr_Occurred() != NULL)                \
      _m_printPyError_withLineNumber();                              \
    else if (_m_attr == NULL )                                      \
      fprintf(stderr,"WARNING: _m_d_getAttr_DECREF: No error occurred,\
 but the returned pointer was still NULL!\n"); \
    else                                                            \
      {                                                             \
        if( !Py##py_type##_Check( _m_attr ) )\
          fprintf(stderr,"WARNING: _m_d_getAttr_DECREF: The value returned by attribute '%s' was not the expected type!\n", name); \
        else                                                            \
          *(c_type_name *)(pvar) = Py##py_type##_AS_##py_c_type(_m_attr);                       \
        Py_DECREF(_m_attr);                                             \
      }                                                                 \
  } while(0)

#define _m_d_setAttr_DECREF( obj, name, function, arg )               \
  do {                                                                   \
	PyObject *val = function(arg);                                  \
    if( val == NULL && PyErr_Occurred() != NULL)                    \
      _m_printPyError_withLineNumber();                              \
    else if (val == NULL )                                          \
      fprintf(stderr,"WARNING: _m_d_setAttr_DECREF: No error occurred,\
 but the returned pointer was still NULL!\n"); \
    else                                                            \
      {                                                             \
        PyObject_SetAttrString( obj, name, val);                    \
        Py_DECREF(val);                                             \
      }                                                             \
  } while(0)

#define _m_d_setStringAttr( obj, name, arg )          \
  do {                                                   \
    PyObject *pyo_str = PyString_FromString( arg );   \
    if (pyo_str == NULL && PyErr_Occurred() != NULL ) \
      _m_printPyError_withLineNumber();               \
    if (PyObject_SetAttrString( obj, name, pyo_str ) == -1 && PyErr_Occurred() != NULL) \
      _m_printPyError_withLineNumber();               \
    Py_XDECREF( pyo_str );                            \
  } while(0)


// Import/instantiate
#define newObject(mod, name) \
  _m_d_newObject( #mod, #name )


// Getters
// NOTE: these three use a different footprint for the _m_getAttr_DECREF, as they need to check
// the python object type.
#define getBoolAttr(obj, name, pvar) _m_d_getAttr_DECREF( obj, #name, pvar, bool, Int,LONG)
#define getLongAttr(obj, name, pvar) _m_d_getAttr_DECREF( obj, #name, pvar, long, Int, LONG)
#define getDoubleAttr(obj, name, pvar) _m_d_getAttr_DECREF( obj, #name, pvar, double, Float, DOUBLE )

#define pingAttr(obj, name) { \
  PyObject *_m_attr = PyObject_GetAttrString( obj, #name );\
  if (_m_attr == NULL && PyErr_Occurred() != NULL )        \
    _m_printPyError_withLineNumber();                       \
  else if (_m_attr == NULL )                               \
    fprintf(stderr,"WARNING: pingAttr: No error occurred,\
 but the returned pointer was still NULL!\n"); \
  else { Py_DECREF(_m_attr); }               \
  }


/*  The following works without ref counting issues as PyList_GET_ITEM
    borrows the references.  Since these macros must be expressions
    and not statements, we cannot raise error conditions here easily,
    thus it is caller's responsibility. */

#define getLongItem(list, index) \
  (PyInt_Check(PyList_GET_ITEM(list, index))?PyInt_AS_LONG(PyList_GET_ITEM(list,index)):-1)

#define getLongItemFromTuple(tuple, index) \
  (PyInt_Check(PyTuple_GET_ITEM(tuple, index))?PyInt_AS_LONG(PyTuple_GET_ITEM(tuple,index)):-1)

// Setters
#define setDoubleAttr(obj, name, arg) _m_d_setAttr_DECREF( obj, #name, PyFloat_FromDouble, (arg))
#define setLongAttr(obj, name, arg) _m_d_setAttr_DECREF( obj, #name, PyLong_FromLong, (arg))

#define setStringAttr(obj, name, arg) _m_d_setStringAttr( obj, #name, (arg ))
// note: any reference this creates is not the responsibility of the caller.

// Testers
#define testLongAttr(obj, name, test, value) _m_d_testLongAttr( obj, #name, #test, value )
#define testBoolAttr(obj, name) _m_d_testLongAttr( obj, #name, "=", 1 )

#define _m_d_pushList( obj, a, b )                                      \
  do {                                                                  \
    PyObject *pyo = a;                                                  \
    if( pyo == NULL  && PyErr_Occurred() != NULL)                       \
      _m_printPyError_withLineNumber();                                 \
    else                                                                \
      {                                                                 \
        if(PyObject_SetAttrString( obj, #b, pyo ) == -1 && PyErr_Occurred() != NULL) \
          _m_printPyError_withLineNumber();                             \
        Py_DECREF(pyo);                                                 \
      }                                                                 \
  }while(0)

#define addResultLine_Energy( obj, energy )                    \
  _m_d_pushList( obj, PyFloat_FromDouble( energy ), add_result_energy)
/* Note: if energy fails to create via PyFloat_FromDouble, it'll be NULL and the
    pushList error checking will catch it. The ref is a new one, but we decref always
    in pushList, as we no longer need the ref and the SetAttrString should cause the owning
    object to have a good ref to it.*/

#define printStatusLine( obj, seed, com_type,time, tag )                    \
  _m_d_pushList( obj, _m_prepStatusTuple( seed, com_type, time,(char *)(tag) ), add_result_status_line)

#define printTrajLine( obj, name, time ) \
  _m_d_pushList( obj, _m_prepTrajTuple( (char *)(name), time ), print_traj_line )

#define printStatusLine_First_Bimolecular( obj,seed,com_type,com_time,frate,tag)  \
  _m_d_pushList( obj, _m_prepStatusFirstTuple( seed, com_type, com_time, frate, (char *)(tag)), add_result_status_line_firststep )

#define printComplexStateLine( obj, seed, id, names, sequence, structure, energy ) \
  _m_d_pushList( obj, _m_prepComplexStateTuple( seed, id, names, sequence, structure, energy ), add_complex_state_line )

// This macro DECREFs the passed obj once it's done with it.
#define pushTransitionInfo( options_obj, obj ) \
  _m_d_pushList( options_obj, obj, add_transition_info )

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

 Functions defined in interface/options.cc

******************************************************/

// Functions
class identList *makeID_list(PyObject *strand_list);
class stopComplexes *getStopComplexList(PyObject *options, int index);
class identList *getID_list(PyObject *options, int index, PyObject *alternate_start = NULL);



/*****************************************************

 #defines for const values used in interface/options.py

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


#define SIMULATION_MODE_NORMAL              0x0010
#define SIMULATION_MODE_FIRST_BIMOLECULAR   0x0030

#define SIMULATION_MODE_PYTHON_NORMAL       0x0040
#define SIMULATION_MODE_PYTHON_FIRST_BI     0x0060

#define SIMULATION_MODE_ENERGY_ONLY         0x0200


// simulation modes are bitwise -> bit 5 is normal mode
//                                 bit 6 is first step mode
//                                 bit 7 is python interface
//                                 bit 10 is compute energy mode only, should not
//                                          be combined with any other flags.
// the following are the bit definitions for tests on those:

#define SIMULATION_MODE_FLAG_NORMAL                    0x0010
#define SIMULATION_MODE_FLAG_FIRST_BIMOLECULAR         0x0020
#define SIMULATION_MODE_FLAG_PYTHON                    0x0040
#define SIMULATION_MODE_FLAG_TRAJECTORY                0x0080
#define SIMULATION_MODE_FLAG_TRANSITION                0x0100

// stopconditions used in ssystem.
// TODO: clean up/add docs.

// normal sim mode stop result flags.
#define STOPRESULT_NORMAL           0x11
#define STOPRESULT_TIME             0x12

// first step mode stop result flags
#define STOPRESULT_FORWARD          0x21
#define STOPRESULT_FTIME            0x22
#define STOPRESULT_REVERSE          0x24

// error states
#define STOPRESULT_ERROR            0x81
#define STOPRESULT_NAN              0x82
#define STOPRESULT_NOMOVES          0x84


#endif
// #ifdef __PYTHON_OPTIONS_H__
