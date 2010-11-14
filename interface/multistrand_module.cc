/*
  Copyright (c) 2009-2010 Caltech. All rights reserved.
  Coded by: Joseph Schaeffer (schaeffer@dna.caltech.edu)

  A simple extension module for python that exposes the
  SimulationSystem object as a createable object that has one method.

*/

#include "Python.h"
#include "structmember.h"

#include "ssystem.h"
/* need C++ for ssystem.h... */
#include "options.h"
#include <string.h>
/* for strcmp */

#include "google/profiler.h"
#include "google/heap-profiler.h"

typedef struct {
  PyObject_HEAD
  SimulationSystem *ob_system;  /* Our one data member, no other attributes. */
  PyObject *options;   
} SimSystemObject;

#define SimSystem_Check(v)  (Py_TYPE(v) == &SimSystem_Type)
/* Should we want to use it later... */

static PyObject *SimSystemObject_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  SimSystemObject *self;

  self = (SimSystemObject *)type->tp_alloc(type,0);
  /* uses the tp_alloc to create the right amount of memory - this is done since we allow sub-classes. */
  if( self != NULL )
    {
      self->ob_system = NULL;
      self->options = NULL;  // will be later set in _init ...
    }
  return (PyObject *)self;
}


static int SimSystemObject_init(SimSystemObject *self, PyObject *args)
{
  if( !PyArg_ParseTuple(args, "O:SimSystem()", &self->options) )
    return -1;

  Py_INCREF( self->options );  /* Will be decreffed in dealloc, or
                                   here if there's a type error. */

  /* check the type */
  if( strcmp(self->options->ob_type->tp_name, "Options") != 0)
    {
      printf("[%s] options name\n",self->options->ob_type->tp_name);
      /* Note that we'll need to change the above once it's packaged nicely. */
      Py_DECREF(self->options);
      PyErr_SetString(PyExc_TypeError,
                      "Must be passed a single Options object.");
      return -1;
    }
  self->ob_system = new SimulationSystem( self->options );
  if( self->ob_system == NULL )  /* something horrible occurred */
    {
      Py_DECREF(self->options);
      PyErr_SetString(PyExc_MemoryError,
			        "Could not create the SimulationSystem [C++] object, possibly memory issues?.");
      return -1;
    }
  return 0;
}

static PyObject *SimSystemObject_start(SimSystemObject *self, PyObject *args)
{
  if( !PyArg_ParseTuple(args, ":start") )
    return NULL;
  
  if( self->ob_system == NULL )
    {
      PyErr_SetString( PyExc_AttributeError,
                       "The associated SimulationSystem [C++] object no longer exists, cannot start the system.");
      return NULL;
    }
  self->ob_system->StartSimulation();

  Py_INCREF(Py_None);
  return Py_None;
}

static int SimSystemObject_traverse( SimSystemObject *self, visitproc visit, void *arg )
{
  Py_VISIT(self->options);
  return 0;
}

static int SimSystemObject_clear( SimSystemObject *self )
{
  Py_CLEAR(self->options);
  return 0;
}

static void SimSystemObject_dealloc( SimSystemObject *self )
{
  if( self->ob_system != NULL )
    {
      delete self->ob_system;
      self->ob_system = NULL;
    }
  
  SimSystemObject_clear( self );

  self->ob_type->tp_free((PyObject *)self);
}

const char docstring_SimSystem[] = "\
Python Wrapper for Multistrand's C++ SimulationSystem object.\n\
\n\
Provides a very very simple interface to the StartSimulation method, to \n\
actually run the simulation. Otherwise fairly boring.\n";

const char docstring_SimSystem_start[] = "\
SimSystem.start( self )\n\
\n\
Start the simulation; only returns when the simulation has been completed. \n\
Information is only returned from the simulation via the Options object it \n\
was created with.\n";

const char docstring_SimSystem_init[] = "\
:meth:`multistrand.system.SimSystem.__init__( self, *args )`\n\
\n\
Initialization of SimSystem object:\n\
\n\
Arguments:\n\
options [type=:class:`multistrand.options.Options`]  -- The options to use for\n\
                                                        this simulation. Is a\n\
                                                        required argument.\n\
\n";

static PyMethodDef SimSystemObject_methods[] = {
  {"__init__", (PyCFunction) SimSystemObject_init, METH_COEXIST | METH_VARARGS,
   PyDoc_STR( docstring_SimSystem_init)},
  {"start", (PyCFunction) SimSystemObject_start, METH_VARARGS,
   PyDoc_STR( docstring_SimSystem_start)},
  { NULL, NULL }  /* Sentinel */
                  /* Note that the dealloc, etc methods are not
                     defined here, they're in the type object's
                     methods table, not the basic methods table. */
};

static PyMemberDef SimSystemObject_members[] = {
  {"options", T_OBJECT_EX, offsetof(SimSystemObject,options), 0,
   "The :class:`multistrand.options.Options` object controlling this simulation system."},
  {NULL} /* Sentinel */
};

static PyTypeObject SimSystem_Type = {
  /* Note that the ob_type field cannot be initialized here. */
  PyVarObject_HEAD_INIT(NULL,0)
  "multistrand.system.SimSystem",    /* tp_name */
  sizeof(SimSystemObject),    /* tp_basicsize */
  0,                          /* tp_itemsize  [it's something that's a relic, should be 0] */
  /* standard method defs are next */
  (destructor)SimSystemObject_dealloc, /* tp_dealloc */
  0,                             /* tp_print [not quite the same as str] */
  0,                             /* tp_getattr see comment out line below. */
  /* (getattrfunc)SimSystem_getattr, /\* tp_getattr - not sure I need this if I define only methods? *\/ */
  0,                              /* tp_setattr */
  0,                              /* tp_compare */
  0,                              /* tp_repr */
  0,                              /* tp_as_number */
  0,                              /* tp_as_sequence */
  0,                              /* tp_as_mapping */
  0,                              /* tp_hash */
  0,                              /* tp_call */
  0,                              /* tp_str */
  0,                              /* tp_getattro */
  0,                              /* tp_setattro */
  0,                              /* tp_as_buffer */
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC, /* tp_flags */
  PyDoc_STR( docstring_SimSystem ),        /* tp_doc */
  (traverseproc)SimSystemObject_traverse,  /* tp_traverse */
  (inquiry) SimSystemObject_clear,         /* tp_clear */
  0,                              /* tp_richcompare */
  0,                              /* tp_weaklistoffset */
  0,                              /* tp_iter */
  0,                              /* tp_iternext */
  SimSystemObject_methods,        /* tp_methods */
  SimSystemObject_members,        /* tp_members */
  0,                              /* tp_getset */
  0,                              /* tp_base */
  0,                              /* tp_dict */
  0,                              /* tp_descr_get */
  0,                              /* tp_descr_set */
  0,                              /* tp_dictoffset */
  (initproc)SimSystemObject_init, /* tp_init */
  0,                              /* tp_alloc */
  SimSystemObject_new,            /* tp_new */
  /* if we do want tp_new, it needs to be set in the module init.
   Or maybe not. Trying it here.*/
  /* 0,                              /\* tp_free *\/ */
  /* 0,                              /\* tp_is_gc *\/ */
};

static PyObject *System_initialize_energymodel( PyObject *self,PyObject *args )
{
  PyObject *options_object = NULL;

  if( !PyArg_ParseTuple(args, "|O:initialize_energy_model( [options])", &options_object) )
    return NULL;
  
  EnergyModel *temp = Loop::GetEnergyModel();

  if( temp != NULL )
    delete temp;

  if (options_object  == NULL || options_object == Py_None )
    Loop::SetEnergyModel( NULL );
  else
    {
      temp = NULL;
      if(  testLongAttr(options_object, parameter_type,=,0) )
        temp = new ViennaEnergyModel( options_object );
      else
        temp = new NupackEnergyModel( options_object );
      Loop::SetEnergyModel( temp );
    }
  Py_INCREF( Py_None );
  return Py_None;
}

static PyObject *System_calculate_energy( PyObject *self,PyObject *args )
{
  SimulationSystem *temp = NULL;
  PyObject *options_object = NULL;
  PyObject *start_state_object = NULL;
  PyObject *energy;
  int typeflag = 0;
  if( !PyArg_ParseTuple(args, "O|Oi:energy(state, options[, energytypeflag])", &start_state_object, &options_object, &typeflag) )
    return NULL;
  if( options_object != NULL )
    Py_INCREF( options_object );
  Py_INCREF( start_state_object );
  
  if( options_object != NULL )
    temp = new SimulationSystem( options_object );
  else
    {
      temp = new SimulationSystem();
      int err = temp->getErrorFlag();
      if( err != 0 )
        {
          PyErr_Format(PyExc_AttributeError,"No energy model available, cannot compute energy. Please pass an options object, or use multistrand.system.initialize_energy_model(...).\n");
          if( options_object != NULL )
            Py_XDECREF( options_object );
          Py_XDECREF( start_state_object );
          return NULL;
        }
    }
  energy = temp->calculateEnergy( start_state_object, typeflag );

  delete temp;

  Py_XDECREF( options_object );
  Py_XDECREF( start_state_object );
  return energy;
}

static PyObject *System_run_system( PyObject *self,PyObject *args )
{
#ifdef PROFILING
  HeapProfilerStart("ssystem_run_system.heap");
#endif

  SimulationSystem *temp = NULL;
  PyObject *options_object = NULL;
  int typeflag = 0;
  if( !PyArg_ParseTuple(args, "O|Oi:run_system( options )", &options_object) )
    return NULL;
  Py_INCREF( options_object );
  
  temp = new SimulationSystem( options_object );
  temp->StartSimulation();

  delete temp;
  temp = NULL;
  
#ifdef PROFILING
  HeapProfilerDump("run_system");
  HeapProfilerStop();
#endif

  Py_XDECREF( options_object );
  Py_INCREF( Py_None );
  return Py_None;

}


static PyMethodDef System_methods[] = {
  {"energy", (PyCFunction) System_calculate_energy, METH_VARARGS,
              PyDoc_STR(" \
energy( start_state, options = None, energy_type = 0)\n\
Computes the energy of the passed state [a list of complexes or resting states], using\
temperature, etc, settings from the options object passed.\n\n\
Parameters\n\
energy_type = 0 [default]: no volume terms included\n\
energy_type = 1: include dG_volume\n\
energy_type = 2: include dG_assoc [NUPACK equivalent]\n\
energy_type = 3: include dG_volume + dG_assoc\n\
\n\
options = None [default]: Use the already initialized energy model.\n\
options = ...: If not none, should be a multistrand.options.Options object, which will be used for initializing the energy model ONLY if there is not one already present.\n")},
  {"initialize_energy_model", (PyCFunction) System_initialize_energymodel, METH_VARARGS,
              PyDoc_STR(" \
initialize_energy_model( options = None )\n\
Initialize the Multistrand module's energy model using the options object given. If a model already exists, this will remove the old model and create a new one - useful for certain parameter changes, but should be avoided if possible. This function is NOT required to use other parts of the module - by default they will create the model if it's not found, or use the one already initialized; this adds control over exactly what model is being used.\n\n\
options [default=None]: when no options object is passed, this removes the old energy model and does not create a new one.\n")},
  {"run_system", (PyCFunction) System_run_system, METH_VARARGS,
              PyDoc_STR(" \
run_system( options )\n\
Run the system defined by the passed in Options object.\n")},
  {NULL}  /*Sentinel*/
};


PyMODINIT_FUNC
initsystem(void)
{
  PyObject *m;
  /* Finalize the simulation system object type */
  if ( PyType_Ready(&SimSystem_Type) < 0)
    return;

  m = Py_InitModule3("system", System_methods, "Base module for holding System objects.");
  if (m == NULL )
    return;

  Py_INCREF( &SimSystem_Type );
  PyModule_AddObject(m, "SimSystem", (PyObject *) &SimSystem_Type );

}
