/*
  Copyright (c) 2009-2010 Caltech. All rights reserved.
  Coded by: Joseph Schaeffer (schaeffer@dna.caltech.edu)

  A simple extension module for python that exposes the
  SimulationSystem object as a createable object that has one method.

*/

#include "Python.h"
#include "ssystem.h"
/* need C++ for ssystem.h... */
#include "options.h"
#include <string.h>
/* for strcmp */

typedef struct {
  PyObject_HEAD
  SimulationSystem *ob_system;  /* Our one data member, no other attributes. */
  //  PyObject  *member_attr;   
} SimSystemObject;

// static PyTypeObject SimSystem_Type; /* Forward decl */

#define SimSystem_Check(v)  (Py_TYPE(v) == &SimSystem_Type)
/* Should we want to use it later... */

static PyObject *SimSystemObject_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  SimSystemObject *self;

  self = (SimSystemObject *)type->tp_alloc(type,0);
  /* uses the tp_alloc to create the right amount of memory - this is done since we allow sub-classes. */
  self->ob_system = NULL;
  return (PyObject *)self;
}


static int SimSystemObject_init(SimSystemObject *self, PyObject *args)
{
  PyObject *options_object;
  if( !PyArg_ParseTuple(args, "O:SimSystem()", &options_object) )
    return -1;

  Py_INCREF( options_object );  /* Will be decreffed in dealloc, or
                                   here if there's a type error. */
  /* check the type */
  if( strcmp(options_object->ob_type->tp_name, "MultistrandOptions") != 0)
    {
      /* Note that we'll need to change the above once it's packaged nicely. */
      Py_DECREF(options_object);
      PyErr_SetString(PyExc_TypeError,
			        "Must be passed a single MultistrandOptions object.");
      return -1;
    }
  self->ob_system = new SimulationSystem( options_object );
  if( self->ob_system == NULL )  /* something horrible occurred */
    {
      Py_DECREF(options_object);
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
  /* Note that we no longer need the _threads variant, since it's not
     thread safe to do that at top level now. */

  Py_INCREF(Py_None);
  return Py_None;
}

static void SimSystemObject_dealloc( SimSystemObject *self )
{
  if( self->ob_system != NULL )
    {
      delete self->ob_system;
      self->ob_system = NULL;
    }
  self->ob_type->tp_free((PyObject *)self);
}


/*static PyObject *SimSystemObject_getattr( SimSystemObject *self, char *name )
{
  
}*/

static PyMethodDef SimSystemObject_methods[] = {
  {"start", (PyCFunction) SimSystemObject_start, METH_VARARGS,
              PyDoc_STR("Start this simulation running. Returns when simulation has been completed.")},
  { NULL, NULL }  /* Sentinel */
                  /* Note that the dealloc, etc methods are not
                     defined here, they're in the type object's
                     methods table, not the basic methods table. */
};


static PyTypeObject SimSystem_Type = {
  /* Note that the ob_type field cannot be initialized here. */
  PyVarObject_HEAD_INIT(NULL,0)
  "SimSystem",    /* tp_name */
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
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,             /* tp_flags */
  PyDoc_STR("Multistrand Simulation System object"),    /* tp_doc */
  0,                              /* tp_traverse */
  0,                              /* tp_clear */
  0,                              /* tp_richcompare */
  0,                              /* tp_weaklistoffset */
  0,                              /* tp_iter */
  0,                              /* tp_iternext */
  SimSystemObject_methods,        /* tp_methods */
  0,                              /* tp_members */
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
