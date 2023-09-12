/*
Multistrand nucleic acid kinetic simulator
Copyright (c) 2008-2023 California Institute of Technology. All rights reserved.
The Multistrand Team (help@multistrand.org)
*/

#ifndef MULTISTRAND_MODULE_H_
#define MULTISTRAND_MODULE_H_

/*
 A simple C++ extension module that exposes the `SimulationSystem` (C++) class
 as a static wrapper type `SimSystem` (Python).
 */

#include "Python.h"
#include "structmember.h"

#include "options.h"
#include "simoptions.h"
#include "ssystem.h"

/* Python class ============================================================= */

// type definition
typedef struct {
    PyObject_HEAD
    SimulationSystem *sys;
} SimSystemObject;

PyObject *SimSystemObject_new(PyTypeObject*, PyObject*, PyObject*);
int       SimSystemObject_init(PyObject*, PyObject*, PyObject*);
void      SimSystemObject_dealloc(PyObject*);

// class utils
PyObject         *SimSystemObject_lookup(PyObject*);
void              SimSystemObject_store(SimSystemObject*, PyObject*);
SimulationSystem *SimSystemObject_construct_sys(PyObject*);

// class methods
PyObject *SimSystemObject_start(SimSystemObject*, PyObject*);
PyObject *SimSystemObject_initialInfo(SimSystemObject*, PyObject*);
PyObject *SimSystemObject_localTransitions(SimSystemObject*, PyObject*);

PyDoc_STRVAR(
  docstring_SimSystem,
  "Python wrapper for the Multistrand `SimulationSystem` class (C++).\n");
PyDoc_STRVAR(
  docstring_SimSystem_init,
  "__init__(options)\n--\n\n"
  "Initialize a `SimulationSystem` (C++) object and its `SimSystem` (Python)\n"
  "wrapper object with the given `options` object, which is used both for\n"
  "configuring the simulation and for storing its results.\n\n"
  "Parameters:\n"
  "- options:  `Options` object\n\n"
  "NOTE:\n"
  "If a `SimSystem` was previously created from the same `Options` object, then\n"
  "the new `SimSystem` will reuse the older instance's `EnergyModel` member (C++)\n"
  "to avoid reloading the parameter file from disk. Therefore, in cases where the\n"
  "`Options` need to change in between different `SimSystem` instances, the user\n"
  "needs to create a new `Options` object.\n");
PyDoc_STRVAR(
  docstring_SimSystem_start,
  "start()\n--\n\n"
  "Start the simulation and block until it has been completed. Information is only\n"
  "returned from the simulation via the `Options` object it was initialized with.\n");
PyDoc_STRVAR(
  docstring_SimSystem_initialInfo,
  "initialInfo()\n--\n\n"
  "Query information about the initial state of the simulation.\n");
PyDoc_STRVAR(
  docstring_SimSystem_localTransitions,
  "localTransitions()\n--\n\n"
  "Traverse all possible transitions given the initial state.\n");

static PyMemberDef SimSystemObject_members[] = {
  { NULL } /* Sentinel */
};
static PyMethodDef SimSystemObject_methods[] = {
  { "__init__", (PyCFunction) SimSystemObject_init,
    METH_COEXIST | METH_VARARGS, docstring_SimSystem_init },
  { "start", (PyCFunction) SimSystemObject_start,
    METH_VARARGS, docstring_SimSystem_start },
  { "initialInfo", (PyCFunction) SimSystemObject_initialInfo,
    METH_VARARGS, docstring_SimSystem_initialInfo },
  { "localTransitions", (PyCFunction) SimSystemObject_localTransitions,
    METH_VARARGS, docstring_SimSystem_localTransitions },
  { NULL, NULL, 0, NULL } /* Sentinel */
};
static PyTypeObject SimSystem_Type = {
  PyVarObject_HEAD_INIT(NULL, 0)
  .tp_name = "multistrand.system.SimSystem",
  .tp_basicsize = sizeof(SimSystemObject), .tp_itemsize = 0,
  /* standard method defs are next */
  .tp_dealloc = (destructor) SimSystemObject_dealloc,
  .tp_vectorcall_offset = 0, /* new slot semantics in Python 3.8 */
  .tp_getattr = 0, .tp_setattr = 0,
  .tp_as_async = 0, /* new slot semantics in Python 3.8 */
  .tp_repr = 0, .tp_as_number = 0, .tp_as_sequence = 0, .tp_as_mapping = 0,
  .tp_hash = 0, .tp_call = 0, .tp_str = 0, .tp_getattro = 0, .tp_setattro = 0,
  .tp_as_buffer = 0,
  .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
  .tp_doc = docstring_SimSystem,
  .tp_traverse = 0, .tp_clear = 0,
  .tp_richcompare = 0, .tp_weaklistoffset = 0, .tp_iter = 0, .tp_iternext = 0,
  .tp_methods = SimSystemObject_methods, .tp_members = SimSystemObject_members,
  .tp_getset = 0, .tp_base = 0, .tp_dict = 0,
  .tp_descr_get = 0, .tp_descr_set = 0, .tp_dictoffset = 0,
  .tp_init = (initproc) SimSystemObject_init, .tp_alloc = 0,
  .tp_new = (newfunc) SimSystemObject_new, .tp_free = 0, .tp_is_gc = 0,
};

/* Python module ============================================================ */

// module functions
PyObject *System_calculate_energy(PyObject*, PyObject*);
PyObject *System_calculate_rate(PyObject*, PyObject*, PyObject*);

PyDoc_STRVAR(
  docstring_calculate_energy,
  "calculate_energy(states, options, energy_type='EnergyType.loop')\n--\n\n"
  "Compute the free energy of the input states, initializing a `SimulationSystem`\n"
  "(C++) object as specified by `options`.\n\n"
  "Parameters:\n"
  "- states:       list of `Complex` objects\n"
  "- options:      `Options` object\n"
  "- energy_type:  one of\n"
  "    - EnergyType.loop:     [default] no volume or association terms included\n"
  "    - EnergyType.volume:   include dG_volume (no clear interpretation)\n"
  "    - EnergyType.complex:  include dG_assoc\n"
  "                           (NUPACK complex microstate energy without symmetry terms)\n"
  "    - EnergyType.tube:     include dG_volume + dG_assoc\n"
  "                           (system state energy summed over complexes)\n");
PyDoc_STRVAR(
  docstring_calculate_rate,
  "calculate_rate(start_energy, end_energy, options, transition_type='TransitionType.unimol')\n--\n\n"
  "Compute the transition rate for the current kinetic model, initializing a\n"
  "`SimulationSystem` (C++) object as specified by `options`.\n\n"
  "Parameters:\n"
  "- {start,end}_energy:  floats, *without* dG_assoc and dG_volume\n"
  "- options:             `Options` object\n"
  "- transition_type:     one of\n"
  "    - TransitionType.unimol:       [default] unimolecular transition\n"
  "    - TransitionType.bimol_join:   bimolecular join (energies are irrelevant)\n"
  "    - TransitionType.bimol_break:  bimolecular break (energies are relevant)\n\n"
  "NOTE:\n"
  "Currently, this legacy interface provides access only to the base rate,\n"
  "ignoring transition context features. This corresponds to the first argument\n"
  "of the `RateEnv` (C++) object.\n");

static PyMethodDef System_methods[] = {
  { "calculate_energy", (PyCFunction) System_calculate_energy,
    METH_VARARGS, docstring_calculate_energy },
  { "calculate_rate", (PyCFunction) System_calculate_rate,
    METH_VARARGS | METH_KEYWORDS, docstring_calculate_rate },
  { NULL, NULL, 0, NULL } /* Sentinel */
};
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    .m_name = "multistrand.system",
    .m_doc = "Python extension module wrapping the Multistrand C++ simulator.",
    .m_size = 0,
    .m_methods = System_methods,
};
PyMODINIT_FUNC PyInit_system(void);

#endif // MULTISTRAND_MODULE_H_
