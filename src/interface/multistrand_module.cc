/*
Multistrand nucleic acid kinetic simulator
Copyright (c) 2008-2023 California Institute of Technology. All rights reserved.
The Multistrand Team (help@multistrand.org)
*/

#include <iostream>
#include <time.h>

#include "multistrand_module.h"

#ifdef PROFILING
#include "google/profiler.h"
#include "google/heap-profiler.h"
#endif


/* Python type definition =================================================== */


static PyObject *SimSystemObject_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
    // Use `tp_alloc`, since we allow sub-classes.
    SimSystemObject *self = (SimSystemObject *) type->tp_alloc(type, 0);
    assert (self != NULL);
    self->sys = NULL;
    return (PyObject *) self;
}

static int SimSystemObject_init(SimSystemObject *self, PyObject *args) {
    PyObject *options = NULL;
    if (!PyArg_ParseTuple(args, "O:__init__", &options))
        return -1;
    if (PSimOptions::checkPythonType(options) != 0)
        return -1;

    // clear previous simulation state
    PSimOptions::clear(options);
    // initialize wrapped object
    self->sys = new SimulationSystem(options);
    // store reference if first instance
    SimSystemObject_store(self, options);
    return 0;
}

static void SimSystemObject_dealloc(SimSystemObject *self) {
    if (self->sys != NULL) {
        delete self->sys;
        self->sys = NULL;
    }
    Py_TYPE(self)->tp_free((PyObject *) self);
}


/* Python class utils ======================================================= */


// Read the reusable instance held by `Options`.
PyObject *SimSystemObject_lookup(PyObject* options) {
    return PyObject_GetAttrString(options, "_reusable_sim_system");
}

// Write the reusable instance held by `Options`.
void SimSystemObject_store(SimSystemObject *self, PyObject *options) {
    if (self->sys->isFirstInstance()) {
        // new value
        Py_INCREF(self);
        PyObject_SetAttrString(options, "_reusable_sim_system", (PyObject *) self);
        // old value
        Py_DECREF(Py_None);
    }
}

/*
** Instantiate a new `SimulationSystem`.
** - If this is the first instance, then wrap it in `SimSystemObject` and store
**   it inside `Options`.
** - If not, then reuse the `EnergyModel` from the previous instance.
 */
SimulationSystem *SimSystemObject_construct_sys(PyObject *options) {
    if (SimSystemObject_lookup(options) == Py_None) {
        SimSystemObject *ssystem = (SimSystemObject *) SimSystem_Type.tp_new(
            &SimSystem_Type, NULL, NULL);
        PyObject *args = PyTuple_Pack(1, options);
        if (SimSystem_Type.tp_init((PyObject *) ssystem, args, NULL) < 0) {
            Py_XDECREF(ssystem);
            Py_DECREF(args);
            PyErr_Print();
            return NULL;
        }
        Py_DECREF(args);
        return ssystem->sys;
    } else {
        return new SimulationSystem(options);
    }
}


/* Python class methods ===================================================== */


static PyObject *SimSystemObject_start(SimSystemObject *self, PyObject *args) {
    if (!PyArg_ParseTuple(args, ":start"))
        return NULL;
    assert (self->sys != NULL);
#ifdef PROFILING
    HeapProfilerStart("ssystem_start.heap");
#endif
    self->sys->StartSimulation();
#ifdef PROFILING
    HeapProfilerDump("start_sim");
    HeapProfilerStop();
#endif
    Py_RETURN_NONE;
}

static PyObject *SimSystemObject_initialInfo(SimSystemObject *self, PyObject *args) {
    if (!PyArg_ParseTuple(args, ":initialInfo"))
        return NULL;
    assert (self->sys != NULL);
    self->sys->initialInfo();
    Py_RETURN_NONE;
}

static PyObject *SimSystemObject_localTransitions(SimSystemObject *self, PyObject *args) {
    if (!PyArg_ParseTuple(args, ":localTransitions"))
        return NULL;
    assert (self->sys != NULL);
    self->sys->localTransitions();
    Py_RETURN_NONE;
}


/* Python module ============================================================ */


static PyObject *System_calculate_energy(PyObject *self, PyObject *args) {
    PyObject *options = NULL, *start_state_object = NULL;
    int energy_type = 0;

    if (!PyArg_ParseTuple(
            args, "OO|i:calculate_energy(state, options[, energy_type=0])",
            &start_state_object, &options, &energy_type))
        return NULL;
    if (PSimOptions::checkPythonType(options) != 0)
        return NULL;
    if (!(0 <= energy_type && energy_type <= 3)) {
        PyErr_SetString(PyExc_TypeError, "Invalid 'energy_type' argument!");
        return NULL;
    }

    SimulationSystem *ssystem = SimSystemObject_construct_sys(options);
    PyObject *energy = ssystem->calculateEnergy(start_state_object, energy_type);

    // free, unless owned by Python
    if (!ssystem->isFirstInstance())
        delete ssystem;
    return energy;
}

static PyObject *System_calculate_rate(PyObject *self, PyObject *args, PyObject *keywds) {
    PyObject *options = NULL;
    double start_energy, end_energy;
    int transition_type = 0;
    static char *kwlist[] = {
        "start_energy", "end_energy", "options", "transition_type", NULL };
    if (!PyArg_ParseTupleAndKeywords(
            args, keywds,
            "ddO|i:calculate_rate(start_energy, end_energy, options[, transition_type=0])",
            kwlist, &start_energy, &end_energy, &options, &transition_type))
        return NULL;
    if (PSimOptions::checkPythonType(options) != 0)
        return NULL;
    if (!(0 <= transition_type && transition_type <= 2)) {
        PyErr_SetString(PyExc_TypeError, "Invalid 'transition_type' argument!");
        return NULL;
    }

    SimulationSystem *ssystem = SimSystemObject_construct_sys(options);
    EnergyModel *energyModel = ssystem->GetEnergyModel();

    double drate = -1.0;
    if (transition_type == 0) // unimol
        drate = energyModel->returnRate(start_energy, end_energy, 0);
    else if (transition_type == 1) // bimol-join
        drate = energyModel->getJoinRate();
    else if (transition_type == 2) // bimol-break
        drate = energyModel->returnRate(start_energy, end_energy, 3);

    // free, unless owned by Python
    if (!ssystem->isFirstInstance())
        delete ssystem;
    return PyFloat_FromDouble(drate);
}

PyMODINIT_FUNC PyInit_system(void) {
    // create `multistrand.system` module
    PyObject *m = PyModule_Create(&moduledef);
    if (m == NULL) {
        PyErr_Print();
        return m;
    }
    // populate module with `SimSystem_Type`
    Py_INCREF(&SimSystem_Type);
    if (PyModule_AddType(m, &SimSystem_Type) < 0) {
        Py_DECREF(&SimSystem_Type);
        Py_DECREF(m);
        PyErr_Print();
        return NULL;
    }
    return m;
}
