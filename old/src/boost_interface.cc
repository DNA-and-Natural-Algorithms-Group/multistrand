/*
  Copyright (c) 2009-2010 Caltech. All rights reserved.
  Coded by: Chris Berlind (cberlind@dna.caltech.edu)
  
  Minor updates by:  Joseph Schaeffer (schaeffer@dna.caltech.edu)

  This file was originally in the KinD repository as multistrand.cc,
  and this version is the child of changeset [aa494e583f82] in that
  repository.

*/

/*
 * Boost.python wrapper for Multistrand.
 */
// For good luck:
#include <stdlib.h>

// Multistrand C++ source files. 
#include "ssystem.h"

// Wrapper
#include <boost/python.hpp>
using namespace boost::python;

BOOST_PYTHON_MODULE(multistrand)
{
  /* SimulationSystem class wrapper */

#ifndef PYTHON_THREADS

  // Choose the correct StartSimulation function
  void (SimulationSystem::*StartSimulation_void)(void) = &SimulationSystem::StartSimulation;
  
  class_<SimulationSystem>("SimSystem", init<PyObject*>())
    .def("start", StartSimulation_void)
  ;

#else

  class_<SimulationSystem>("SimSystem", init<PyObject*>())
    .def("start", &SimulationSystem::StartSimulation_threads)
  ;

#endif //PYTHON_THREADS

}
