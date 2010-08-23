/*
   Copyright (c) 2007-2008 Caltech. All rights reserved.
   Coded by: Joseph Schaeffer (schaeffer@dna.caltech.edu)
*/
 
/* SimulationSystem class header. This is the main object which controls the entire simulated system. */ 

#ifndef __SSYSTEM_H__
#define __SSYSTEM_H__

#include "energymodel.h"
#include "loop.h"
#include "scomplex.h"
#include "scomplexlist.h"

/* Defining PYTHON_THREADS will allow multithreading in Python to work properly when
   simulations are run in separate Python threads.
   
   Defining PYTHON_THREADS is only necessary when the Boost.Python interface is being
   used along with Python threading.
   
   Currently, this symbol only affects SimulationSystem::StartSimulation_threads(void).
*/
#ifdef PYTHON_THREADS
#include <boost/python.hpp>
#endif


class SimulationSystem
{
 public: 

  //  SimulationSystem( char *filename ); // standard constructor, filename is so we can load an energy model. Future constructors will handle command-line processing and sother input files.
  SimulationSystem( PyObject *system_options );

  //  SimulationSystem( Options &globalOpts );


#ifndef PYTHON_THREADS
	// we do not want this when using the python interface, as it requires the parsing code.
  SimulationSystem( int argc, char **argv );
#endif


  ~SimulationSystem( void ); 
  //  int LoadSystem( FILE *instream );
  //  int ResetSystem( void );
  //  int StartSimulation( int input_flags, int num_sims, double simtime );


  void StartSimulation( void );
#ifdef PYTHON_THREADS
  void StartSimulation_threads( void );
#endif


 private:
  void SimulateTrajectory( void );
  void StartSimulation_First_Bimolecular( void );
  void SimulationLoop( long r_seed );
  void SimulationLoop_First_Bimolecular( long r_seed, double *completiontime, int *completiontype, double *frate, char **tag );
  void InitializeSystem( void );

  EnergyModel *dnaEnergyModel;

  StrandComplex *startState; 
  SComplexList *complexList;

  PyObject *system_options;

  bool boltzmann_sampling;
  bool use_fixed_random_seed;
  bool initial_trajectory;
  long simulation_mode;
  long simulation_count_remaining;
};

#endif
