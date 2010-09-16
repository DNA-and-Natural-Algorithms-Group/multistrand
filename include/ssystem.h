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

class SimulationSystem
{
 public: 
  SimulationSystem( PyObject *system_options );

  // the following constructor is probably defunct now.
  SimulationSystem( int argc, char **argv );

  ~SimulationSystem( void ); 

  void StartSimulation( void );

 private:
  void StartSimulation_First_Bimolecular( void );
  void SimulationLoop( void );
  void SimulationLoop_First_Bimolecular( double *completiontime, int *completiontype, double *frate, char **tag );
  void InitializeSystem( void );

  void InitializeRNG( void );
  void generateNextRandom( void );

  EnergyModel *dnaEnergyModel;

  StrandComplex *startState; 
  SComplexList *complexList;

  PyObject *system_options;

  long current_seed;
  bool initial_trajectory;
  long simulation_mode;
  long simulation_count_remaining;
};

#endif
