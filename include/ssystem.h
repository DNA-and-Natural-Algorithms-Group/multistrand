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
  SimulationSystem( void  );

  // the following constructor is probably defunct now.
  SimulationSystem( int argc, char **argv );

  ~SimulationSystem( void ); 

  void StartSimulation( void );

  PyObject *calculateEnergy( PyObject *start_state, int typeflag );
  int getErrorFlag( void );

 private:
  void StartSimulation_Standard( void );
  void StartSimulation_FirstStep( void );
  void StartSimulation_Trajectory( void );
  void StartSimulation_Transition( void );

  void SimulationLoop_Standard( void );
  void SimulationLoop_FirstStep( void );
  void SimulationLoop_Trajectory( void );
  void SimulationLoop_Transition( int ointerval, double otime );

  void InitializeSystem( PyObject *alternate_start = NULL);

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
