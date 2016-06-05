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
#include <vector>

typedef std::vector<bool> boolvector;
typedef std::vector<bool>::iterator boolvector_iterator;

class SimulationSystem
{
 public: 
  SimulationSystem( PyObject *system_options );
  SimulationSystem( void  );

  // the following constructor is probably defunct now.
  SimulationSystem( int argc, char **argv );

  ~SimulationSystem( void ); 

  void StartSimulation( void );
  void InfoInitial ( void );

  PyObject *calculateEnergy( PyObject *start_state, int typeflag );
  int getErrorFlag( void );

  int InitializeSystem( PyObject *alternate_start = NULL);

 private:
  void StartSimulation_Standard( void );
  void StartSimulation_FirstStep( void );
  void StartSimulation_Trajectory( void );
  void StartSimulation_Transition( void );

  void SimulationLoop_Standard( void );
  void SimulationLoop_FirstStep( void );
  void SimulationLoop_Trajectory( long ointerval, double otime );
  void SimulationLoop_Transition( void );



  void InitializeRNG( void );
  void generateNextRandom( void );

  // helper function for sending current state to Python side
  void dumpCurrentStateToPython( void );
  void sendTrajectory_CurrentStateToPython( double current_time );
  void sendTransitionStateVectorToPython( boolvector transition_states, double current_time);


  EnergyModel *dnaEnergyModel;

  StrandComplex *startState; 
  SComplexList *complexList;

  PyObject *system_options;
  SimOptions *system_options_wrapper;

  long current_seed;
  bool initial_trajectory;
  long simulation_mode;
  long simulation_count_remaining;
};

#endif
