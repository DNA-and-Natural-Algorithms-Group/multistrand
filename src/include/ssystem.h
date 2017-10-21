/*
Copyright (c) 2017 California Institute of Technology. All rights reserved.
Multistrand nucleic acid kinetic simulator
help@multistrand.org
*/

/* SimulationSystem class header. This is the main object which controls the entire simulated system. */

#ifndef __SSYSTEM_H__
#define __SSYSTEM_H__

#include <vector>
#include <unordered_map>
#include <iostream>
#include <string>

#include "energymodel.h"
#include "scomplexlist.h"

typedef std::vector<bool> boolvector;
typedef std::vector<bool>::iterator boolvector_iterator;

// structs


// simply some simulation datapoints
struct SimTimer{

	SimTimer(SimOptions& myOptions);

	void advanceTime(void);
	bool checkForNewNucleotide(void);

public:
	double rchoice = 0.0;
	double rate = 0.0;
	double stime = 0.0;
	double maxsimtime = 0.0;
	double last_trajectory_time = 0.0;
	long stopcount = 0;
	long stopoptions = 0;

	const int initialCT = 5;
	int nuclAdded = 0;

private:

	// temporary flags:
	const bool cotranscriptional = false;
	const double delayCT = 0.001; // 1ms delay between adding nucleotides.



};



class SimulationSystem {
public:
	SimulationSystem(SimOptions* options);
	SimulationSystem(PyObject* system_options);
	SimulationSystem(void);

	// helper method for constructors
	void construct(void);

	~SimulationSystem(void);

	void StartSimulation(void);
	void InitialInfo(void);	// printing function
	void printTransition(double); // printing function

	PyObject *calculateEnergy(PyObject *start_state, int typeflag);
	int isEnergymodelNull(void);

private:
	void StartSimulation_Standard(void);
	void StartSimulation_FirstStep(void);
	void StartSimulation_Trajectory(void);
	void StartSimulation_Transition(void);

	void SimulationLoop_Standard(void);
	void SimulationLoop_FirstStep(void);
	void SimulationLoop_Trajectory(void);
	void SimulationLoop_Transition(void);

	int InitializeSystem(PyObject *alternate_start = NULL);

	void InitializeRNG(void);
	void generateNextRandom(void);
	void finalizeRun(void);
	void finalizeSimulation(void);

	// helper function for sending current state to Python side
	void dumpCurrentStateToPython(void);
	void sendTrajectory_CurrentStateToPython(double current_time, int arrType = -77);
	void sendTransitionStateVectorToPython(boolvector transition_states, double current_time);

	void exportTime(double& simTime, double& lastExportTime);
	void exportInterval(double simTime, int period, int arrType = -88);
	void exportTrajState(double simTime, double* lastExportTime, int period);

	void printAllMoves(void);

	EnergyModel* energyModel = NULL;
	StrandComplex *startState = NULL;
	SComplexList *complexList = NULL;
	SimOptions *simOptions = NULL;

	PyObject *system_options = NULL;

	long current_seed = NULL;
	long simulation_mode;
	long simulation_count_remaining;

	//bool triggers for output
	bool exportStatesTime = false;
	bool exportStatesInterval = false;

	// counters for timeouts and no-move initial states.
	int noInitialMoves = 0;
	int timeOut = 0;


};

#endif
