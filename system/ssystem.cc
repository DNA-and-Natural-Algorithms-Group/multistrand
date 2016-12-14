/*
 h Copyright (c) 2007-2010 Caltech. All rights reserved.
 Coded by: Joseph Schaeffer (schaeffer@dna.caltech.edu)
 */

#include "options.h"
#include "ssystem.h"
#include "simoptions.h"

#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <vector>
#include <iostream>

#ifdef PROFILING
#include "google/profiler.h"
#include "google/heap-profiler.h"
#endif

int noInitialMoves = 0;
int timeOut = 0;

SimulationSystem::SimulationSystem(PyObject *system_o) {

//	std::cout << "Initializing SimulationSystem w python object \n";

	system_options = system_o;
	simOptions = new PSimOptions(system_o);

	construct();
	initialPrint();

}

SimulationSystem::SimulationSystem(SimOptions* options) {

	system_options = NULL;
	simOptions = options;

	construct();
	initialPrint();

}

void SimulationSystem::construct(void) {

// We no longer need the below line; we are guaranteed that options
// will have a good reference for the lifetime of our object, as the
// controlling wrapper in multistrand_module.cc grabs the reference.

	simulation_mode = simOptions->getSimulationMode();
	simulation_count_remaining = simOptions->getSimulationCount();

	if (Loop::GetEnergyModel() == NULL) {
		energyModel = NULL;
		energyModel = new NupackEnergyModel(simOptions->getPythonSettings());
		Loop::SetEnergyModel(energyModel);
	} else {
		energyModel = Loop::GetEnergyModel();
	}

	startState = NULL;
	complexList = NULL;

	// move these to sim_settings
	exportStatesInterval = (simOptions->getOInterval() >= 0);
	exportStatesTime = (simOptions->getOTime() >= 0);

}

void SimulationSystem::initialPrint(void) {

//	if (energyModel->useArrhenius() && (current_seed == 777)) {
	if (energyModel->useArrhenius()) {

		energyModel->printPrecomputedArrRates();

	}

}

SimulationSystem::SimulationSystem(void) {

	std::cout<< "Initializing SimulationSystem \n";

	simulation_mode = -1;
	simulation_count_remaining = -1;

	if (Loop::GetEnergyModel() == NULL) {
		energyModel = NULL;
	} else {
		energyModel = Loop::GetEnergyModel();
	}

	system_options = NULL;
	simOptions = NULL;
	startState = NULL;
	complexList = NULL;
}

int SimulationSystem::isEnergymodelNull(void) {

	return (energyModel == NULL);

}

SimulationSystem::~SimulationSystem(void) {
	if (complexList != NULL)
		delete complexList;
	complexList = NULL;

// the remaining members are not our responsibility, we null them out
// just in case something thread-unsafe happens.

	energyModel = NULL;
	simOptions = NULL;
	startState = NULL;
}

void SimulationSystem::StartSimulation(void) {

	InitializeRNG();

	if (simulation_mode & SIMULATION_MODE_FLAG_FIRST_BIMOLECULAR) {
		StartSimulation_FirstStep();
	} else if (simulation_mode & SIMULATION_MODE_FLAG_TRAJECTORY) {
		StartSimulation_Trajectory();
	} else if (simulation_mode & SIMULATION_MODE_FLAG_TRANSITION) {
		StartSimulation_Transition();
	} else
		StartSimulation_Standard();

	finalizeSimulation();

}

void SimulationSystem::StartSimulation_FirstStep(void) {

	while (simulation_count_remaining > 0) {
		if (InitializeSystem() != 0)
			return;

		SimulationLoop_FirstStep();
		finalizeRun();

	}

}

void SimulationSystem::StartSimulation_Standard(void) {

	while (simulation_count_remaining > 0) {
		if (InitializeSystem() != 0)
			return;

		SimulationLoop_Standard();
		finalizeRun();

	}

}

void SimulationSystem::StartSimulation_Transition(void) {

	while (simulation_count_remaining > 0) {
		if (InitializeSystem() != 0) {
			return;
		}

		SimulationLoop_Transition();
		finalizeRun();

	}
}

void SimulationSystem::StartSimulation_Trajectory(void) {

	while (simulation_count_remaining > 0) {

		if (InitializeSystem() != 0) {

			cout << "system not initialized; returning \n";
			return;
		}

		SimulationLoop_Trajectory();
		finalizeRun();

	}
}

void SimulationSystem::finalizeRun(void) {

	simulation_count_remaining--;
	pingAttr(system_options, increment_trajectory_count);

	generateNextRandom();
}

void SimulationSystem::finalizeSimulation(void) {

	if (noInitialMoves > 0) {

		cout << "No initial moves for this first step simulation x" << noInitialMoves << "\n";
	}

	if (timeOut > 0) {

		cout << "time-out detected x" << timeOut << "\n";

	}

	// display size of statespace if used
	if (SimOptions::countStates) {

		cout << "Explored " << countMap.size() << " states. \n";

	}

	cout << flush;
}

void SimulationSystem::SimulationLoop_Standard(void) {

	double rchoice, rate, stime, ctime;
	rchoice = rate = stime = ctime = 0.0;

	bool checkresult = false;
	class stopComplexes *traverse = NULL, *first = NULL;

	double maxsimtime = simOptions->getMaxSimTime();
	long stopcount = simOptions->getStopCount();
	long stopoptions = simOptions->getStopOptions();

	complexList->initializeList();

	rate = complexList->getTotalFlux();

	do {

		rchoice = rate * drand48();
		stime += (log(1. / (1.0 - drand48())) / rate);

		// 1.0 - drand as drand returns in the [0.0, 1.0) range, we need a (0.0,1.0] range.
		// see notes below in First Step mode.

		if (stime < maxsimtime) {
			// Why check here? Because we want to report the final state
			// as the one we were in before transitioning past the maximum
			// time, rather than the one after it. This is not entirely
			// obvious - to see why, look at a max time simulation for
			// finding an equilibrium distribution. What state are we
			// likely to observe after the time expires? Let's look at the
			// previous state - the more stable the state, the less total
			// rate out and thus the expected time for transition is
			// higher. So in our most likely case for what state we'll see
			// /after/ the max sim time is hit is actually the (more)
			// unstable state it transitioned to!

			// FD: Mathematically it is also the correct thing to do,
			// FD: when we remember the memoryless property of the Markov chain

			(void) complexList->doBasicChoice(rchoice, stime);

			///Add the state to the hashmap counter
			this->countState(complexList);

			rate = complexList->getTotalFlux();

			if (stopoptions) {

				if (stopcount <= 0) {
					simOptions->stopResultError(current_seed);
					return;
				}

				checkresult = false;
				first = simOptions->getStopComplexes(0);
				checkresult = complexList->checkStopComplexList(first->citem);
				traverse = first;

				while (traverse->next != NULL && !checkresult) {
					traverse = traverse->next;
					checkresult = complexList->checkStopComplexList(traverse->citem);
				}
				// Note: we cannot delete first here if checkresult != 0,
				// as traverse->tag may be needed. It will get checked at
				// that point and deleted once traverse->tag is used,
				// later.
				if (!checkresult) {
					delete first;
				}
			}
		}
	} while (stime < maxsimtime && !checkresult);

	if (stime == NAN) {

		simOptions->stopResultNan(current_seed);

	} else if (checkresult) {

		dumpCurrentStateToPython();
		simOptions->stopResultNormal(current_seed, stime, traverse->tag);
		delete first;

	} else { // stime >= maxsimtime

		dumpCurrentStateToPython();
		simOptions->stopResultTime(current_seed, maxsimtime);

	}
}

void SimulationSystem::countState(SComplexList* complexList) {

	if (SimOptions::countStates) {

		string myComplex = complexList->toString();

		if (countMap.count(myComplex) == 0) {

			countMap[myComplex] = 1;

		} else {

			countMap[myComplex]++;

		}

	}

}

void SimulationSystem::SimulationLoop_Trajectory() {

	double rchoice, rate, stime, last_trajectory_time;
	rchoice = rate = 0.0;

	double maxsimtime = simOptions->getMaxSimTime();
	long stopcount = simOptions->getStopCount();
	long stopoptions = simOptions->getStopOptions();

	bool stopFlag = false;
	long current_state_count = 0;
	class stopComplexes *traverse = NULL, *first = NULL;

	complexList->initializeList();
	rate = complexList->getTotalFlux();

// We start at the beginning of time.
	stime = 0.0;

// The last time we gave the output state.
	last_trajectory_time = 0.0;

	if (stopoptions) {
		if (stopcount <= 0) {
			simOptions->stopResultError(current_seed);
			return;
		}
		first = simOptions->getStopComplexes(0);
	}

	// write the initial state:
	if (exportStatesInterval) {
		exportInterval(stime, current_state_count);
	}

	do {

		rchoice = rate * drand48();
		stime += (log(1. / (1.0 - drand48())) / rate);
		// 1.0 - drand as drand returns in the [0.0, 1.0) range, we need a (0.0,1.0] range.
		// see notes below in First Step mode.

		if (exportStatesTime) {
			exportTime(stime, &last_trajectory_time);
		}

		int ArrMoveType = complexList->doBasicChoice(rchoice, stime);
		rate = complexList->getTotalFlux();
		current_state_count += 1;

		if (exportStatesInterval) {
			exportInterval(stime, current_state_count, ArrMoveType);
		}

		if (stopoptions) {
			stopFlag = false;
			stopFlag = complexList->checkStopComplexList(first->citem);
			traverse = first;
			while (traverse->next != NULL && !stopFlag) {
				traverse = traverse->next;
				stopFlag = complexList->checkStopComplexList(traverse->citem);
			}
		}

	} while (stime < maxsimtime && !stopFlag);

	if (stime == NAN) {

		simOptions->stopResultNan(current_seed);

	} else if (stopFlag) {

		simOptions->stopResultNormal(current_seed, stime, traverse->tag);

	} else {

		simOptions->stopResultTime(current_seed, stime);

	}

	if (first != NULL) {
		delete first;
	}
}

void SimulationSystem::SimulationLoop_Transition(void) {

	double rchoice, rate, stime, ctime;

	rchoice = rate = stime = ctime = 0.0;

	double maxsimtime, otime;
	maxsimtime = otime = -1.0;

	bool checkresult = false;
	bool stopFlag = false;
	bool state_changed = false;
	long stopcount = 0;
	class stopComplexes *traverse = NULL, *first = NULL;

	long sMode = simOptions->getSimulationMode();
	long ointerval = simOptions->getOInterval();
	long stopoptions = simOptions->getStopOptions();
	stopcount = simOptions->getStopCount();
	maxsimtime = simOptions->getMaxSimTime();
	otime = simOptions->getOTime();

	if (stopcount <= 0 || !stopoptions) {
		// this simulation mode MUST have some stop conditions set.
		simOptions->stopResultError(current_seed);
		return;
	}

// figure out which stop entries should cause us to halt, update a bool vector to
// have true in the indices corresponding to which stop states are halting states.

	boolvector stop_entries;
	boolvector transition_states;
	stop_entries.resize(stopcount, false);
	transition_states.resize(stopcount, false);

	complexList->initializeList();

	first = simOptions->getStopComplexes(0);
	traverse = first;
	checkresult = false;
	for (int idx = 0; idx < stopcount; idx++) {
		if (strstr(traverse->tag, "stop:") == traverse->tag)
			stop_entries[idx] = true;

		checkresult = complexList->checkStopComplexList(traverse->citem);

		transition_states[idx] = checkresult;
		traverse = traverse->next;
	}
	delete first;
	sendTransitionStateVectorToPython(transition_states, stime);
// start

	rate = complexList->getTotalFlux();
	state_changed = false;
	stopFlag = false;
	do {

		rchoice = rate * drand48();

		stime += (log(1. / (1.0 - drand48())) / rate);
		// 1.0 - drand as drand returns in the [0.0, 1.0) range, we need a (0.0,1.0] range.
		// see notes below in First Step mode.

		if (stime < maxsimtime) {
			// See note in SimulationLoop_Standard

			complexList->doBasicChoice(rchoice, stime);
			rate = complexList->getTotalFlux();

			// check if our transition state membership vector has changed
			checkresult = false;
			first = simOptions->getStopComplexes(0);
			traverse = first;
			for (int idx = 0; idx < stopcount; idx++) {

				checkresult = complexList->checkStopComplexList(traverse->citem);

				if (checkresult && stop_entries[idx] == true) {
					// multiple stop states could suddenly be true, we add
					// a status line entry for the first one found.
					if (!stopFlag) {
						simOptions->stopResultNormal(current_seed, stime, traverse->tag);
					}

					stopFlag = true;
				}

				if (!state_changed && transition_states[idx] != checkresult) {
					state_changed = true;
				}

				transition_states[idx] = checkresult;
				traverse = traverse->next;
			}
			delete first; // we can do this now as we no longer need to
						  // save the stoplist until the loop exits, due
						  // to moving printStatusLine to immediately upon
						  // finding the stopping condition.
			if (state_changed) {
				sendTransitionStateVectorToPython(transition_states, stime);
				state_changed = false;
			}
		}
	} while (stime < maxsimtime && !stopFlag);

	if (stime == NAN) {

		simOptions->stopResultNan(current_seed);

	} else if (stopFlag) {

		dumpCurrentStateToPython();

	} else { // stime >= maxsimtime

		dumpCurrentStateToPython();
		simOptions->stopResultTime(current_seed, maxsimtime);

	}

}

void SimulationSystem::SimulationLoop_FirstStep(void) {
	double rchoice, rate, stime = 0.0;
	bool stopFlag = false;
	double last_trajectory_time = 0.0;

	class stopComplexes *traverse = NULL, *first = NULL;
	double frate = 0.0;

	double maxsimtime = simOptions->getMaxSimTime();
	long stopcount = simOptions->getStopCount();
	long stopoptions = simOptions->getStopOptions();
	long ointerval = simOptions->getOInterval();

	long current_state_count = 0;

	complexList->initializeList();

	rate = complexList->getJoinFlux();

// scomplexlist returns a 0.0 rate if there was a single complex in
// the system, and a -1.0 rate if there are exactly 0 join moves. So
// the 0.0 rate should probably be caught, though if you use a
// single complex system for a starting state it's probably
// deserved.

	if (rate == 0.0) { // no initial moves

		noInitialMoves++;

		simOptions->stopResultBimolecular("NoMoves", current_seed, 0.0, 0.0,
		NULL);
		return;
	}

	rchoice = rate * drand48();

	int ArrMoveType = complexList->doJoinChoice(rchoice);

	if (exportStatesInterval) {
		exportInterval(stime, current_state_count, ArrMoveType);
	}

// store the forward rate used for the initial step so we can record it.
	frate = rate * energyModel->getJoinRate_NoVolumeTerm() / energyModel->getJoinRate();

// rate is the total flux across all join moves - this is exactly equal to total_move_count *
// dnaEnergyModel->getJoinRate()

// This join rate is the dG_volume * bimolecular scaling constant
// used for forward transitions.  What we actually need is the
// bimolecular scaling constant * total move count. (dG volume is
// the volume dependent term that is not actually related to the
// 'collision' rate, but rather the volume we are simulating.

// Begin normal steps.
	rate = complexList->getTotalFlux();
	do {

		rchoice = rate * drand48();
		stime += (log(1. / (1.0 - drand48())) / rate);

		// trajectory output via outputtime option
		if (exportStatesTime) {
			exportTime(stime, &last_trajectory_time);
		}

		int ArrMoveType = complexList->doBasicChoice(rchoice, stime);
		rate = complexList->getTotalFlux();
		current_state_count++;

		if (exportStatesInterval) {
			exportInterval(stime, current_state_count, ArrMoveType);
		}

		if (stopcount > 0 && stopoptions) {
			stopFlag = false;
			first = simOptions->getStopComplexes(0);
			traverse = first;
			stopFlag = complexList->checkStopComplexList(traverse->citem);
			while (traverse->next != NULL && !stopFlag) {
				traverse = traverse->next;
				stopFlag = complexList->checkStopComplexList(traverse->citem);
			}
			if (!stopFlag && first != NULL)
				delete first;
		}

	} while (stime < maxsimtime && !stopFlag);

	if (stopFlag) {
		dumpCurrentStateToPython();
		if (strcmp(traverse->tag, "REVERSE") == 0)
			simOptions->stopResultBimolecular("Reverse", current_seed, stime, frate, traverse->tag);
		else
			simOptions->stopResultBimolecular("Forward", current_seed, stime, frate, traverse->tag);
		delete first;
	} else {
		timeOut++;
		dumpCurrentStateToPython();
		simOptions->stopResultBimolecular("FTime", current_seed, stime, frate,
		NULL);
	}

}

///////////////////////////////////////////////////////////
// void dumpCurrentStateToPython( void );				  //
// 													  //
// Helper function to send current state to python side. //
///////////////////////////////////////////////////////////
void SimulationSystem::dumpCurrentStateToPython(void) {
	int id;
	char *names, *sequence, *structure;
	double energy;
	SComplexListEntry *temp;
	temp = complexList->getFirst();
	while (temp != NULL) {
		temp->dumpComplexEntryToPython(&id, &names, &sequence, &structure, &energy);
		printComplexStateLine(simOptions->getPythonSettings(), current_seed, id, names, sequence, structure, energy);
		temp = temp->next;
	}
}

/////////////////////////////////////////////////////////////////////////////////////
// void sendTransitionStateVectorToPython( boolvector transition_states );		   //
// 																				   //
// Helper function to prepare a Python list object containing the bool information //
//  about which transition states we are in.									   //
/////////////////////////////////////////////////////////////////////////////////////

void SimulationSystem::sendTransitionStateVectorToPython(boolvector transition_states, double current_time) {
	PyObject *mylist = PyList_New((Py_ssize_t) transition_states.size());
// we now have a new reference here that we'll need to DECREF.

	if (mylist == NULL)
		return; // TODO: Perhaps raise an exception to the Python side here that
				// we couldn't pass the information back...

	boolvector_iterator it;
	Py_ssize_t index = 0;

	for (it = transition_states.begin(); it < transition_states.end(); it++) {
		if (*it) // bool value was true
		{
			Py_INCREF(Py_True);
			PyList_SET_ITEM(mylist, index, Py_True);
			// ownership of this reference to Py_True has now been stolen by PyList_SET_ITEM.
		} else {
			Py_INCREF(Py_False);
			PyList_SET_ITEM(mylist, index, Py_False);
			// ownership of this reference to Py_False has now been stolen by PyList_SET_ITEM.
		}
		index++;
	}

	PyObject *transition_tuple = Py_BuildValue("dO", current_time, mylist);
// we now have a new reference to transition_tuple.  note that our
// reference to mylist has NOT been stolen by this call [it was
// increffed in the call itself, so we can decref now with no
// worries]
	Py_DECREF(mylist);
// we now have no references to mylist directly owned [though transition tuple has one via BuildValue]

	pushTransitionInfo(system_options, transition_tuple);

// transition_tuple has been decreffed by this macro, so we no longer own any references to it

}

///////////////////////////////////////////////////////////
// void sendTrajectory_CurrentStateToPython( void );	  //
// 													  //
// Helper function to send current state to python side. //
///////////////////////////////////////////////////////////

void SimulationSystem::sendTrajectory_CurrentStateToPython(double current_time, int arrType) {
	int id;
	char *names, *sequence, *structure;
	double energy;
	SComplexListEntry *temp;
	temp = complexList->getFirst();
	while (temp != NULL) {

		temp->dumpComplexEntryToPython(&id, &names, &sequence, &structure, &energy);

		pushTrajectoryComplex(system_options, current_seed, id, names, sequence, structure, energy);

		temp = temp->next;
	}

	pushTrajectoryInfo(system_options, current_time);
	pushTrajectoryInfo2(system_options, arrType);

}

// FD: OK to have alternate_start = NULL
int SimulationSystem::InitializeSystem(PyObject *alternate_start) {
	class StrandComplex *tempcomplex;
	class identList *id;

	simOptions->generateComplexes(alternate_start, current_seed);

// FD: Somehow, check if complex list is pre-populated.
	startState = NULL;
	if (complexList != NULL)
		delete complexList;

	complexList = new SComplexList(energyModel);

// FD: this is the python - C interface
	for (unsigned int i = 0; i < simOptions->myComplexes->size(); i++) {

		char* tempSequence = copyToCharArray(simOptions->myComplexes->at(i).sequence);
		char* tempStructure = copyToCharArray(simOptions->myComplexes->at(i).structure);

		id = simOptions->myComplexes->at(i).list;

		tempcomplex = new StrandComplex(tempSequence, tempStructure, id);

		startState = tempcomplex;
		complexList->addComplex(tempcomplex);

	}

	return 0;
}

void SimulationSystem::InitializeRNG(void) {

	FILE *fp = NULL;

	if (simOptions->useFixedRandomSeed()) {
		current_seed = simOptions->getInitialSeed();
	} else {
		if ((fp = fopen("/dev/urandom", "r")) != NULL) { // if urandom exists, use it to provide a seed
			long deviceseed;
			(void) fread(&deviceseed, sizeof(long), 1, fp);

			current_seed = deviceseed;
			fclose(fp);
		} else // use the possibly flawed time as a seed.
		{
			current_seed = time(NULL);
		}
	}
// now initialize this generator using our random seed, so that we can reproduce as necessary.
	srand48(current_seed);
}

void SimulationSystem::generateNextRandom(void) {
	current_seed = lrand48();
	srand48(current_seed);
}

PyObject *SimulationSystem::calculateEnergy(PyObject *start_state, int typeflag) {
	double *values = NULL;
	PyObject *retval = NULL;

// calc based on current state, do not clean up anything.
	if (start_state != Py_None) {
		InitializeSystem(start_state);
		complexList->initializeList();
	}

	values = complexList->getEnergy(typeflag); // NUPACK energy output : bimolecular penalty, no Volume term.

	retval = PyTuple_New(complexList->getCount());
// New Reference, we return it.
// The complex list is a linked list and new items are added at the head; so we need to reverse the resulting list to get the data back out.
	for (int loop = complexList->getCount() - 1; loop >= 0; loop--)
		PyTuple_SET_ITEM(retval, loop, PyFloat_FromDouble(values[loop]));
// the reference from PyFloat_FromDouble is immediately stolen by PyTuple_SET_ITEM.

	delete[] values;

	return retval;
}

void SimulationSystem::printTransition(double input) {

	cout << "Using RNG =" << input;
	Move* myMove = startState->getChoice(&input);
	cout << ", we selected move: \n " << myMove->toString(simOptions->energyOptions) << " \n ";

}

void SimulationSystem::exportTime(double simTime, double* lastExportTime) {

	if (simTime - *lastExportTime > simOptions->getOTime()) {

		*lastExportTime += simOptions->getOTime();
		sendTrajectory_CurrentStateToPython(*lastExportTime);

	}

}

void SimulationSystem::exportInterval(double simTime, int transitionCount, int arrType) {

	if ((transitionCount % simOptions->getOInterval()) == 0) {

		sendTrajectory_CurrentStateToPython(simTime, arrType);

	}

}

void SimulationSystem::printAllMoves() {

	complexList->initializeList();

	// also generate the half contexts
	complexList->updateOpenInfo();

	complexList->printComplexList();

	startState->printAllMoves();

}

// FD: a simple peak into the initial state
void SimulationSystem::InitialInfo(void) {

	if (InitializeSystem() != 0) {
		return;
	}

	printAllMoves();

	// print info on bimolecular rates
	double biRate = complexList->getJoinFlux();

	cout << "\n";
	cout << "JoinFlux is " << biRate << "\n";

//	complexList->doJoinChoice(0.99999 * biRate);

	biRate = biRate / energyModel->applyPrefactors(energyModel->getJoinRate(), loopMove, loopMove);

	cout << "\n";
	cout << "#inner-nucleotide joins is " << biRate;
	cout << "\n";

//	for (int i = 0; i < MOVETYPE_SIZE; i++) {
//
//		for (int j = 0; j < MOVETYPE_SIZE; j++) {
//
//			int primeVal = moveutil::getPrimeCode((MoveType) i, (MoveType) j);
//
//			string desc = moveutil::primeToDesc(primeVal);
//
//			cout << primeVal << " : \" " << desc << " \n" << " \" ," << "\n";
//
//		}
//
//	}

}

