/*
 Copyright (c) 2017 California Institute of Technology. All rights reserved.
 Multistrand nucleic acid kinetic simulator
 help@multistrand.org
 */

#include "options.h"
#include "ssystem.h"
#include "simoptions.h"
#include "statespace.h"

#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <vector>
#include <iostream>

SimulationSystem::SimulationSystem(PyObject *system_o) {

	system_options = system_o;
	simOptions = new PSimOptions(system_o);

	construct();
	energyModel->writeConstantsToFile();

}

SimulationSystem::SimulationSystem(SimOptions* options) {

	system_options = NULL;
	simOptions = options;

	construct();
	energyModel->writeConstantsToFile();

}

void SimulationSystem::construct(void) {

// We no longer need the below line; we are guaranteed that options
// will have a good reference for the lifetime of our object, as the
// controlling wrapper in multistrand_module.cc grabs the reference.

	simulation_mode = simOptions->getSimulationMode();
	simulation_count_remaining = simOptions->getSimulationCount();

	energyModel = new NupackEnergyModel(simOptions->getPythonSettings());
	Loop::SetEnergyModel(energyModel);

// move these to sim_settings
	exportStatesInterval = (simOptions->getOInterval() > 0);
	exportStatesTime = (simOptions->getOTime() >= 0);

	builder = Builder(simOptions);

}

SimulationSystem::SimulationSystem(void) {

	simulation_mode = -1;
	simulation_count_remaining = -1;

	if (Loop::GetEnergyModel() != NULL) {

		energyModel = Loop::GetEnergyModel();

	}

}

int SimulationSystem::isEnergymodelNull(void) {

	return (energyModel == NULL);

}

SimulationSystem::~SimulationSystem(void) {

	if (complexList != NULL) {
		delete complexList;
	}
	complexList = NULL;

// the remaining members are not our responsibility, we null them out
// just in case something thread-unsafe happens.

	if (energyModel != NULL) {
		delete energyModel;
	}

	if (simOptions->myComplexes != NULL) {
		delete simOptions->myComplexes;
	}

	if (simOptions != NULL) {
		delete simOptions;
	}

	if (complexList != NULL) {
		delete complexList;
	}

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

	// also ensure the builder does not remember the previous state
	builder.lastState = ExportData();

}

void SimulationSystem::finalizeSimulation(void) {

	if (noInitialMoves > 0 and simOptions->verbosity) {

		cout << "No initial moves x" << noInitialMoves << "   (set verbosity = 0 to suppress) \n";
	}

	if (timeOut > 0 and simOptions->verbosity) {

		cout << "time-out detected x" << timeOut << "   (set verbosity = 0 to suppress) \n";

	}

	if (simOptions->statespaceActive) {

		builder.writeToFile();

	}

}

void SimulationSystem::SimulationLoop_Standard(void) {

	SimTimer myTimer(*simOptions);
	stopComplexes *traverse = NULL, *first = NULL;

	bool checkresult = false;

	complexList->initializeList();
	myTimer.rate = complexList->getTotalFlux();

	do {

		myTimer.advanceTime();

		if (myTimer.stime < myTimer.maxsimtime) {
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

			(void) complexList->doBasicChoice(myTimer);

			myTimer.rate = complexList->getTotalFlux();

			if (myTimer.stopoptions) {

				if (myTimer.stopcount <= 0) {
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
	} while (myTimer.stime < myTimer.maxsimtime && !checkresult);

	if (myTimer.stime == NAN) {

		simOptions->stopResultNan(current_seed);

	} else if (checkresult) {

		dumpCurrentStateToPython();
		simOptions->stopResultNormal(current_seed, myTimer.stime, traverse->tag);
		delete first;

	} else { // stime >= maxsimtime

		dumpCurrentStateToPython();
		simOptions->stopResultTime(current_seed, myTimer.maxsimtime);

	}
}

void SimulationSystem::SimulationLoop_Trajectory() {

	SimTimer myTimer(*simOptions);
	stopComplexes *traverse = NULL, *first = NULL;

	bool stopFlag = false;
	long current_state_count = 0;

	complexList->initializeList();
	myTimer.rate = complexList->getTotalFlux();

	if (myTimer.stopoptions) {
		if (myTimer.stopcount <= 0) {
			simOptions->stopResultError(current_seed);
			return;
		}
		first = simOptions->getStopComplexes(0);
	}

	// write the initial state:
	if (exportStatesInterval) {
		exportInterval(myTimer.stime, current_state_count);
	}

	do {

		myTimer.advanceTime();

		if (debugTraces) {
			cout << "Printing my complexlist! *************************************** \n";
			cout << complexList->toString() << endl;
		}

		if (exportStatesTime) {
			exportTime(myTimer.stime, myTimer.last_trajectory_time);
		}

		double ArrMoveType = complexList->doBasicChoice(myTimer);
		myTimer.rate = complexList->getTotalFlux();
		current_state_count += 1;

		if (exportStatesInterval) {
			exportInterval(myTimer.stime, current_state_count, ArrMoveType);
		}

		if (myTimer.stopoptions) {

			stopFlag = false;
			stopFlag = complexList->checkStopComplexList(first->citem);
			traverse = first;
			while (traverse->next != NULL && !stopFlag) {
				traverse = traverse->next;
				stopFlag = complexList->checkStopComplexList(traverse->citem);
			}
		}

	} while (myTimer.stime < myTimer.maxsimtime && !stopFlag);

	if (myTimer.stime == NAN) {

		simOptions->stopResultNan(current_seed);

	} else if (stopFlag) {

		simOptions->stopResultNormal(current_seed, myTimer.stime, traverse->tag);
		// now export the tag to the builder as well
		builder.stopResultNormal(myTimer.stime, string(traverse->tag));

	} else {

		simOptions->stopResultTime(current_seed, myTimer.stime);

	}

	if (first != NULL) {
		delete first;
	}
}

void SimulationSystem::SimulationLoop_Transition(void) {

	SimTimer myTimer(*simOptions);
	stopComplexes *traverse = NULL, *first = NULL;

	bool checkresult = false;
	bool stopFlag = false;
	bool state_changed = false;

	if (myTimer.stopcount <= 0 || !myTimer.stopoptions) {
		// this simulation mode MUST have some stop conditions set.
		simOptions->stopResultError(current_seed);
		return;
	}

// figure out which stop entries should cause us to halt, update a bool vector to
// have true in the indices corresponding to which stop states are halting states.

	boolvector stop_entries;
	boolvector transition_states;
	stop_entries.resize(myTimer.stopcount, false);
	transition_states.resize(myTimer.stopcount, false);

	complexList->initializeList();

	first = simOptions->getStopComplexes(0);
	traverse = first;
	checkresult = false;

	for (int idx = 0; idx < myTimer.stopcount; idx++) {

		if (strstr(traverse->tag, "stop:") == traverse->tag)
			stop_entries[idx] = true;

		checkresult = complexList->checkStopComplexList(traverse->citem);

		transition_states[idx] = checkresult;
		traverse = traverse->next;
	}
	delete first;
	sendTransitionStateVectorToPython(transition_states, myTimer.stime);
// start

	myTimer.rate = complexList->getTotalFlux();
	state_changed = false;
	stopFlag = false;
	do {

		myTimer.advanceTime();

		if (myTimer.stime < myTimer.maxsimtime) {
			// See note in SimulationLoop_Standard

			complexList->doBasicChoice(myTimer);
			myTimer.rate = complexList->getTotalFlux();

			// check if our transition state membership vector has changed
			first = simOptions->getStopComplexes(0);
			checkresult = false;
			traverse = first;

			for (int idx = 0; idx < myTimer.stopcount; idx++) {

				checkresult = complexList->checkStopComplexList(traverse->citem);

				if (checkresult && stop_entries[idx]) {
					// multiple stop states could suddenly be true, we add
					// a status line entry for the first one found.
					if (!stopFlag) {
						simOptions->stopResultNormal(current_seed, myTimer.stime, traverse->tag);
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
				sendTransitionStateVectorToPython(transition_states, myTimer.stime);
				state_changed = false;
			}
		}
	} while (myTimer.stime < myTimer.maxsimtime && !stopFlag);

	if (myTimer.stime == NAN) {

		simOptions->stopResultNan(current_seed);

	} else if (stopFlag) {

		dumpCurrentStateToPython();

	} else { // stime >= maxsimtime

		dumpCurrentStateToPython();
		simOptions->stopResultTime(current_seed, myTimer.maxsimtime);

	}

}

void SimulationSystem::SimulationLoop_FirstStep(void) {

	SimTimer myTimer(*simOptions);
	stopComplexes *traverse = NULL, *first = NULL;

	bool stopFlag = false;

	double frate = 0.0;

	long ointerval = simOptions->getOInterval();
	long current_state_count = 0;

	complexList->initializeList();

	myTimer.rate = complexList->getJoinFlux();

	// if the toggle is set, export the initial state with arrType equal to flux
	// and timestamp -1
	if (simOptions->getPrintIntialFirstStep()) {
		exportInterval(-1.0, current_state_count, myTimer.rate);
	}

// scomplexlist returns a 0.0 rate if there was a single complex in
// the system, and a -1.0 rate if there are exactly 0 join moves. So
// the 0.0 rate should probably be caught, though if you use a
// single complex system for a starting state it's probably
// deserved.

	if (myTimer.rate == 0.0) { // no initial moves

		noInitialMoves++;

		simOptions->stopResultFirstStep(current_seed, 0.0, 0.0, result_type::STR_NOINITIAL.c_str());
		return;
	}

	myTimer.advanceTime(); // select an rchoice
	myTimer.stime = 0.0; // but reset the jump in time, because first step mode.

	int ArrMoveType = complexList->doJoinChoice(myTimer);

	if (exportStatesInterval) {
		exportInterval(myTimer.stime, current_state_count, ArrMoveType);
	}

// store the forward rate used for the initial step so we can record it.
	frate = myTimer.rate * energyModel->getJoinRate_NoVolumeTerm() / energyModel->getJoinRate();

// This join rate is the dG_volume * bimolecular scaling constant
// used for forward transitions.  What we actually need is the
// bimolecular scaling constant * total move count. (dG volume is
// the volume dependent term that is not actually related to the
// 'collision' rate, but rather the volume we are simulating.

// Begin normal steps.
	myTimer.rate = complexList->getTotalFlux();

	do {

		myTimer.advanceTime();

		if (debugTraces) {
			cout << "Printing my complexlist! *************************************** \n";
			cout << complexList->toString() << endl;
		}

		// trajectory output via outputtime option
		if (exportStatesTime) {
			exportTime(myTimer.stime, myTimer.last_trajectory_time);
		}

		int ArrMoveType = complexList->doBasicChoice(myTimer);

		myTimer.rate = complexList->getTotalFlux();
		current_state_count++;

		if (exportStatesInterval) {
			exportInterval(myTimer.stime, current_state_count, ArrMoveType);
		}

		if (myTimer.stopcount > 0 && myTimer.stopoptions) {

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

	} while (myTimer.stime < myTimer.maxsimtime && !stopFlag);

	if (stopFlag) {
		dumpCurrentStateToPython();
		simOptions->stopResultFirstStep(current_seed, myTimer.stime, frate, traverse->tag);
		delete first;
	} else {
		timeOut++;
		dumpCurrentStateToPython();
		simOptions->stopResultFirstStep(current_seed, myTimer.stime, frate, result_type::STR_TIMEOUT.c_str());
	}

}

///////////////////////////////////////////////////////////
// void dumpCurrentStateToPython( void );				  //
// 													  //
// Helper function to send current state to python side. //
///////////////////////////////////////////////////////////
void SimulationSystem::dumpCurrentStateToPython(void) {

	SComplexListEntry *temp = complexList->getFirst();
	ExportData data;

	while (temp != NULL) {

		temp->dumpComplexEntryToPython(data);
		printComplexStateLine(simOptions->getPythonSettings(), current_seed, data);

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

void SimulationSystem::sendTrajectory_CurrentStateToPython(double current_time, double arrType) {

	ExportData data;
	ExportData mergedData;

	SComplexListEntry *temp = complexList->getFirst();

	while (temp != NULL) {

		temp->dumpComplexEntryToPython(data);

		if (!simOptions->statespaceActive) {
			pushTrajectoryComplex(system_options, current_seed, data);
		}

		temp = temp->next;

		mergedData.merge(data);

	}

	// for now, keep exporting the state to the regular interface too.
	if (simOptions->statespaceActive) {

		builder.addState(mergedData, arrType);

	} else {

		pushTrajectoryInfo(system_options, current_time);
		pushTrajectoryInfo2(system_options, arrType);

	}

}

// FD: OK to have alternate_start = NULL
int SimulationSystem::InitializeSystem(PyObject *alternate_start) {

	StrandComplex *tempcomplex;
	identList *id;

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

	if (utility::debugTraces) {

		cout << "Done initializing!" << endl;

	}

	return 0;
}

void SimulationSystem::InitializeRNG(void) {

	FILE *fp = NULL;

	if (simOptions->useFixedRandomSeed()) {
		current_seed = simOptions->getSeed();
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

void SimulationSystem::exportTime(double& simTime, double& lastExportTime) {

	if (simTime - lastExportTime > simOptions->getOTime()) {

		lastExportTime += simOptions->getOTime();
		sendTrajectory_CurrentStateToPython(lastExportTime);

	}

}

void SimulationSystem::exportInterval(double simTime, int transitionCount, double arrType) {

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
void SimulationSystem::initialInfo(void) {

	if (InitializeSystem() != 0) {
		return;
	}

	printAllMoves();

	// print info on bimolecular rates
	double biRate = complexList->getJoinFlux();

	cout << "\n";
	cout << "JoinFlux is " << biRate << "\n";

	biRate = biRate / energyModel->applyPrefactors(energyModel->getJoinRate(), loopMove, loopMove);

	cout << "\n";
	cout << "#nucleotide joins is " << biRate << endl;
	cout << "joinrate is " << energyModel->applyPrefactors(energyModel->getJoinRate(), loopMove, loopMove) << " /s" << endl;
	cout << "join concentration is " << energyModel->simOptions->energyOptions->getJoinConcentration() << endl;

}

// FD: Build the statespace by taking N transitions, where N is the number of transitions out of the state
// If you call this function, the system will assume the active statespace toggle is set.
void SimulationSystem::localTransitions(void) {

	assert(simOptions->statespaceActive);
	energyModel->inspection = true;

	InitializeRNG(); // the output dir will be '0' if unset
	InitializeSystem();
	complexList->initializeList();
	complexList->updateOpenInfo();

	uint16_t N = complexList->getMoveCount();
	uint16_t collisions = round(complexList->getJoinFlux());

	for (uint16_t i = 0; i < (N + collisions); i++) {

		InitializeSystem();
		complexList->initializeList();
		complexList->updateOpenInfo();
		complexList->getTotalFlux();	 // required to set joinrate

		SimTimer myTimer(*simOptions);
		myTimer.rchoice = i + 0.01;

		// export the initial state
		if (exportStatesInterval) {
			exportInterval(myTimer.stime, 0, 8888);
		}

		// move the state using the i-th transition
		int ArrMoveType = complexList->doBasicChoice(myTimer);

		// export the state after the transition
		if (exportStatesInterval) {
			exportInterval(myTimer.stime, 1, ArrMoveType);
		}

		finalizeRun();

	}

	finalizeSimulation();

}

